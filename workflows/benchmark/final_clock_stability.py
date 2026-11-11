"""Validate and integrate all existing group bootstrap records and stable levels."""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from final_clock_common import *

def validate_record(ctx,path,meta,tag,b):
    marker=path.with_suffix('.complete.json')
    if not source_valid(ctx,marker):raise RuntimeError('Missing/corrupt bootstrap cache: '+str(path))
    r=joblib.load(path);eps=ctx.cfg['stability']['coefficient_epsilon']
    ix=np.asarray(r['indices'],int);coef=np.asarray(r['coefficient'],float)
    if len(ix)!=len(coef) or len(np.unique(ix))!=len(ix) or not np.isfinite(coef).all():
        raise ValueError('Invalid bootstrap coefficients '+str(path))
    if r['bootstrap']!=b:raise ValueError('Wrong bootstrap identifier')
    draw=np.asarray(r['draw_groups']).astype(str)
    expected=np.random.default_rng(ctx.cfg['seed']+b).choice(meta.group_id.unique(),meta.group_id.nunique(),replace=True).astype(str)
    if not np.array_equal(draw,expected):raise ValueError('Bootstrap seed/group draw mismatch '+str(path))
    counts=Counter(draw);fit=Counter(map(str,r['fit_ids']))
    if fit!=Counter({str(row.sample_id):counts[str(row.group_id)] for row in meta.itertuples() if counts[str(row.group_id)]}):
        raise ValueError('Donor/sample multiplicity mismatch '+str(path))
    info=dict(lineage=ctx.lineage,scope=tag,bootstrap=b,seed=ctx.cfg['seed']+b,status='success',
        source='reused_checked_actual_bootstrap',path=str(path.resolve()),sha256=file_hash(path),
        n_fit_rows=len(r['fit_ids']),n_draw_groups=len(draw),n_unique_draw_groups=len(counts),
        n_selected=int((np.abs(coef)>eps).sum()),n_exposed=len(ix),alpha=r['params']['alpha'],l1_ratio=r['params']['l1_ratio'],
        solver=r.get('solver','source_record'),draw_group_counts_json=json.dumps(dict(counts),sort_keys=True),
        fit_ids_source='joblib record fit_ids; duplicates retained',independent_donor_verified=False)
    return info,ix,np.where(np.abs(coef)>eps,coef,0)

def run(ctx,smoke=False):
    ds=dataset(ctx);fullmeta=ds['meta'];p=len(ds['schema']['features']);B=ctx.options['n_bootstrap']
    if B!=ctx.cfg['stability']['n_bootstrap']:raise ValueError('Requested bootstrap count differs from reusable experiment')
    outputs=[];runs=[];countrows=[];pins=[]
    def save(rel,frame):
        target=ctx.path(rel);write_csv(target,frame);outputs.append(target)
    for tag in (['full_train'] if smoke else ['outer1','outer2','outer3','outer4','outer5','full_train']):
        base=ctx.path('stability/'+tag);definition=json.loads((base/'bootstrap_definition.json').read_text())
        ids=set(definition['fit_rows']);meta=fullmeta[fullmeta.sample_id.isin(ids)].copy()
        assert set(meta.sample_id)==ids and len(meta)==len(definition['fit_rows'])
        if tag!='full_train':
            info=ctx.json(f'models/K{ctx.K}/ElasticNet/fold{tag[-1]}/model.json')
            assert ids==set(info['train_ids']) and not ids&set(info['valid_ids'])
        paths=[base/'bootstrap'/f'boot_{b:04d}.joblib' for b in range(3 if smoke else B)]
        # Bootstrap fits are already complete: parallel I/O/hash validation uses
        # small records and shares the single read-only matrix, not 48 copies.
        with ThreadPoolExecutor(max_workers=ctx.options['io_workers']) as pool:
            checked=list(pool.map(lambda q:validate_record(ctx,q[1],meta,tag,q[0]),enumerate(paths)))
        if smoke:
            target=ctx.path('reports/final_clock_bootstrap_smoke.json')
            write_json(target,dict(status='passed',smoke_only=True,production_complete=False,scope=tag,
                checked_records=len(checked),B_requested=B,source='existing real bootstrap records; not newly fitted 3 bootstraps'))
            return [target]
        selected=np.zeros(p,np.int64);exposed=np.zeros(p,np.int64);positive=np.zeros(p,np.int64);negative=np.zeros(p,np.int64);total=np.zeros(p,float)
        for row,ix,coef in checked:
            assert ix.min()>=0 and ix.max()<p
            runs.append(row);pins.append(Path(row['path']));exposed[ix]+=1
            selected[ix]+=(coef!=0);positive[ix]+=(coef>0);negative[ix]+=(coef<0);total[ix]+=coef
        table=pd.read_csv(base/'bootstrap_all_features.csv',low_memory=False)
        assert table.feature.tolist()==ds['schema']['features']
        np.testing.assert_array_equal(selected,table.selected_count)
        np.testing.assert_array_equal(exposed,table.exposure_count)
        np.testing.assert_array_equal(positive,table.positive_count);np.testing.assert_array_equal(negative,table.negative_count)
        np.testing.assert_allclose(selected/B,table.selection_frequency,atol=1e-12,rtol=1e-12)
        np.testing.assert_allclose(total/B,table.coefficient_mean,atol=1e-7,rtol=1e-5)
        for th in [ctx.options['threshold'],ctx.options['sensitivity_threshold']]:
            s=table[table.selection_frequency>=th]
            countrows.append(dict(lineage=ctx.lineage,scope=tag,threshold=th,n_successful_bootstraps=B,
                n_stable_Celltype_Gene=len(s),n_stable_Celltypes=s.Celltype.nunique(),n_unique_stable_Genes=s.Gene.nunique()))
        LOG.info('%s %s: validated %d/%d actual bootstrap draws, coefficient counts and means',ctx.lineage,tag,len(checked),B)
    save('stability/bootstrap_runs_manifest.csv',pd.DataFrame(runs))
    coverage=ctx.csv('features/train/celltype_coverage.csv').set_index('Celltype')
    table['structural_subtype_absent_fraction']=1-table.Celltype.map(coverage.n_samples)/len(fullmeta)
    x=ds['X'].copy();x.eliminate_zeros()
    table['zero_expression_fraction']=1-np.asarray(x.getnnz(axis=0)).ravel()/len(fullmeta)
    table['exposure_fraction']=table.exposure_count/B
    save('stability/bootstrap_feature_selection_frequency.csv',table)
    stable=table[table.selection_frequency>=ctx.options['threshold']].copy()
    old=ctx.csv('stability/stable_Celltype_Gene.csv')
    assert set(stable.feature)==set(old.feature)
    save('stability/stable_celltype_gene_features.csv',stable)
    gene_rows=[];cell_rows=[]
    for name,q in stable.groupby('Celltype',sort=True):
        c=coverage.loc[name];eligible=table[table.Celltype==name]
        denominator=int(eligible.candidate.sum())
        cell_rows.append(dict(Celltype=name,n_stable_features=len(q),n_unique_stable_genes=q.Gene.nunique(),
            candidate_features=denominator,stable_fraction_of_candidates=len(q)/max(denominator,1),
            coefficient_sum=q.coefficient_mean.sum(),coefficient_absolute_sum=q.coefficient_mean.abs().sum(),
            coefficient_mean=q.coefficient_mean.mean(),selection_frequency_mean=q.selection_frequency.mean(),
            n_cells=int(c.n_cells),n_samples=int(c.n_samples),stable_per_1000_cells=1000*len(q)/max(c.n_cells,1),
            stable_genes=';'.join(sorted(set(q.Gene)))))
    cells=pd.DataFrame(cell_rows).sort_values(['stable_fraction_of_candidates','coefficient_absolute_sum'],ascending=False)
    for gene,q in stable.groupby('Gene',sort=True):
        pos=int((q.coefficient_mean>0).sum());neg=int((q.coefficient_mean<0).sum())
        gene_rows.append(dict(Gene=gene,n_stable_celltypes=q.Celltype.nunique(),Celltypes=';'.join(sorted(set(q.Celltype))),
            max_selection_frequency=q.selection_frequency.max(),mean_selection_frequency=q.selection_frequency.mean(),
            coefficient_sum=q.coefficient_mean.sum(),coefficient_absolute_sum=q.coefficient_mean.abs().sum(),
            positive_celltypes=pos,negative_celltypes=neg,direction='mixed' if pos and neg else 'positive_association' if pos else 'negative_association',
            direction_consistent_across_celltypes=not(pos and neg)))
    genes=pd.DataFrame(gene_rows)
    save('stability/stable_celltypes_ranked.csv',cells)
    save('stability/stable_genes_unique.csv',genes)
    save('stability/stable_genes_by_celltype.csv',stable[['Celltype','Gene','feature','selection_frequency','coefficient_mean','sign_consistency']])
    save('stability/stable_feature_counts_summary.csv',pd.DataFrame(countrows))
    sensitivity=pd.DataFrame(countrows)
    save('stability/bootstrap_threshold_sensitivity.csv',sensitivity)
    target=ctx.path('stability/final_clock_stability_definition.json')
    write_json(target,dict(status='complete',successful_bootstraps=B,scopes=6,total_checked_records=len(runs),
        n_new_bootstrap_fits=0,frequency_denominator='successful group bootstrap fits, verified all300 per scope',
        primary_threshold=ctx.options['threshold'],sensitivity_threshold=ctx.options['sensitivity_threshold'],
        coefficient_epsilon=ctx.cfg['stability']['coefficient_epsilon'],
        feature_mapping='explicit schema Celltype×Gene product, never splitting strings on underscore',
        K=ctx.K,preprocessing='variance/MI/scaler refit within each bootstrap; source inner-tuned fixed HP',
        full_stability_role='exploratory full-development explanation; not used to select outer validation features',
        group_unit='GSE/GSM proxy; true donor identity unavailable',
        missing_note='structural absence is inferred from observed subtype coverage; expression zeros are not all missingness'))
    outputs.append(target)
    # Keep source record hashes in the completion marker to validate real reuse.
    return outputs+pins
