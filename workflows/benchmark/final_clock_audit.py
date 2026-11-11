"""Actual file/model/cohort/figure audit before Phase3/4 integration."""
from concurrent.futures import ThreadPoolExecutor
import importlib.metadata,json,sys
from pathlib import Path
import joblib,numpy as np,pandas as pd
from final_clock_common import *

def run(ctx,smoke=False):
    outputs=[];records=[];models=[]
    def save(rel,frame):
        p=ctx.path(rel);write_csv(p,frame);outputs.append(p)
    for marker in ['features/train/complete.json','features/combined/complete.json']:
        if not source_valid(ctx,ctx.path(marker)):raise RuntimeError('Invalid source feature checkpoint '+marker)
    ds=dataset(ctx);train=ds['meta'];external=ctx.csv('features/combined/samples.csv')
    rank=ctx.csv('benchmark/model_comparison.csv');en=rank[rank.model=='ElasticNet'].sort_values(['R2_mean','MAE_mean','K'],ascending=[False,True,True])
    if int(en.iloc[0].K)!=ctx.K:raise RuntimeError('Cached stability K differs from best mean internal EN R2')
    obs=ctx.json('audit/obs_schema.json')[0]['obs_columns']
    donorlike=[c for c in obs if any(s in c.lower() for s in ['donor','patient','individual','subject','participant'])]
    if donorlike and ctx.cfg['fields'].get('donor') is None:
        raise RuntimeError('Potential donor identifiers require review: '+str(donorlike))
    rows=[];known_gsm=set(train.loc[~train.GSM.isin(['','nan','none']),'GSM'])
    inventory=ctx.json('audit/inventory.json')
    for item in inventory[1:]:
        m=ctx.csv('features/'+item['dataset']+'/samples.csv')
        for r in m.itertuples():
            overlap=r.sample_id in set(train.sample_id) or r.group_id in set(train.group_id) or (r.GSM in known_gsm)
            if overlap:raise RuntimeError('Development/external overlap '+r.sample_id)
            rows.append(dict(dataset=item['dataset'],sample_id=r.sample_id,group_id=r.group_id,GSM=r.GSM,
                overlaps_development=False,overlaps_combined=r.sample_id in set(external.sample_id) or r.GSM in set(external.GSM),
                independent_donor_verified=False,grouping='GSE/GSM proxy; no donor column available'))
    save('audit/final_clock_sample_overlap_audit.csv',pd.DataFrame(rows))
    for path in sorted(ctx.path('models').glob('K*/ElasticNet/fold*/model.json'))+sorted(ctx.path('stability').glob('outer*/threshold_*/model.json'))+[ctx.path('models/final_model/model.json')]:
        m=json.loads(path.read_text());est=path.parent/'estimator.joblib';pre=path.parent/'preprocessor.joblib'
        models.append(dict(model=m['model'],K=m['K'],fold=m.get('fold','full_development'),path=str(est.resolve()),
            preprocessor=str(pre.resolve()),n_features=m['effective_K'],n_train=len(m.get('train_ids',m.get('preprocessing_fit_ids',[]))),
            alpha=m.get('params',{}).get('alpha',''),l1_ratio=m.get('params',{}).get('l1_ratio',''),exists=est.exists(),
            role='historical_global_winner' if path.parent.name=='final_model' else 'stable_outer_train' if 'threshold_' in str(path) else 'benchmark_EN'))
    save('audit/existing_model_inventory.csv',pd.DataFrame(models))
    figure=[]
    out=ROOT/'Stromal_Combined_Figure/outputs'
    manifest=json.loads((out/'source_manifest_5.json').read_text())
    for p in ['Stromal_clock_A_F_ft30_all_5.png','Stromal_individual_aging_clock_Fibroblast_Endothelial_5.png']:
        for scope in ['internal','external','biological_interpretation']:
            figure.append(dict(figure=str(out/p),figure_sha256=file_hash(out/p),lineage=ctx.lineage,scope=scope,
                source=str(ctx.path('benchmark/cv_predictions.csv' if scope=='internal' else
                    'external/external_ft30_all_folds_predictions.csv' if scope=='external' else 'stability/stable_Celltype_Gene.csv')),
                model_identity='ElasticNet source fold5, same-fold internal-selected K; not final full-development' if scope!='biological_interpretation' else 'shared full-training bootstrap EN',
                external_scope='FT30 ALL100 includes adapt30; NOT final full-development no-adapt',
                selected_final_model=False))
    save('audit/figure_source_audit.csv',pd.DataFrame(figure))
    refs=[]
    for key,path in ctx.cfg['reference_paths'].items():
        p=Path(path)
        try:
            text=p.read_text()
            if p.suffix=='.ipynb':text='\n\n'.join(''.join(c['source']) for c in json.loads(text)['cells'] if c['cell_type']=='code')
            target=ctx.path(f'audit/final_clock_reference_{key}.source.txt');write_text(target,text);outputs.append(target)
            refs.append(dict(reference=key,path=str(p),status='read',sha256=file_hash(p),source_characters=len(text)))
        except Exception as e:refs.append(dict(reference=key,path=str(p),status='unavailable',reason=str(e)))
    save('audit/final_clock_references.csv',pd.DataFrame(refs))
    pins=[ctx.path(p) for p in ['features/train/pseudobulk.npz','features/train/samples.csv','features/train/schema.json',
        'features/combined/pseudobulk.npz','benchmark/cv_fold_metrics.csv','benchmark/cv_predictions.csv',
        'stability/stable_vs_original_EN_cv.csv','stability/full_train/bootstrap_all_features.csv',
        'stability/stable_Celltype_Gene.csv','stability/stable_Celltype.csv','stability/stable_Gene.csv',
        'stability/selection_decision.json','enrichment/pathway_enrichment.csv','enrichment/aging_sets_overlap.csv']]
    pins += [Path(m['path']) for m in models]+[Path(m['preprocessor']) for m in models]
    def info(p):return dict(path=str(p.resolve()),size=p.stat().st_size,sha256=file_hash(p),reuse='source_read_only')
    with ThreadPoolExecutor(max_workers=ctx.options['io_workers']) as pool:records=list(pool.map(info,pins))
    save('audit/reused_artifacts_manifest.csv',pd.DataFrame(records))
    p=ctx.path('audit/final_clock_source_pins.json');write_json(p,records);outputs.append(p)
    decision=ctx.json('stability/selection_decision.json');stable=ctx.csv('stability/stable_Celltype_Gene.csv')
    text=f'''# Phase3/4 actual gap audit: {ctx.lineage}

Development: {len(train)} samples, {train.group_id.nunique()} conservative GSE/GSM groups. Combined external: {len(external)} samples. No donor-like obs columns were found; donor independence is unavailable, not invented. Cached development/external sample/GSM/group overlap was checked and none remained. GSE representations overlapping combined are flagged, not independent confirmations.

Existing internal EN K={ctx.K}, highest mean fivefold Valid R2. Global benchmark winner was {decision['internal_winner']}; existing models/final_model is NOT a final EN clock. Stable results already exist: {len(stable)} Celltype_Gene, {stable.Celltype.nunique()} Celltypes, {stable.Gene.nunique()} unique Genes. Existing 300×(5 outer+full development) bootstrap records will be validated, not blindly assumed complete.

Already done / reuse: 13-model benchmark; sample-level sparse pseudobulk and explicit schema; group-wise fold models and inner-tuned EN parameters; bootstrap/stability tables; original/stable outer comparisons; versioned GMT/HGNC/provided aging sets and measured associations/enrichment.

Missing / repair: full-development original EN and stable EN side-by-side bundles; explicit final-EN decision and portable inference; final-EN no-adapt external and final-origin FT30; cross-lineage/direction enrichment; figures using the actual final clock instead of old best-fold FT30. All stages use new final_clock filenames in EXISTING directories, respecting the user's no-new-directory requirement. Existing benchmark/model artifacts remain pinned read-only.

Internal K is selected/descriptive; outer stability and preprocessing were trained only on outer TRAIN. Full-development stability is explanatory, never used as if it were unseen-fold feature selection. Prior external EN-best is exploratory and does not select K or the stable-model decision.

Figure upstreams, model inventory, cohort overlap, read references and source SHA256 manifests are saved beside this report. MAD is mean absolute AGE deviation; MedAE is median prediction error. Clinical/biological-age validation is not implied.
'''
    p=ctx.path('reports/phase3_4_gap_audit.md');write_text(p,text);outputs.append(p)
    env=[]
    for pkg in ['numpy','pandas','scipy','scikit-learn','celer','anndata','joblib','matplotlib','nbformat','psutil']:
        env.append(pkg+'=='+importlib.metadata.version(pkg))
    p=ctx.path('audit/final_clock_environment.txt');write_text(p,sys.executable+'\n'+'\n'.join(env)+'\n');outputs.append(p)
    LOG.info('%s audit saved: K=%d, %d development / %d combined; stable=%d',ctx.lineage,ctx.K,len(train),len(external),len(stable))
    return outputs
