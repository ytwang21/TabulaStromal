#!/usr/bin/env python
"""Expression-only implementation of fig3.ipynb cells 19--29, 34 and 39.

Uses existing pseudobulk caches. Each scope freezes F-test top8000/scaling,
tunes initial EN, resamples groups 300 times, selects frequency>=0.70, and
retunes reduced EN. Holdout and full-development fits are separate experiments.
"""
import os
for _key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[_key]='1'
from pathlib import Path
import argparse,copy,json,logging,subprocess,sys,time,warnings,re,functools
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import VarianceThreshold,SelectKBest,f_regression
from sklearn.linear_model import ElasticNet,ElasticNetCV
from sklearn.model_selection import GroupShuffleSplit,KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score,mean_absolute_error,median_absolute_error
from joblib import Parallel,delayed,parallel_config
from threadpoolctl import threadpool_limits
from final_clock_common import (ROOT,ROOTS,Context,dataset,make_bundle,bundle_predict,
    load_final_clock,predict_age,now,file_hash,digest,write_json,write_csv,write_text,dump,fit_checked)
from tabulastromal.storage import process_lock

REFERENCE=Path('external_inputs/fig3.ipynb')
PREFIX='fig3_reference'
SPEC=dict(protocol='fig3-fixed-pool-group-bootstrap-v1',reference_cells=[19,22,23,24,25,26,27,28,29,34,39],
    k_best=8000,variance_threshold=1e-5,n_bootstrap=300,stable_threshold=0.70,seed=42,
    initial_l1=[0.5,0.7,0.9,0.95,0.99],final_l1=[0.1,0.5,0.7,0.9,0.95,0.99,1.0],
    cv_splits=5,n_alphas=100,eps=0.001,tol=1e-4,bootstrap_unit='GSE|GSM proxy group',
    holdout_fraction_groups=0.2,bootstrap_support='coef != 0 exactly',
    external_policy='primary: no adaptation, predict ALL100',MAD='median absolute prediction error')
LOG=logging.getLogger(PREFIX)
CV_SOLVER=os.environ.get('FIG3_REFERENCE_CV_SOLVER','celer')

@functools.lru_cache(maxsize=1)
def reference_identity():
    """Pin reviewed algorithm cells; unrelated later notebook edits do not invalidate fits."""
    import nbformat
    original=json.loads((ROOT/'fig3_reference_source_audit.json').read_text())
    snapshot=ROOT/'fig3_reference_source_extract.py'
    assert file_hash(snapshot)==original['extract_sha256'],'Reviewed source snapshot was modified'
    parts=re.split(r'# REFERENCE CELL (\d+)\n',snapshot.read_text())
    old={int(parts[i]):parts[i+1].rstrip() for i in range(1,len(parts),2)}
    current=nbformat.read(REFERENCE,as_version=4)
    changed=[i for i,src in old.items() if i>=len(current.cells) or src!=current.cells[i].source.rstrip()]
    if changed:raise RuntimeError(f'Reference algorithm cells changed; explicit review required: {changed}')
    return original['sha256']

def metrics(y,p):
    y=np.asarray(y,float);p=np.asarray(p,float)
    assert y.ndim==p.ndim==1 and y.shape==p.shape and len(y)>0
    assert np.isfinite(y).all() and np.isfinite(p).all()
    e=y-p;sse=float(e@e);sst=float(((y-y.mean())**2).sum())
    m=dict(n=len(y),R2=1-sse/sst if sst>0 and len(y)>1 else np.nan,
        MAE=float(np.abs(e).mean()),MAD=float(np.median(np.abs(e))),MedAE=float(np.median(np.abs(e))),
        Age_MAD=float(np.abs(y-y.mean()).mean()),RMSE=float(np.sqrt(np.mean(e**2))),
        PCC=float(pearsonr(y,p).statistic) if len(y)>1 and np.ptp(y)>0 and np.ptp(p)>0 else np.nan,
        SSE=sse,SST=sst)
    np.testing.assert_allclose(m['R2'],r2_score(y,p),rtol=1e-11,atol=1e-11)
    np.testing.assert_allclose([m['MAE'],m['MAD']],[mean_absolute_error(y,p),median_absolute_error(y,p)])
    return m

def group_cv(groups,shuffle):
    groups=np.asarray(groups);unique=pd.unique(groups)
    cv=KFold(n_splits=5,shuffle=shuffle,random_state=42 if shuffle else None)
    out=[]
    for a,b in cv.split(unique):
        train=np.flatnonzero(np.isin(groups,unique[a]));valid=np.flatnonzero(np.isin(groups,unique[b]))
        assert not set(groups[train])&set(groups[valid]);out.append((train,valid))
    return out

def cached(path,signature):
    try:
        p=Path(path);m=json.loads(p.with_suffix(p.suffix+'.json').read_text())
        return m['signature']==signature and m['sha256']==file_hash(p)
    except (OSError,KeyError,ValueError):return False

def cache(path,value,signature):
    dump(path,value);write_json(str(path)+'.json',dict(signature=signature,sha256=file_hash(path),created_at=now()))

def fit_cv(X,y,groups,l1,n_jobs,label):
    LOG.info('TUNE %s n=%d p=%d l1=%s',label,len(y),X.shape[1],l1)
    # Identical alpha path/default tolerance to sklearn1.3 ElasticNetCV in reference.
    records=[]
    limits=[100,500,2500] if CV_SOLVER=='celer' else [5000,25000,125000]
    for limit in limits:
        options=dict(l1_ratio=l1,cv=group_cv(groups,shuffle=label!='initial'),
            n_alphas=100,eps=.001,max_iter=limit,n_jobs=n_jobs,tol=1e-4)
        if CV_SOLVER=='celer':
            from celer import ElasticNetCV as FastElasticNetCV
            options['tol']=1e-7  # stricter accuracy for fast path solver
            model=FastElasticNetCV(**options,max_epochs=50000)
        else:model=ElasticNetCV(**options,random_state=42)
        start=time.time()
        with threadpool_limits(1),warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always');model.fit(np.asfortranarray(X,dtype=float),y)
        messages=[str(w.message) for w in ws if issubclass(w.category,ConvergenceWarning)]
        records.append(dict(solver=CV_SOLVER,solver_tol=model.tol,max_iter=limit,seconds=time.time()-start,n_convergence_warnings=len(messages),examples=messages[:3]))
        if not messages:break
        LOG.warning('%s: %d convergence warnings, retry with higher iteration limit',label,len(messages))
    else:raise RuntimeError(f'Unconverged CV {label}: {records[-1]}')
    model.fig3_reference_fit_audit_=records
    LOG.info('TUNED %s alpha=%.9g l1=%.3g gap=%.5g',label,model.alpha_,model.l1_ratio_,model.dual_gap_)
    return model

def bootstrap_task(i,X,y,groups,alpha,l1,path,signature):
    if cached(path,signature):return str(path)
    unique=pd.unique(groups);draw=np.random.RandomState(i).choice(len(unique),len(unique),replace=True)
    blocks=[np.flatnonzero(groups==g) for g in unique]
    ix=np.concatenate([blocks[j] for j in draw])
    model=fit_checked(ElasticNet(alpha=alpha,l1_ratio=l1,max_iter=5000,random_state=42),X[ix],y[ix])
    result=dict(iteration=i,coef=model.coef_,intercept=float(model.intercept_),draw_groups=unique[draw].tolist(),
        drawn_rows=len(ix),drawn_group_count=len(draw),dual_gap=float(model.dual_gap_),n_iter=int(model.n_iter_),
        alpha=alpha,l1_ratio=l1,warnings=getattr(model,'final_clock_fit_warnings_',[]))
    cache(path,result,signature)
    return str(path)

def subset_scaler(scaler,mask):
    result=copy.deepcopy(scaler)
    for field in ['mean_','var_','scale_']:setattr(result,field,getattr(scaler,field)[mask])
    result.n_features_in_=int(np.sum(mask));return result

def feature_table(schema,indices):
    n=len(schema['genes'])
    rows=[dict(feature=schema['features'][j],Celltype=schema['celltypes'][j//n],Gene=schema['genes'][j%n],source_column_index=int(j)) for j in indices]
    return pd.DataFrame(rows)

def run_scope(ctx,ds,scope,fitrows,validrows,args,signature):
    p=lambda section,suffix:ctx.path(f'{section}/{PREFIX}_{scope}_{suffix}')
    meta=ds['meta'];y=meta.Age.to_numpy(float);groups=meta.group_id.to_numpy(str)
    train=meta.iloc[fitrows];g=groups[fitrows];yy=y[fitrows]
    scope_sig=digest(dict(signature=signature,scope=scope,fit_ids=train.sample_id.tolist()))
    prepath=p('checkpoints','preprocessor.joblib')
    if cached(prepath,scope_sig):pre=joblib.load(prepath)
    else:
        LOG.info('%s %s selecting frozen F-test top8000 from %d × %d',ctx.lineage,scope,len(fitrows),ds['X'].shape[1])
        raw=ds['X'][fitrows].astype(np.float64)
        vt=VarianceThreshold(threshold=SPEC['variance_threshold']);filtered=vt.fit_transform(raw)
        skb=SelectKBest(f_regression,k=min(8000,filtered.shape[1]));selected=skb.fit_transform(filtered,yy)
        ix=np.flatnonzero(vt.get_support())[skb.get_support()]
        scaler=StandardScaler().fit(selected.toarray())
        pre=dict(indices=ix,scaler=scaler,fit_ids=train.sample_id.tolist(),fit_groups=sorted(set(g)),
            variance_pass_indices=np.flatnonzero(vt.get_support()),f_scores=skb.scores_,f_pvalues=skb.pvalues_)
        cache(prepath,pre,scope_sig)
    assert pre['fit_ids']==train.sample_id.tolist()
    ix=pre['indices'];scaler=pre['scaler'];raw=ds['X'][:,ix].toarray();X=scaler.transform(raw[fitrows])
    initialpath=p('models','initial_EN.joblib')
    if cached(initialpath,scope_sig):initial=joblib.load(initialpath)
    else:
        initial=fit_cv(X,yy,g,SPEC['initial_l1'],args.cv_jobs,'initial');cache(initialpath,initial,scope_sig)
    paths=[p('stability',f'boot_{i:04d}.joblib') for i in range(300)]
    LOG.info('%s %s bootstrap cached=%d/300 workers=%d',ctx.lineage,scope,sum(cached(q,scope_sig) for q in paths),args.workers)
    with parallel_config(backend='loky',inner_max_num_threads=1):
        Parallel(n_jobs=args.workers,batch_size=1,verbose=10)(delayed(bootstrap_task)(i,X,yy,g,
            float(initial.alpha_),float(initial.l1_ratio_),paths[i],scope_sig) for i in range(300))
    boots=[joblib.load(q) for q in paths];coef=np.stack([b['coef'] for b in boots]);support=coef!=0
    table=feature_table(ds['schema'],ix);table['selection_count']=support.sum(axis=0)
    table['selection_frequency']=support.mean(axis=0);table['coef_mean']=coef.mean(axis=0)
    table['positive_count']=(coef>0).sum(axis=0);table['negative_count']=(coef<0).sum(axis=0)
    table['sign_consistency']=np.maximum(table.positive_count,table.negative_count)/np.maximum(table.selection_count,1)
    stable=table.selection_count.to_numpy()>=210
    assert np.array_equal(stable,table.selection_frequency.to_numpy()>=.7)
    if not stable.any():raise RuntimeError('No stable features at pre-specified0.70; threshold will not be changed using validation data')
    table['stable']=stable;table['in_fixed_candidate_pool']=True
    write_csv(p('stability','selection_frequency.csv'),table.sort_values(['selection_frequency','feature'],ascending=[False,True]))
    write_csv(p('stability','stable_features.csv'),table[stable].sort_values(['selection_frequency','feature'],ascending=[False,True]))
    for field,label in [('Celltype','celltypes'),('Gene','genes')]:
        counts=table[stable].groupby(field).agg(n_stable_features=('feature','size'),mean_frequency=('selection_frequency','mean')).reset_index()
        write_csv(p('stability',f'stable_{label}.csv'),counts.sort_values(['n_stable_features',field],ascending=[False,True]))
    write_csv(p('stability','bootstrap_manifest.csv'),pd.DataFrame([dict(iteration=i,path=str(q.resolve()),sha256=file_hash(q),
        n_drawn_rows=b['drawn_rows'],n_drawn_groups=b['drawn_group_count'],dual_gap=b['dual_gap'],n_iter=b['n_iter'],
        n_nonzero=int(np.count_nonzero(b['coef'])),n_warnings=len(b['warnings'])) for i,(q,b) in enumerate(zip(paths,boots))]))
    LOG.info('%s %s STABLE %d features, %d celltypes, %d genes',ctx.lineage,scope,stable.sum(),table[stable].Celltype.nunique(),table[stable].Gene.nunique())
    rows=[];predictions=[];bundles={}
    for name,mask in [('initial_EN',np.ones(len(ix),bool)),('all8000_EN',np.ones(len(ix),bool)),('stable_EN',stable)]:
        modelpath=p('models',f'{name}.joblib')
        if name=='initial_EN':model=initial
        elif cached(modelpath,scope_sig):model=joblib.load(modelpath)
        else:
            model=fit_cv(X[:,mask],yy,g,SPEC['final_l1'],args.cv_jobs,name);cache(modelpath,model,scope_sig)
        bundle=make_bundle(ctx,model,subset_scaler(scaler,mask),ix[mask],ds['schema'],f'{ctx.lineage}_{PREFIX}_{scope}_{name}_{scope_sig[:10]}')
        bundle.update(K=8000,training_scope=scope,training_ids=train.sample_id.tolist(),training_groups=sorted(set(g)),
            protocol=SPEC,reference_notebook=str(REFERENCE),reference_sha256=reference_identity(),experiment_signature=signature,
            final_alpha=float(model.alpha_),final_l1_ratio=float(model.l1_ratio_),scope_signature=scope_sig)
        bundlepath=p('models',f'{name}_bundle.joblib');dump(bundlepath,bundle);bundles[name]=str(bundlepath.resolve())
        # Source arrays retain float64 precision; public bundle inference rounds input
        # to float32, exactly matching source cache dtype (its values are float32).
        pp=bundle_predict(bundle,raw[:,mask]);np.testing.assert_allclose(model.predict(scaler.transform(raw)[:,mask]),pp,atol=2e-5,rtol=2e-6)
        np.testing.assert_allclose(pp,bundle_predict(load_final_clock(bundlepath),raw[:,mask]),atol=1e-12)
        parts=[('Train',fitrows)]+([('Valid',validrows)] if len(validrows) else [])
        for split,rr in parts:
            rows.append(dict(lineage=ctx.lineage,scope=scope,model=name,split=split,n_features=int(mask.sum()),
                alpha=float(model.alpha_),l1_ratio=float(model.l1_ratio_),model_path=str(bundlepath.resolve()),**metrics(y[rr],pp[rr])))
            q=meta.iloc[rr].copy();q['prediction']=pp[rr];q['scope']=scope;q['model']=name;q['split']=split;predictions.append(q)
    write_csv(p('stability','metrics.csv'),pd.DataFrame(rows));write_csv(p('stability','predictions.csv'),pd.concat(predictions,ignore_index=True))
    write_csv(p('stability','sample_split.csv'),meta.assign(split=np.where(np.isin(np.arange(len(meta)),fitrows),'Train','Valid')))
    audit=dict(signature=scope_sig,status='complete',n_bootstrap=300,n_candidate_features=len(ix),n_stable_features=int(stable.sum()),
        n_stable_celltypes=int(table[stable].Celltype.nunique()),n_stable_genes=int(table[stable].Gene.nunique()),
        fit_rows=len(fitrows),validation_rows=len(validrows),fit_ids=train.sample_id.tolist(),valid_ids=meta.iloc[validrows].sample_id.tolist(),
        no_group_overlap=not bool(set(g)&set(groups[validrows])),initial_alpha=float(initial.alpha_),initial_l1=float(initial.l1_ratio_),
        initial_fit_audit=initial.fig3_reference_fit_audit_,bundles=bundles,created_at=now())
    write_json(p('checkpoints','complete.json'),audit)
    return audit

def external(ctx,signature):
    ds=dataset(ctx,'combined');meta=ds['meta'];y=meta.Age.to_numpy(float);out=[];scores=[]
    trainmeta=pd.read_csv(ctx.path('features/train/samples.csv'))
    assert not set(meta.sample_id)&set(trainmeta.sample_id)
    assert not set(meta.group_id)&set(trainmeta.group_id)
    for name in ['initial_EN','all8000_EN','stable_EN']:
        path=ctx.path(f'models/{PREFIX}_full_{name}_bundle.joblib');bundle=load_final_clock(path)
        assert bundle['experiment_signature']==signature
        pred=predict_age(ds['path'],path)
        assert pred.sample_id.tolist()==meta.sample_id.tolist()
        q=meta.copy();q['prediction']=pred.age_pred.to_numpy();q['model']=name;q['split']='External_ALL100_no_adaptation'
        q['scope']='full';out.append(q)
        scores.append(dict(lineage=ctx.lineage,scope='full',split=q['split'].iloc[0],model=name,n_features=len(bundle['selected_features']),
            model_path=str(path.resolve()),**metrics(y,q.prediction)))
    write_csv(ctx.path(f'external/{PREFIX}_combined_predictions.csv'),pd.concat(out,ignore_index=True))
    write_csv(ctx.path(f'external/{PREFIX}_combined_metrics.csv'),pd.DataFrame(scores))

def run_lineage(args,lineage):
    global CV_SOLVER
    CV_SOLVER=args.cv_solver
    ctx=Context(lineage)
    logging.basicConfig(level=logging.INFO,format='%(asctime)s | %(levelname)s | %(message)s',handlers=[
        logging.StreamHandler(),logging.FileHandler(ctx.path(f'logs/{PREFIX}.log'))],force=True)
    with process_lock(ctx.path(f'checkpoints/{PREFIX}.lock')):
        ds=dataset(ctx);inputs=[ds['path']/f for f in ['pseudobulk.npz','samples.csv','schema.json']]
        inputs += [ctx.path('features/combined')/f for f in ['pseudobulk.npz','samples.csv','schema.json']]
        signature=digest(dict(spec=SPEC,source=[dict(path=str(p),sha256=file_hash(p)) for p in inputs],reference_sha256=reference_identity()))
        write_json(ctx.path(f'checkpoints/{PREFIX}_complete.json'),dict(status='running',signature=signature,
            started_at=now(),note='Cached fits are validated; this marker becomes complete only after external output is regenerated.'))
        write_json(ctx.path(f'audit/{PREFIX}_protocol.json'),dict(spec=SPEC,signature=signature,inputs=[dict(path=str(p),sha256=file_hash(p)) for p in inputs],
            reference_sha256=reference_identity(),reference_current_file_sha256=file_hash(REFERENCE),reviewed_cells_unchanged=True,
            script_sha256=file_hash(__file__),started_at=now(),
            deviations=['Expression Celltype_Gene only; no ratios','Group-aware20% holdout, tuning CV and bootstrap to retain repeated GSM together',
                'Convergence warnings trigger max_iter increase, never hyperparameter/threshold changes',
                'Reported generalization is untouched20% holdout, not selection-conditioned inner CV',
                'Full-development selection and model refit separately, for genuine future/external prediction',
                'MAD=median absolute prediction error to match reference; Age_MAD stored separately'],
            cv_solver_for_uncached_fits=CV_SOLVER,bootstrap_solver='sklearn ElasticNet (unchanged)'))
        splitter=GroupShuffleSplit(n_splits=1,test_size=.2,random_state=42)
        fitrows,validrows=next(splitter.split(ds['X'],ds['meta'].Age,groups=ds['meta'].group_id))
        eval_result=run_scope(ctx,ds,'holdout',fitrows,validrows,args,signature)
        final_result=run_scope(ctx,ds,'full',np.arange(len(ds['meta'])),np.array([],dtype=int),args,signature)
        external(ctx,signature)
        write_json(ctx.path(f'checkpoints/{PREFIX}_complete.json'),dict(status='complete',signature=signature,
            holdout=eval_result,full=final_result,completed_at=now()))
    LOG.info('%s ALL COMPLETE',lineage)

def render_report():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'font.size':10})
    fig,axes=plt.subplots(4,2,figsize=(13,18),layout='constrained')
    summaries=[];allmetrics=[];files=[]
    for col,lineage in enumerate(['Fibro','Endo']):
        ctx=Context(lineage);title={'Fibro':'Fibroblast','Endo':'Endothelial'}[lineage]
        hold=ctx.csv(f'stability/{PREFIX}_holdout_metrics.csv');ext=ctx.csv(f'external/{PREFIX}_combined_metrics.csv')
        m=pd.concat([hold,ext],ignore_index=True);allmetrics.append(m)
        stable=ctx.csv(f'stability/{PREFIX}_full_stable_features.csv');freq=ctx.csv(f'stability/{PREFIX}_full_selection_frequency.csv')
        ax=axes[0,col];ax.hist(freq.selection_frequency,bins=np.linspace(0,1,31),color='#4677AA')
        ax.axvline(.7,color='#D55E00',ls='--');ax.set(xlabel='Bootstrap selection frequency (300 runs)',ylabel='Features',title=f'{title}: {len(stable)} stable features')
        valid=hold[hold.split=='Valid'].set_index('model').loc[['all8000_EN','stable_EN']]
        ax=axes[1,col];bars=ax.bar(['All8000 EN','Stable EN'],valid.R2,color=['#4677AA','#D55E00']);ax.axhline(0,color='gray',lw=.7)
        ax.bar_label(bars,fmt='%.3f',padding=3);ax.set(ylabel='R²',title='Untouched internal20% validation (group split)')
        hp=ctx.csv(f'stability/{PREFIX}_holdout_predictions.csv');ep=ctx.csv(f'external/{PREFIX}_combined_predictions.csv')
        for row,frame,split in [(2,hp[hp.split=='Valid'],'Internal Valid'),(3,ep,'External ALL100; no adaptation')]:
            ax=axes[row,col]
            texts=[]
            for name,color,label in [('all8000_EN','#4677AA','All8000'),('stable_EN','#D55E00','Stable')]:
                part=frame[frame.model==name];ax.scatter(part.Age,part.prediction,s=16,alpha=.65,color=color,label=label)
                s=metrics(part.Age,part.prediction);texts.append(f'{label}: R²={s["R2"]:.3f}, MAE={s["MAE"]:.2f}, MAD={s["MAD"]:.2f}')
            lo=min(frame.Age.min(),frame.prediction.min(),0);hi=max(frame.Age.max(),frame.prediction.max(),100)
            ax.plot([lo,hi],[lo,hi],ls='--',c='gray',lw=1);ax.set(xlabel='Actual age (years)',ylabel='Predicted age (years)',title=split)
            ax.legend(loc='lower right');ax.text(.03,.97,'\n'.join(texts),transform=ax.transAxes,va='top',fontsize=8,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
        completed=ctx.json(f'checkpoints/{PREFIX}_complete.json')
        summaries.append(dict(lineage=lineage,n_stable_features=len(stable),n_stable_celltypes=stable.Celltype.nunique(),n_stable_genes=stable.Gene.nunique(),
            final_model=completed['full']['bundles']['stable_EN'],validation_stable_R2=float(valid.loc['stable_EN','R2']),
            validation_all8000_R2=float(valid.loc['all8000_EN','R2']),validation_stable_MAE=float(valid.loc['stable_EN','MAE']),
            validation_all8000_MAE=float(valid.loc['all8000_EN','MAE'])))
        for row in range(4):axes[row,col].text(-.12,1.07,chr(65+row*2+col),transform=axes[row,col].transAxes,fontweight='bold',fontsize=15)
    output=ROOT/'Stromal_Combined_Figure/outputs'
    for suffix in ['png','pdf']:
        p=output/f'Stromal_fig3_reference_bootstrap.{suffix}';fig.savefig(p,dpi=350,bbox_inches='tight');files.append(p)
    plt.close(fig)
    write_csv(output/f'{PREFIX}_model_comparison.csv',pd.concat(allmetrics,ignore_index=True))
    write_csv(ROOT/f'{PREFIX}_results_summary.csv',pd.DataFrame(summaries))
    lines=['# fig3.ipynb reference Bootstrap experiment','',f'Reviewed reference: `{REFERENCE}`; original SHA256 `{reference_identity()}`. Algorithm cells are checked against the saved source extract; unrelated added cells/outputs do not invalidate checkpoints.','',
        'Read cells19--29 plus complete functions34/39. Frozen training-only VarianceThreshold>1e-5 → SelectKBest(f_regression,8000) → StandardScaler → initial ElasticNetCV(5 folds, l1=.5/.7/.9/.95/.99) →300 replacement bootstraps with fixed alpha/l1 and frozen pool/scaler →coef!=0 frequency>=.70 →stable-only ElasticNetCV(l1=.1/.5/.7/.9/.95/.99/1).',
        '', 'GSE|GSM proxy groups replace row-based splitting/resampling, keeping repeated GSM in one unit. True donor IDs are unavailable. Grouping is not an independent GSE holdout. No cell proportions are used. Holdout20% is never used to select features/hyperparameters. Reference cross_val_score after whole-training feature selection is selection-conditioned; it is not reported as unbiased validation.','',
        'The stable final model is reselected/refitted using ALL development samples. Its external predictions use ALL combined external samples without adaptation, to assess transfer. No external age is used in fitting/feature selection. This is a new experiment, not a replacement for13-model benchmark results. Earlier FINAL_v2 outputs remain distinct.','',
        '**MAD here follows reference notebook cell30: median(abs(actual-prediction)); MedAE is identical. Age_MAD is the separate age-dispersion statistic used by historical reports.** MAE, R², PCC and MAD always use the same complete sample pairs.','',
        'For each lineage, `holdout` files concern the independent internal20% validation experiment; `full` files concern the deployment model. Their stable lists differ by design. `all8000_EN` uses the same final tuning grid as stable_EN for a paired comparison; `initial_EN` uses the narrower initial grid.','',
        'All stored models are self-contained joblib bundles, with explicit Celltype/Gene mapping and frozen selected scaler. Load with `final_clock_common.load_final_clock` and predict with `final_clock_common.predict_age`; this handles named pseudobulk or raw h5ad.','']
    lines += ['CPU acceleration: uncached ElasticNetCV fits use the locally vendored Celer solver for the same ElasticNet objective,100-alpha path and folds, with stricter solver tolerance1e-7. Previously converged sklearn fits are reused. Bootstrap coefficients are still fitted with sklearn ElasticNet exactly as before. Solver/tolerance are recorded in each tuning model fig3_reference_fit_audit_. `--cv-solver sklearn` is available for new uncached fits; it does not invalidate valid mathematical-protocol checkpoints.','']
    for s in summaries:
        root=ROOTS[s['lineage']];stable=pd.read_csv(root/f'stability/{PREFIX}_full_stable_features.csv')
        lines += [f'## {s["lineage"]}', '',f'Final stable model: `{s["final_model"]}`',
            f'Stable Celltype_Gene={s["n_stable_features"]}; Celltypes={s["n_stable_celltypes"]}; Genes={s["n_stable_genes"]}.',
            f'Holdout R² all8000/stable={s["validation_all8000_R2"]:.6f}/{s["validation_stable_R2"]:.6f}; MAE={s["validation_all8000_MAE"]:.6f}/{s["validation_stable_MAE"]:.6f}.','',
            f'Full lists: `{root}/stability/{PREFIX}_full_stable_features.csv`, `_stable_celltypes.csv`, `_stable_genes.csv`.',
            f'External per-sample ages: `{root}/external/{PREFIX}_combined_predictions.csv`.',
            f'Bootstrap300 checkpoints/scope: `{root}/stability/{PREFIX}_{{holdout,full}}_boot_0000.joblib` … `0299.joblib`.',
            f'Log: `{root}/logs/{PREFIX}.log`.','', 'Stable celltypes: '+', '.join(sorted(stable.Celltype.unique())),'']
    lines += ['## Re-run / resume','','```bash',f'cd {ROOT}', 'mkdir -p logs',
        'nohup external_inputs/python -u run_fig3_reference_bootstrap.py --lineage both --workers 24 --cv-jobs 25 >> logs/fig3_reference_both.nohup.log 2>&1 < /dev/null &',
        'tail -f logs/fig3_reference_both.nohup.log', '```','',
        'Each scope caches preprocessing, initial/final tuning and all300 bootstrap fits with input/reference signatures and SHA256 checks. Resume is automatic. Separate advisory locks prevent duplicate writers. Two lineage processes run concurrently; bootstrap workers limit BLAS to one thread. This is sklearn CPU ElasticNet, so an unverified GPU estimator is not substituted.','']
    write_text(ROOT/'README_FIG3_REFERENCE_BOOTSTRAP.md','\n'.join(lines))
    write_json(ROOT/f'{PREFIX}_figure_complete.json',dict(status='complete',created_at=now(),outputs=[dict(path=str(p),sha256=file_hash(p)) for p in files]))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--lineage',choices=['both','Fibro','Endo'],default='both')
    ap.add_argument('--workers',type=int,default=24);ap.add_argument('--cv-jobs',type=int,default=25);ap.add_argument('--report-only',action='store_true')
    ap.add_argument('--cv-solver',choices=['celer','sklearn'],default=CV_SOLVER)
    args=ap.parse_args()
    if args.report_only:render_report();return
    if args.lineage!='both':run_lineage(args,args.lineage);return
    handles=[];children=[]
    try:
        with process_lock(ROOT/f'logs/{PREFIX}_both.lock'):
            for lineage in ['Fibro','Endo']:
                handle=open(ROOT/f'logs/{PREFIX}_{lineage}.launcher.log','a',buffering=1);handles.append(handle)
                child=subprocess.Popen([sys.executable,'-u',str(Path(__file__).resolve()),'--lineage',lineage,'--workers',str(args.workers),'--cv-jobs',str(args.cv_jobs),'--cv-solver',args.cv_solver],stdout=handle,stderr=subprocess.STDOUT)
                children.append(child);print(f'{lineage} PID={child.pid}; log={handle.name}',flush=True)
            failures=[]
            for lineage,child in zip(['Fibro','Endo'],children):
                rc=child.wait()
                if rc:failures.append((lineage,rc))
            if failures:raise RuntimeError(f'Lineage failures (check launcher logs): {failures}')
            render_report();print('BOTH COMPLETE',flush=True)
    finally:
        for handle in handles:handle.close()

if __name__=='__main__':main()
