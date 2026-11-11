#!/usr/bin/env python3
"""Bounded, training-study-only search over source-supported global marker lists.

No Leiden ID enters the search. Scores are always unweighted means of actual log
expression, with the same gene list used for every cell/cluster in a lineage.
Holdout studies are evaluated only after the selected list is frozen.
"""
from __future__ import annotations
import itertools
import logging
from pathlib import Path
import time
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse
import signature_pipeline as sp

VERSION='global_combinations_trainGSE_20261106'


def bank_and_group_scores(cfg,lin):
    c=sp.cache(cfg,lin);units=pd.read_parquet(c/'gene_group_units.parquet')
    ag=np.load(c/'gene_group_aggregates.npz');genes=np.load(c/'marker_genes.npy').tolist();lookup={g:i for i,g in enumerate(genes)}
    expression=ag['sums']/np.maximum(ag['counts'][:,None],1)
    detection=ag['detected']/np.maximum(ag['counts'][:,None],1)
    units['n_cells']=ag['counts']
    order=np.load(c/'final_score_names.npy').tolist()
    baseline=sp.read_json(c/'frozen_final_markers.json')['marker_sets']
    specs=sp.read_json(c/'score_specs_indexed.json')
    evidence=pd.read_csv(sp.root(cfg)/f'data/marker_selection_evidence_{lin}.csv')
    train=units['split'].eq('train').to_numpy()&(units['n_cells']>=cfg['min_cells_per_group_subtype']).to_numpy()
    pairs={s:[] for s in order}
    for g,frame in units[train].groupby('GSE',observed=True):
        indices={str(r.Celltype_general):i for i,r in frame.iterrows()}
        for s in order:
            if s in indices and len(indices)>1:pairs[s].append((indices[s],[i for k,i in indices.items() if k!=s]))
    bank={s:{} for s in order}
    def add(s,gg,source):
        gg=tuple(sorted(set(gg)&set(genes)))
        if not gg:return
        if gg not in bank[s]:bank[s][gg]={'genes':list(gg),'sources':[]}
        if source not in bank[s][gg]['sources']:bank[s][gg]['sources'].append(source)
    for item in specs:
        if item['mapping']=='projected' and item['key'] in bank:
            add(item['key'],item['genes'],item['candidate_id'])
    for s in order:
        add(s,baseline[s],'previous_frozen_training_QC_final')
        # Fixed, declared variants of the previous TRAINING gene-QC ranking.
        ev=evidence[evidence.subtype.eq(s)&evidence.passes_prespecified_QC].sort_values(
            ['min_pairwise_median_effect','positive_fraction','median_effect','gene'],ascending=[False,False,False,True])
        for n in [3,5,8]:add(s,ev.head(n).gene.tolist(),f'train_gene_QC_top_{n}')
        passed=set(ev.gene)
        for item in list(bank[s].values()):
            trimmed=[g for g in item['genes'] if g in passed]
            if len(trimmed)>=2:add(s,trimmed,'training_QC_intersection:'+','.join(item['sources']))
    score_bank={};qc=[];choices={};baseline_choice={}
    for s in order:
        entries=[];vectors=[]
        target_indices=np.flatnonzero(train&units.Celltype_general.eq(s).to_numpy())
        for key,item in bank[s].items():
            ii=[lookup[g] for g in item['genes']];v=expression[:,ii].mean(1)
            diff=np.array([v[ti]-v[oi].mean() for ti,oi in pairs[s]])
            cov=float(detection[target_indices][:,ii].mean())
            sm=sp.paired_summary(diff)
            eligible=(sm['n_pairs']>=10 and cov>=.05 and sm['median_effect']>=.02 and sm['positive_fraction']>=.60)
            is_baseline=set(item['genes'])==set(baseline[s])
            row={'lineage':lin,'subtype':s,'genes':';'.join(item['genes']),'n_genes':len(ii),
                 'sources':';'.join(item['sources']),'training_mean_gene_detection_GSE_equal':cov,
                 **sm,'eligible':eligible,'is_previous_baseline':is_baseline,
                 'QC_scope':'training GSE only; set-level paired mean-log contrast vs other available subtypes'}
            qc.append(row)
            if eligible or is_baseline:
                # Baseline is always available as an explicit comparator; failed QC is flagged.
                if is_baseline:baseline_choice[s]=len(entries)
                entries.append({**item,'training_QC_passed':eligible,'training_set_effect':sm,'training_detection':cov})
                vectors.append(v)
        if not entries:raise ValueError('No usable signatures for '+s)
        choices[s]=entries;score_bank[s]=np.column_stack(vectors)
    pd.DataFrame(qc).to_csv(sp.root(cfg)/f'data/marker_combination_candidate_QC_{lin}.csv',index=False)
    sp.dump(sp.root(cfg)/f'data/marker_combination_candidates_{lin}.json',
            {'order':order,'baseline_indices':baseline_choice,
             'index_meaning':'zero-based position in each subtype list; trial tuple follows order',
             'choices':choices})
    return units,order,baseline,choices,score_bank,tuple(baseline_choice[s] for s in order)


def evaluation_context(units,order,split):
    use=units.split.eq(split).to_numpy()&units.n_cells.ge(20).to_numpy()
    frame=units.loc[use].copy();y=pd.Categorical(frame.Celltype_general,categories=order).codes
    gcodes,gnames=pd.factorize(frame.GSE,sort=True);gcounts=np.bincount(gcodes)
    weights=1/(len(gnames)*gcounts[gcodes])
    return use,frame,y,gcodes,gnames,weights


def evaluate_matrix(v,context):
    use,frame,y,gcodes,gnames,weights=context;k=v.shape[1]
    rank=np.argsort(v,axis=1);pred=rank[:,-1]
    top=v[np.arange(len(v)),pred];second=v[np.arange(len(v)),rank[:,-2]]
    good=np.isfinite(top)&np.isfinite(second)&((top-second)>1e-10)
    match=(pred==y)&good
    recall=np.array([match[y==j].mean() if np.any(y==j) else 0 for j in range(k)])
    return {'GSE_equal_top1':float(weights@match),'balanced_unit_accuracy':float(recall.mean()),
            'min_subtype_recall':float(recall.min()),'n_GSE':len(gnames),'n_GSE_subtype_units':len(frame)},pred,match


def prepare(cfg,lin,resume=True):
    c=sp.cache(cfg,lin);r=sp.root(cfg);dest=r/f'marker_sets/{sp.LABELS[lin]}_subtype_markers_figure.json'
    scorepath=c/'figure_mean.npy';done=c/'figure_combinations_complete.json'
    if resume and done.exists() and dest.exists() and scorepath.exists():
        old=sp.read_json(done)
        doc=sp.read_json(dest)
        frozen=sp.read_json(c/'figure_markers_frozen_before_holdout.json')
        digest=hashlib.sha256(json.dumps(doc['marker_sets'],sort_keys=True).encode()).hexdigest()
        if (old['source_fingerprint']==sp.fingerprint(cfg['inputs'][lin]) and old['seed']==cfg['seed']
            and old.get('version')==VERSION and old['marker_sets_sha256']==digest
            and doc['marker_sets']==frozen['marker_sets']):
            return doc
    units,order,baseline,choices,bank,base=bank_and_group_scores(cfg,lin)
    trainctx=evaluation_context(units,order,'train');trainuse=trainctx[0]
    training={s:bank[s][trainuse] for s in order}
    history=[];memo={}
    def assess(combo):
        if combo in memo:return memo[combo]
        matrix=np.column_stack([training[s][:,j] for s,j in zip(order,combo)])
        met,_,_=evaluate_matrix(matrix,trainctx)
        ng=sum(len(choices[s][j]['genes']) for s,j in zip(order,combo))
        # Training classification criteria follow gene-coverage and cross-study-effect QC.
        objective=met['GSE_equal_top1']+.15*met['balanced_unit_accuracy']-.0001*ng
        rank=(objective,met['min_subtype_recall'],-ng)
        memo[combo]=(rank,met)
        history.append({'lineage':lin,'trial':len(history)+1,'combination_indices':';'.join(map(str,combo)),
                        **met,'training_objective':objective,'n_genes_total':ng,
                        'holdout_used':False,'Leiden_used':False})
        return rank,met
    sizes=[len(choices[s]) for s in order];total=int(np.prod(sizes));best=base;assess(base)
    if total<=10000:
        for combo in itertools.product(*(range(n) for n in sizes)):
            if assess(combo)[0]>assess(best)[0]:best=combo
        stop=f'exhausted prespecified candidate bank ({total} global combinations)'
    else:
        rng=np.random.default_rng(cfg['seed']);starts=[base]
        # Uniform source-dictionary families are included as deterministic starts.
        sources=sorted({src for s in order for item in choices[s] for src in item['sources'] if ':' in src})
        for source in sources:
            combo=[]
            for s in order:
                found=[j for j,x in enumerate(choices[s]) if source in x['sources']]
                if not found:break
                combo.append(found[0])
            if len(combo)==len(order):starts.append(tuple(combo))
        for _ in range(12):starts.append(tuple(int(rng.integers(n)) for n in sizes))
        for start in dict.fromkeys(starts):
            cur=start;assess(cur)
            for sweep in range(8):
                previous=cur
                for j,s in enumerate(order):
                    trials=[cur[:j]+(q,)+cur[j+1:] for q in range(sizes[j])]
                    cur=max(trials,key=lambda x:assess(x)[0])
                if cur==previous:break
            if assess(cur)[0]>assess(best)[0]:best=cur
        stop='bounded multi-start coordinate search; <=8 sweeps/start; stop at training objective stability, not at100%Leiden agreement'
    marker_sets={s:choices[s][j]['genes'] for s,j in zip(order,best)}
    doc={'lineage':sp.LABELS[lin],'version':VERSION,'annotation_column':'Celltype_general',
         'marker_sets':marker_sets,'scoring_method':'unweighted mean of matched genes in audited log-normalized X',
         'selection':{'scope':'training GSE only','unit':'GSE×reference subtype, >=20 cells',
           'candidate_counts':dict(zip(order,sizes)),'candidate_product':total,'trials_evaluated':len(history),
           'objective':'training GSE_equal_top1 +0.15 balanced unit accuracy -0.0001 total genes; tie break minimum class recall',
           'candidate_QC':'training paired studies>=10; mean per-gene target detection>=5%; paired median effect>=0.02; positive-study fraction>=60%; previous baseline explicitly retained as comparator',
           'stop_rule':stop,'Leiden_ID_used':False,'holdout_used':False,'seed':cfg['seed']},
         'selected_sources':{s:choices[s][j] for s,j in zip(order,best)},
         'source_fingerprint':sp.fingerprint(cfg['inputs'][lin]),
         'limitations':['Supervised selection uses marker-informed reference labels. Holdout isolates this search, not original annotation circularity.',
           'The existing holdout GSE were already examined in earlier baseline analyses; this is retrospective grouped evaluation, not a fresh blinded test.',
           'Different raw gene-set means are uncalibrated; no guarantee of all-cluster agreement or independent biological identity.',
           'One global list per subtype for every cell; no label-dependent score terms, cluster-specific gene lists, or removal of discordant clusters.']}
    # Frozen before holdout evaluation or any Leiden aggregation.
    sp.dump(c/'figure_markers_frozen_before_holdout.json',doc)
    pd.DataFrame(history).to_csv(r/f'data/marker_combination_search_{lin}.csv',index=False)
    comparison=[];status={}
    for split in ['train','holdout']:
        ctx=evaluation_context(units,order,split)
        for name,combo in [('previous_final',base),('selected_global_combination',best)]:
            matrix=np.column_stack([bank[s][ctx[0],j] for s,j in zip(order,combo)])
            met,pred,match=evaluate_matrix(matrix,ctx)
            groupmatches=np.bincount(ctx[3],weights=match)/np.bincount(ctx[3])
            lo,hi=sp.bootstrap_group_ci(groupmatches,cfg['seed'],cfg['bootstrap_replicates'])
            comparison.append({'lineage':lin,'split':split,'marker_version':name,**met,'GSE_bootstrap95_low':lo,'GSE_bootstrap95_high':hi})
            if split=='holdout' and name=='selected_global_combination':
                for j,s in enumerate(order):
                    take=ctx[2]==j;v=matrix[take];marg=v[:,j]-np.max(np.delete(v,j,axis=1),axis=1)
                    low,high=sp.bootstrap_group_ci(marg,cfg['seed'],cfg['bootstrap_replicates'])
                    n=int(take.sum());fraction=float(match[take].mean()) if n else np.nan
                    status[s]={'n_holdout_GSE':n,'GSE_unit_top1_fraction':fraction,
                               'own_minus_best_other_mean':float(marg.mean()) if n else np.nan,
                               'own_margin_bootstrap95CI':[low,high],
                               'status':'insufficient_evidence' if n<10 else ('supported_relative_to_reference' if fraction>=.65 and low>0 else 'ambiguous')}
    doc['holdout_validation']=status;doc['comparison']=comparison
    sp.dump(dest,doc);pd.DataFrame(comparison).to_csv(r/f'data/marker_combination_train_holdout_{lin}.csv',index=False)
    # Score all actual cells using expression alone after list freezing.
    genes=np.load(c/'marker_genes.npy').tolist();lookup={g:i for i,g in enumerate(genes)}
    x=np.load(c/'marker_expression.npy',mmap_mode='r')
    w=np.zeros((len(genes),len(order)),np.float32)
    for j,s in enumerate(order):
        for g in marker_sets[s]:w[lookup[g],j]=1/len(marker_sets[s])
    w=sparse.csr_matrix(w);output=np.lib.format.open_memmap(scorepath,mode='w+',dtype=np.float32,shape=(len(x),len(order)))
    for st in range(0,len(x),65536):
        en=min(len(x),st+65536);output[st:en]=(w.T@x[st:en].T).T
    output.flush();np.save(c/'figure_score_names.npy',np.asarray(order))
    sp.dump(done,{'source_fingerprint':doc['source_fingerprint'],'seed':cfg['seed'],'version':VERSION,
                  'n_cells':len(x),'marker_sets_sha256':hashlib.sha256(json.dumps(marker_sets,sort_keys=True).encode()).hexdigest(),
                  'computed_from':'unweighted cached log expression; no annotation columns used in cell scoring'})
    logging.info('%s combination search: %s trials, train %.3f→%.3f, holdout %.3f→%.3f; no Leiden optimization',
                 lin,len(history),comparison[0]['GSE_equal_top1'],comparison[1]['GSE_equal_top1'],comparison[2]['GSE_equal_top1'],comparison[3]['GSE_equal_top1'])
    return doc


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'config.yaml'))
    p.add_argument('--lineage',choices=['fibro','endo','both'],default='both');p.add_argument('--resume',action='store_true');a=p.parse_args()
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s');cfg=sp.read_config(a.config)
    for lin in ['fibro','endo'] if a.lineage=='both' else [a.lineage]:prepare(cfg,lin,resume=a.resume)
