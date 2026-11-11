"""Measured feature annotation, universe-aware ORA and covariate-adjusted associations."""
from pathlib import Path
import json,logging,re
import numpy as np
import pandas as pd
from scipy.stats import hypergeom,t as student_t
from .storage import write_csv,write_json,write_text,file_hash
from .data import load_dataset
LOG=logging.getLogger(__name__)


def bh(p):
    p=np.asarray(p,float);out=np.full(len(p),np.nan);valid=np.flatnonzero(np.isfinite(p))
    if not len(valid):return out
    order=valid[np.argsort(p[valid])];adj=p[order]*len(order)/np.arange(1,len(order)+1)
    out[order]=np.minimum(1,np.minimum.accumulate(adj[::-1])[::-1]);return out


def hgnc_index(path):
    if not path or not Path(path).is_file():return None,{}
    table=pd.read_csv(path,sep='\t',dtype=str).fillna('');approved=set(table.loc[table.status=='Approved','symbol']) if 'status' in table else set(table.symbol)
    aliases={}
    for row in table.itertuples():
        if row.symbol not in approved:continue
        for field in ['alias_symbol','prev_symbol']:
            for alias in str(getattr(row,field,'')).split('|'):
                if alias:aliases.setdefault(alias,set()).add(row.symbol)
    return approved,aliases


def qc_sets(ctx):
    source=json.loads(Path(ctx.cfg['enrichment']['provided_sets']).read_text());approved,aliases=hgnc_index(ctx.cfg['enrichment']['hgnc_tsv'])
    ids={}
    if approved is not None:
        hgnc=pd.read_csv(ctx.cfg['enrichment']['hgnc_tsv'],sep='\t',dtype=str).fillna('')
        ids=dict(zip(hgnc.symbol,hgnc.hgnc_id))
    qc=[];clean={}
    for name,genes in source['sets'].items():
        clean[name]=set();seen=set()
        for original in genes:
            symbol=str(original);mapped=None;status=''
            if symbol in seen:status='duplicate_original'
            elif re.search(r'\s',symbol):status='invalid_whitespace_requires_manual_review'
            elif approved is None:
                status='HGNC_unavailable_unverified_exact_symbol';mapped=symbol if re.fullmatch(r'[A-Z0-9][A-Z0-9.\-]*',symbol) else None
                if mapped is None:status='mixed_case_or_nonstandard_requires_review'
            elif symbol in approved:mapped=symbol;status='approved_exact'
            elif symbol.upper() in approved:mapped=symbol.upper();status='case_mapped_approved'
            elif symbol in aliases and len(aliases[symbol])==1:mapped=next(iter(aliases[symbol]));status='unique_HGNC_alias'
            elif symbol in aliases:status='ambiguous_alias'
            else:status='unmapped'
            duplicate_mapped=mapped in clean[name] if mapped else False
            if mapped:clean[name].add(mapped)
            qc.append(dict(set=name,original=symbol,mapped_symbol=mapped or '',hgnc_id=ids.get(mapped,''),duplicate_mapped_symbol=duplicate_mapped,status=status,source_label='provided_unverified',source_path=source['source_path']))
            seen.add(symbol)
    write_csv(ctx.root/'enrichment/provided_gene_sets_qc.csv',pd.DataFrame(qc));write_json(ctx.root/'enrichment/provided_gene_sets_cleaned.json',{k:sorted(v) for k,v in clean.items()})
    limitations=[]
    if approved is None:limitations.append('HGNC mapping file unavailable: exact uppercase symbols are only syntax-checked, not officially approved; mixed-case/ambiguous terms are not silently corrected.')
    return clean,limitations


def load_libraries(ctx):
    base=Path(ctx.cfg['enrichment']['gmt_dir']);libs={};provenance=[]
    manifest=base/'download_manifest.json';records=json.loads(manifest.read_text()) if manifest.exists() else []
    by_version={r['Version']:r for r in records if 'Version' in r}
    for path in sorted(base.glob('*.gmt')):
        records={}
        for line in path.read_text().splitlines():
            parts=line.split('\t')
            if len(parts)>=3:records[parts[0]]=set(g.split(',')[0] for g in parts[2:] if g)
        libs[path.stem]=records;record=by_version.get(path.stem,{})
        sha=file_hash(path)
        if record.get('SHA256') and record['SHA256']!=sha:raise ValueError(f'GMT differs from download provenance: {path}')
        provenance.append(dict(library=path.stem,path=str(path),sha256=sha,version=path.stem,source=record.get('Source','local_file_source_unverified'),downloaded=record.get('Downloaded_UTC',''),n_terms=len(records)))
    write_csv(ctx.root/'enrichment/gene_set_provenance.csv',pd.DataFrame(provenance))
    return libs


def ora(query,universe,sets,scope,library,min_size=1,max_size=100000):
    universe=set(universe);query=set(query)&universe;N=len(universe);n=len(query);rows=[]
    if not N or not n:return rows
    for name,genes in sets.items():
        eligible=set(genes)&universe;M=len(eligible)
        if M<min_size or M>max_size:continue
        overlap=query&eligible;k=len(overlap)
        # Include zero-overlap terms in multiplicity correction; do not preselect significant hits.
        rows.append(dict(scope=scope,library=library,term=name,source_label='provided_unverified' if library=='provided_aging_sets' else 'versioned_local_GMT',
                         p_value=float(hypergeom.sf(k-1,N,M,n)),overlap_count=k,query_size=n,set_size=M,universe_size=N,
                         GeneRatio=f'{k}/{n}',genes=';'.join(sorted(overlap))))
    if rows:
        fdr=bh([r['p_value'] for r in rows])
        for r,q in zip(rows,fdr):r['FDR']=q
    return rows


def partial_association(x,age,meta):
    n=len(age)
    if n<8:return dict(n=n,beta_age=np.nan,partial_r=np.nan,p_value=np.nan,reason='fewer_than_8_observed_samples')
    covars=[];used=[]
    for c in ['tissue','GSE','sex','disease_global']:
        if c in meta and meta[c].nunique()>1:
            covars.append(pd.get_dummies(meta[c].fillna('missing').astype(str),drop_first=True,dtype=float).to_numpy());used.append(c)
    C=np.column_stack([np.ones(n)]+covars);rank=np.linalg.matrix_rank(C);df=n-rank-1
    if df<3:return dict(n=n,beta_age=np.nan,partial_r=np.nan,p_value=np.nan,reason='insufficient_residual_degrees_of_freedom',covariates=';'.join(used))
    xr=x-C@np.linalg.lstsq(C,x,rcond=None)[0];ar=age-C@np.linalg.lstsq(C,age,rcond=None)[0]
    if np.std(xr)<1e-10 or np.std(ar)<1e-10:return dict(n=n,beta_age=np.nan,partial_r=np.nan,p_value=np.nan,reason='constant_after_covariate_adjustment',covariates=';'.join(used))
    r=float(np.corrcoef(xr,ar)[0,1]);stat=r*np.sqrt(df/max(1-r*r,1e-15))
    return dict(n=n,beta_age=float(np.dot(ar,xr)/np.dot(ar,ar)),partial_r=r,p_value=float(2*student_t.sf(abs(stat),df)),
                residual_df=df,covariates=';'.join(used),reason='')


def associations(ctx,stable,ds):
    names=ds['schema']['features'];lookup={n:i for i,n in enumerate(names)};ctmap={c:i for i,c in enumerate(ds['schema']['celltypes'])};ix=ds['bag_index']
    rows=[];directions=[];train_sign={}
    for r in stable.itertuples():
        j=lookup[r.feature];obsrows=np.unique(ix['sample'][ix['celltype']==ctmap[r.Celltype]])
        meta=ds['meta'].iloc[obsrows].reset_index(drop=True);x=ds['X'][obsrows,j].toarray().ravel()
        stat=partial_association(x,meta.Age.to_numpy(float),meta);rows.append(dict(feature=r.feature,Celltype=r.Celltype,Gene=r.Gene,**stat));train_sign[r.feature]=np.sign(stat.get('beta_age',np.nan))
    frame=pd.DataFrame(rows,columns=['feature','Celltype','Gene','n','beta_age','partial_r','p_value','reason','covariates','residual_df'])
    frame['FDR']=bh(frame.p_value);write_csv(ctx.root/'enrichment/stable_gene_age_association.csv',frame)
    for it in json.loads((ctx.root/'audit/inventory.json').read_text())[1:]:
        ext=load_dataset(ctx.root/'features'/it['dataset']);m=ext['meta'];bi=ext['bag_index']
        for r in stable.itertuples():
            observed=np.unique(bi['sample'][bi['celltype']==ctmap[r.Celltype]])
            for gse in sorted(m.GSE.unique()):
                selected=observed[m.GSE.iloc[observed].to_numpy()==gse]
                mm=m.iloc[selected].reset_index(drop=True);xx=ext['X'][selected,lookup[r.feature]].toarray().ravel()
                stat=partial_association(xx,mm.Age.to_numpy(float),mm)
                sign=np.sign(stat.get('beta_age',np.nan));ts=train_sign[r.feature]
                directions.append(dict(feature=r.feature,Celltype=r.Celltype,Gene=r.Gene,dataset=it['dataset'],GSE=gse,method=it['method'],
                    source_kind=it['role'],n=len(selected),**{k:v for k,v in stat.items() if k!='n'},training_sign=ts,
                    consistent_direction=bool(sign==ts) if np.isfinite(sign) and np.isfinite(ts) else np.nan,
                    evidence_note='A/B and integration representations may overlap; not independent pooled replicates'))
    direction=pd.DataFrame(directions,columns=['feature','Celltype','Gene','dataset','GSE','method','source_kind','n','beta_age','partial_r','p_value','reason','training_sign','consistent_direction','evidence_note'])
    direction['FDR']=bh(direction.p_value);write_csv(ctx.root/'enrichment/independent_gse_direction.csv',direction)
    return frame,direction


def enrichment(ctx):
    table=pd.read_csv(ctx.root/'stability/full_train/bootstrap_all_features.csv');stable=pd.read_csv(ctx.root/'stability/stable_Celltype_Gene.csv');sets,limitations=qc_sets(ctx)
    approved,aliases=hgnc_index(ctx.cfg['enrichment']['hgnc_tsv'])
    gene_qc=[];mapping={}
    for gene in sorted(table.Gene.unique()):
        if approved is None:mapped=gene;status='unverified_HGNC_unavailable'
        elif gene in approved:mapped=gene;status='approved_exact'
        elif gene in aliases and len(aliases[gene])==1:mapped=next(iter(aliases[gene]));status='unique_HGNC_alias'
        else:mapped=None;status='unmapped_or_ambiguous_excluded_from_enrichment'
        mapping[gene]=mapped;gene_qc.append(dict(original_gene=gene,mapped_symbol=mapped or '',status=status))
    write_csv(ctx.root/'enrichment/feature_gene_hgnc_qc.csv',pd.DataFrame(gene_qc))
    for frame in [table,stable]:
        frame['OriginalGene']=frame.Gene
        frame['Gene']=frame.Gene.map(mapping)
    libraries=load_libraries(ctx)
    if not libraries:limitations.append('No local versioned GMT libraries: pathway enrichment not completed.')
    extra=ctx.cfg['enrichment']['stromal_sasp_evidence_csv']
    if extra and Path(extra).is_file():
        evidence=pd.read_csv(extra)
        required={'symbol','cell_context','reference','evidence','category'}
        if not required.issubset(evidence):raise ValueError(f'Stromal SASP evidence needs columns {required}')
        if evidence[list(required)].isna().any().any():raise ValueError('Incomplete stromal SASP evidence')
        sets['SASP_stromal_literature_curated']=set(evidence.symbol);write_csv(ctx.root/'enrichment/SASP_stromal_literature_curated_evidence.csv',evidence)
    else:limitations.append('No verified stromal SASP extension evidence supplied: no genes or citations were fabricated.')
    pathway=[];overlap=[];universe_rows=[]
    for scope in ['global']+sorted(table.Celltype.unique()):
        eligible=table if scope=='global' else table[table.Celltype==scope]
        # Have an opportunity to enter: expression variance eligible OR exposed in a bootstrap.
        eligible=eligible[(eligible.candidate.astype(bool))|(eligible.exposure_count>0)]
        query=stable if scope=='global' else stable[stable.Celltype==scope]
        universe=set(eligible.Gene.dropna());genes=set(query.Gene.dropna())
        universe_rows.extend(dict(scope=scope,Gene=g) for g in sorted(universe))
        overlap.extend(ora(genes,universe,sets,scope,'provided_aging_sets'))
        for lib,terms in libraries.items():
            pathway.extend(ora(genes,universe,terms,scope,lib,ctx.cfg['enrichment']['min_set_size'],ctx.cfg['enrichment']['max_set_size']))
    columns=['scope','library','term','source_label','p_value','FDR','overlap_count','query_size','set_size','universe_size','GeneRatio','genes']
    enrich=pd.DataFrame(pathway,columns=columns);ov=pd.DataFrame(overlap,columns=columns)
    write_csv(ctx.root/'enrichment/pathway_enrichment.csv',enrich);write_csv(ctx.root/'enrichment/aging_sets_overlap.csv',ov);write_csv(ctx.root/'enrichment/analysis_universe.csv',pd.DataFrame(universe_rows,columns=['scope','Gene']))
    ds=load_dataset(ctx.root/'features/train');assoc,direction=associations(ctx,stable,ds)
    if not len(stable):limitations.append('No features reach the prespecified stability threshold; no stable-gene enrichment or mechanism can be inferred.')
    if len(direction) and not direction.beta_age.notna().any():limitations.append('No external GSE has sufficient adjusted association degrees of freedom; cross-GSE direction verification is incomplete.')
    write_json(ctx.root/'enrichment/limitations.json',limitations)
    significant=enrich[(enrich.FDR<.05)&(enrich.overlap_count>0)].sort_values('FDR') if len(enrich) else enrich
    text=f'''# Biological interpretation: {ctx.lineage}

Stable EN pairs: {len(stable)}; unique genes: {stable.Gene.nunique()}; observed subtypes: {stable.Celltype.nunique()}.
These are predictive associations, not established causal aging drivers. If the internally selected model is not EN, this analysis is explicitly exploratory.
Universe is the lineage/subtype genes with training expression-variance eligibility or observed bootstrap topK exposure, recorded in analysis_universe.csv.
BH FDR is calculated over all size-eligible terms, including zero-overlap terms, within each scope/library. Cross-scope exploratory testing is not a confirmatory experiment.
All eight supplied aging lists retain provided_unverified provenance; QC does not convert them into official complete signatures.

## Observed pathway evidence

'''
    if len(significant):
        for r in significant.head(20).itertuples():text+=f'- {r.scope}: {r.library} / {r.term}; FDR={r.FDR:.3g}; genes={r.genes}; overlap {r.GeneRatio}; background n={r.universe_size}. Hypothesis: these measured predictive genes may reflect this annotated process; causal evidence and targeted literature review are required.\n'
    else:text+='No size-eligible pathway has BH FDR <0.05, or no stable query genes are available. No mechanism is asserted.\n'
    text+=f'\n## Confounding and replication\n\nTraining adjusted associations: {len(assoc)} tested, {int((assoc.FDR<.05).sum())} with FDR<0.05. Per-GSE representations: {len(direction)} rows. Covariates, degrees of freedom and untestable reasons are saved for every feature. No single-GSE or duplicate-method row is treated as independent confirmation.\n'
    text+='\n## Limitations\n\n'+'\n'.join('- '+s for s in limitations)+'\n'
    write_text(ctx.root/'reports/biological_interpretation.md',text)
    return [ctx.root/'enrichment'/f for f in ['pathway_enrichment.csv','aging_sets_overlap.csv','provided_gene_sets_qc.csv','stable_gene_age_association.csv','independent_gse_direction.csv','limitations.json']]+[ctx.root/'reports/biological_interpretation.md']
