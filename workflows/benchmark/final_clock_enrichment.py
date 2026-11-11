"""Reuse audited ORA; add direction/shared/specific scopes with measured universes."""
from concurrent.futures import ThreadPoolExecutor
from scipy.stats import hypergeom
from tabulastromal.biology import bh
from final_clock_common import *

def ora(query,universe,sets,scope,library,min_size=1,max_size=100000):
    universe=set(universe);query=set(query)&universe;N=len(universe);n=len(query)
    if not n or not N:return pd.DataFrame()
    rows=[]
    for name,genes in sets.items():
        eligible=genes&universe;M=len(eligible)
        if not min_size<=M<=max_size:continue
        hits=query&eligible
        rows.append(dict(scope=scope,library=library,term=name,overlap_count=len(hits),query_size=n,set_size=M,
            universe_size=N,GeneRatio=f'{len(hits)}/{n}',genes=';'.join(sorted(hits)),
            source_label='provided_unverified' if library=='provided_aging_sets' else 'versioned_local_GMT'))
    frame=pd.DataFrame(rows)
    if len(frame):
        frame['p_value']=hypergeom.sf(frame.overlap_count.to_numpy()-1,N,frame.set_size.to_numpy(),n)
        frame['FDR']=bh(frame.p_value)
    return frame

def run(ctx,smoke=False):
    outputs=[];stable=ctx.csv('stability/stable_celltype_gene_features.csv')
    source=ctx.csv('enrichment/pathway_enrichment.csv');aging=ctx.csv('enrichment/aging_sets_overlap.csv')
    provenance=ctx.csv('enrichment/gene_set_provenance.csv');libs={};versions={};pins=[]
    for row in provenance.itertuples():
        assert file_hash(row.path)==row.sha256
        terms={}
        for line in Path(row.path).read_text().splitlines():
            q=line.split('\t')
            if len(q)>=3:terms[q[0]]=set(g.split(',')[0] for g in q[2:] if g)
        libs[row.library]=terms;versions[row.library]=row.sha256;pins.append(Path(row.path))
    provided=ctx.json('enrichment/provided_gene_sets_cleaned.json')
    required={'SASP','SenMayo','CellAge','SigRS','SenUp','AgingAtlas','ASIG','GenAge'}
    assert required<=set(provided)
    libs['provided_aging_sets']={k:set(v) for k,v in provided.items()}
    versions['provided_aging_sets']=file_hash(ctx.path('enrichment/provided_gene_sets_cleaned.json'))
    # Independently verify saved hypergeometric P and BH, including zero hits.
    for frame in [source,aging]:
        expected=hypergeom.sf(frame.overlap_count.to_numpy()-1,frame.universe_size.to_numpy(),frame.set_size.to_numpy(),frame.query_size.to_numpy())
        np.testing.assert_allclose(frame.p_value,expected,atol=1e-12,rtol=1e-8)
        for _,part in frame.groupby(['scope','library']):np.testing.assert_allclose(part.FDR,bh(part.p_value),atol=1e-12,rtol=1e-8)
    universe_frame=ctx.csv('enrichment/analysis_universe.csv');universes={k:set(q.Gene) for k,q in universe_frame.groupby('scope')}
    qc=ctx.csv('enrichment/feature_gene_hgnc_qc.csv');mapping=dict(zip(qc.original_gene,qc.mapped_symbol))
    stable['MappedGene']=stable.Gene.map(mapping).replace('',np.nan)
    own=set(stable.MappedGene.dropna())
    other_lineage='Endo' if ctx.lineage=='Fibro' else 'Fibro';other=Context(other_lineage)
    otherstable=other.csv('stability/stable_celltype_gene_features.csv');oq=other.csv('enrichment/feature_gene_hgnc_qc.csv')
    omap=dict(zip(oq.original_gene,oq.mapped_symbol));othergenes=set(otherstable.Gene.map(omap).replace('',np.nan).dropna())
    ou=other.csv('enrichment/analysis_universe.csv');otheruniverse=set(ou.loc[ou.scope=='global','Gene'])
    shared=own&othergenes;own_universe=universes['global']
    queries={
        'positive_coefficient_associations':(set(stable.loc[stable.coefficient_mean>0,'MappedGene'].dropna()),own_universe),
        'negative_coefficient_associations':(set(stable.loc[stable.coefficient_mean<0,'MappedGene'].dropna()),own_universe),
        'shared_stable_genes':(shared,own_universe&otheruniverse),
        'lineage_specific_stable_genes':(own-othergenes,own_universe)}
    sets_positive,sets_negative=queries['positive_coefficient_associations'][0],queries['negative_coefficient_associations'][0]
    queries['mixed_coefficient_associations']=(sets_positive&sets_negative,own_universe)
    arguments=[]
    for scope,(query,universe) in queries.items():
        for lib,terms in libs.items():
            provided_lib=lib=='provided_aging_sets'
            arguments.append((query,universe,terms,scope,lib,1 if provided_lib else ctx.cfg['enrichment']['min_set_size'],
                              100000 if provided_lib else ctx.cfg['enrichment']['max_set_size']))
    with ThreadPoolExecutor(max_workers=ctx.options['io_workers']) as pool:
        frames=list(pool.map(lambda a:ora(*a),arguments))
    extra=pd.concat([f for f in frames if len(f)],ignore_index=True)
    pathway=pd.concat([source,extra[extra.library!='provided_aging_sets']],ignore_index=True)
    overlap=pd.concat([aging,extra[extra.library=='provided_aging_sets']],ignore_index=True)
    for f in [pathway,overlap]:
        f['database_version']=f.library;f['database_sha256']=f.library.map(versions)
        f['FDR_family']='all size-eligible terms including zero overlap, within scope/library'
    def save(rel,frame):
        p=ctx.path(rel);write_csv(p,frame);outputs.append(p)
    save('enrichment/enrichment_GO.csv',pathway[pathway.library.str.startswith('GO_')])
    save('enrichment/enrichment_Reactome.csv',pathway[pathway.library.str.startswith('Reactome')])
    save('enrichment/enrichment_KEGG.csv',pathway[pathway.library.str.startswith('KEGG')])
    save('enrichment/final_clock_pathway_enrichment.csv',pathway)
    save('enrichment/aging_gene_set_overlap.csv',overlap)
    save('enrichment/final_clock_gene_set_provenance.csv',provenance)
    genes=[]
    # The stable-feature lists remain raw schema symbols. Enrichment uses the
    # explicit existing HGNC mapping; never silently parse or relabel features.
    for g in sorted(own|othergenes):
        fibro=g in (own if ctx.lineage=='Fibro' else othergenes);endo=g in (own if ctx.lineage=='Endo' else othergenes)
        genes.append(dict(Gene=g,Fibroblast_stable=fibro,Endothelial_stable=endo,
            pattern='shared' if fibro and endo else 'Fibroblast_specific' if fibro else 'Endothelial_specific',
            measurable_in_both=g in own_universe&otheruniverse,
            interpretation='stability under different sample sizes/opportunity, not proven biological specificity'))
    save('enrichment/shared_vs_specific_stable_genes.csv',pd.DataFrame(genes))
    save('enrichment/final_clock_gene_symbol_qc.csv',qc[qc.original_gene.isin(stable.Gene)])
    query_manifest=[]
    for scope,(query,universe) in queries.items():
        query_manifest.append(dict(scope=scope,n_query=len(query&universe),n_universe=len(universe),
            genes=';'.join(sorted(query&universe)),universe=';'.join(sorted(universe)),
            status='tested' if query&universe else 'no_query_genes_no_test'))
    save('enrichment/final_clock_enrichment_inputs_manifest.csv',pd.DataFrame(query_manifest))
    for filename in ['stable_gene_age_association.csv','independent_gse_direction.csv','provided_gene_sets_qc.csv']:
        save('enrichment/final_clock_'+filename,ctx.csv('enrichment/'+filename));pins.append(ctx.path('enrichment/'+filename))
    significant=pathway[(pathway.scope=='global')&(pathway.FDR<.05)&(pathway.overlap_count>0)].sort_values('FDR')
    sigaging=overlap[(overlap.scope=='global')&(overlap.FDR<.05)&(overlap.overlap_count>0)].sort_values('FDR')
    counts=pathway.assign(significant=pathway.FDR.lt(.05)&pathway.overlap_count.gt(0)).groupby(['scope','library']).agg(
        n_tests=('term','size'),n_significant=('significant','sum')).reset_index()
    save('enrichment/final_clock_enrichment_summary.csv',counts)
    text=[f'# Final clock biological interpretation: {ctx.lineage}',
        f'Stable raw Celltype_Gene={len(stable)}, Celltypes={stable.Celltype.nunique()}, raw unique Genes={stable.Gene.nunique()}; HGNC-mapped unique genes={len(own)}. Shared mapped genes={len(shared)}.',
        'Original global/subtype GO BP/MF/CC, Reactome, KEGG and aging-set tests were verified and reused. New coefficient-direction, mixed-direction, shared and lineage-specific tests were actually computed.',
        'Universes are the measured lineage/subtype genes eligible by variance or observed bootstrap exposure; shared uses the intersection of lineage universes. Gene symbols are mapped using the existing audited HGNC table.',
        'BH correction includes zero-overlap terms within each library/scope. All scopes and direction strata are exploratory; no causal pro-aging/anti-aging claim follows from EN signs. All eight user aging sets remain provided_unverified, not official validated signatures.',
        f'Global significant pathway terms (FDR<0.05): {len(significant)}. Global significant provided aging sets: {len(sigaging)}.',
        '## Measured top global pathway overlaps']
    if len(significant):
        for r in significant.head(20).itertuples():text.append(f'- {r.library} / {r.term}: FDR={r.FDR:.4g}; overlap={r.GeneRatio}; genes={r.genes}; universe={r.universe_size}.')
    else:text.append('No global pathway meets FDR<0.05. No enrichment is inferred from nonsignificant visual ranking.')
    text += ['## Aging sets']
    for r in overlap[overlap.scope=='global'].sort_values('term').itertuples():text.append(f'- {r.term}: FDR={r.FDR:.4g}, overlap={r.GeneRatio}, genes={r.genes}.')
    text += ['## Celltype-specific and covariate evidence',
        'Per-celltype queries use the corresponding stable genes, not celltype names as gene sets. Full subtype tests and their own candidate universes are retained in the enrichment CSVs.',
        'Existing actually computed age-expression partial associations and per-GSE direction tables are reused with tissue/GSE/sex/disease covariates where identifiable. Sample-size, missing covariates and residual-DF failures remain explicit. Primary final-clock prediction strata by GSE/tissue/sex/disease are in external/final_clock_prediction_covariate_sensitivity.csv.',
        'Combined/GSE and integration representations can overlap and are not independent replications. Verified donor IDs are unavailable. Differences in subtype coverage and candidate opportunities can explain apparent lineage specificity.',
        'No new literature evidence, DOI or PMID is asserted. Gene→pathway→celltype links in these tables are measured overlaps and hypotheses, not literature-validated mechanisms.']
    p=ctx.path('reports/final_clock_biological_interpretation.md');write_text(p,'\n\n'.join(text)+'\n');outputs.append(p)
    p=ctx.path('reports/final_clock_enrichment_validation.json');write_json(p,dict(status='passed',global_significant_pathways=len(significant),
        global_significant_aging_sets=len(sigaging),shared_mapped_stable_genes=len(shared),n_mapped_genes=len(own),
        libraries=sorted(libs),provided_aging_sets_status='provided_unverified',source_p_and_BH_verified=True,
        new_scopes=list(queries),background='lineage/subtype candidate or bootstrap-exposed genes; shared intersection'))
    outputs.append(p)
    return outputs+pins+[ctx.path('enrichment/analysis_universe.csv'),ctx.path('enrichment/feature_gene_hgnc_qc.csv'),
        ctx.path('enrichment/provided_gene_sets_cleaned.json'),ctx.path('enrichment/pathway_enrichment.csv'),ctx.path('enrichment/aging_sets_overlap.csv')]
