"""Final-clock figures and reports, distinct from historical best-fold FT30 graphs."""
import string,textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from final_clock_common import *
BLUE='#32688f';LIGHT='#a9c9df';GRAY='#8b969f'

def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,
        'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})

def read_data(ctx):
    def read(p):return pd.read_csv(ctx.path(p),low_memory=False)
    d=dict(ctx=ctx,decision=versions(ctx),cv=read('stability/final_clock_internal_cv_metrics.csv'),
        oof=read('stability/final_clock_internal_oof_predictions.csv'),external=read('external/final_clock_external_combined_predictions.csv'),
        protocols=read('external/final_clock_external_protocol_comparison.csv'),stable=read('stability/stable_celltype_gene_features.csv'),
        cells=read('stability/stable_Celltype.csv'),genes=read('stability/stable_genes_unique.csv'),
        pathways=read('enrichment/final_clock_pathway_enrichment.csv'),aging=read('enrichment/aging_gene_set_overlap.csv'),
        association=read('enrichment/final_clock_stable_gene_age_association.csv'),shared=read('enrichment/shared_vs_specific_stable_genes.csv'),
        counts=read('stability/stable_feature_counts_summary.csv'))
    freq=pd.read_csv(ctx.path('stability/bootstrap_feature_selection_frequency.csv'),usecols=['selection_frequency'])
    hist,edges=np.histogram(freq.selection_frequency,bins=np.linspace(0,1,26));d['hist']=(hist,edges)
    d['top_features']=d['stable'].assign(strength=d['stable'].coefficient_mean.abs()).sort_values(['strength','feature'],ascending=[False,True]).head(50)
    d['top_cells']=d['cells'].sort_values(['stable_fraction_of_candidates','raw_stable_features','Celltype'],ascending=[False,False,True]).head(50)
    d['top_genes']=d['genes'].sort_values(['max_selection_frequency','coefficient_absolute_sum','Gene'],ascending=[False,False,True]).head(50)
    return d

def note(ax,text):
    ax.text(.5,.5,text,ha='center',va='center',transform=ax.transAxes,wrap=True);ax.set_xticks([]);ax.set_yticks([])

def cv_panel(ax,d):
    a=d['cv'][d['cv'].model_version=='original_EN'].set_index('fold')
    b=d['cv'][d['cv'].model_version=='stable_EN_0.70'].set_index('fold')
    for fold in a.index:ax.plot([0,1],[a.loc[fold,'MAE'],b.loc[fold,'MAE']],color='#b9c0c5',lw=.8,zorder=1)
    for x,q,c in [(0,a,BLUE),(1,b,GRAY)]:
        ax.bar(x,q.MAE.mean(),.45,color=c,alpha=.3,zorder=0)
        ax.scatter(np.repeat(x,len(q)),q.MAE,s=30,color=c,zorder=3)
        ax.text(x,q.MAE.max()+.7,f'Mean={q.MAE.mean():.2f}\nR²={q.R2.mean():.3f}',ha='center',fontsize=9)
    ax.set_xticks([0,1],['Original EN','Stable EN (0.70)']);ax.set_ylabel('Outer Valid MAE (years)')
    ax.set_ylim(0,max(a.MAE.max(),b.MAE.max())*1.25)
    ax.set_xlabel(f'Locked K={d["ctx"].K}; paired fivefold, conditional on selected K',fontsize=8)

def bootstrap_panel(ax,d):
    hist,edges=d['hist'];ax.bar(edges[:-1],hist,width=np.diff(edges),align='edge',color=LIGHT)
    ax.set_yscale('log');ax.set_ylabel('Celltype_Gene count (log scale)');ax.set_xlabel('Selection frequency across 300 group bootstraps')
    ax.axvline(.70,color=BLUE,ls='--');ax.axvline(.75,color=GRAY,ls=':')
    f=d['counts'][d['counts'].scope=='full_train']
    counts={float(r.threshold):int(r.n_stable_Celltype_Gene) for r in f.itertuples()}
    ax.text(.97,.95,f'300/300 valid full-development bootstraps\n≥0.70: {counts[.7]} stable pairs\n≥0.75: {counts[.75]} stable pairs\n{d["stable"].Celltype.nunique()} celltypes; {d["stable"].Gene.nunique()} unique genes',
            transform=ax.transAxes,ha='right',va='top',fontsize=9)

def scatter_panel(ax,d):
    selected=d['decision']['selected_version'];q=d['oof'][d['oof'].model_version==selected]
    for label,y,p,color,marker in [('Internal OOF (fold-trained)',q.Age,q.prediction,BLUE,'o'),
                                  ('Final full-dev → external no-adapt',d['external'].age_true,d['external'].age_pred,GRAY,'^')]:
        s=score(y,p);ax.scatter(y,p,s=17,alpha=.65,color=color,marker=marker,
            label=f'{label}\nn={len(y)}; R²={s["R2"]:.3f}; MAE={s["MAE"]:.2f}; PCC={s["PCC"]:.3f}')
    vals=np.r_[q.Age,q.prediction,d['external'].age_true,d['external'].age_pred];lo,hi=vals.min()-2,vals.max()+2
    ax.plot([lo,hi],[lo,hi],'--',color='#adb4b9',lw=.8)
    ax.set_xlabel('Chronological age (years)');ax.set_ylabel('Predicted age (years)');ax.legend(fontsize=8,frameon=False,loc='upper left')
    ax.text(.98,.03,f'Final family: {selected}; K={d["ctx"].K}\nPrimary external: no target labels used',transform=ax.transAxes,ha='right',fontsize=8)

def horizontal_numbers(ax,values):
    for i,v in enumerate(values):ax.annotate(f'{v:.3g}',(v,i),xytext=(3 if v>=0 else -3,0),textcoords='offset points',ha='left' if v>=0 else 'right',va='center',fontsize=6)
    ax.margins(x=.18,y=.012)

def coefficient_panel(ax,d):
    s=d['top_features'].iloc[::-1]
    ax.barh(s.feature,s.coefficient_mean,color=np.where(s.coefficient_mean>=0,BLUE,GRAY));ax.axvline(0,color='.7',lw=.6)
    ax.tick_params(axis='y',labelsize=7);ax.set_xlabel('Mean standardized bootstrap EN coefficient',fontsize=8)
    horizontal_numbers(ax,s.coefficient_mean)

def celltype_panel(ax,d):
    s=d['top_cells'].iloc[::-1];ax.barh(s.Celltype,s.stable_fraction_of_candidates,color=LIGHT)
    horizontal_numbers(ax,s.stable_fraction_of_candidates);ax.tick_params(axis='y',labelsize=7)
    ax.set_xlabel('Stable / variance-eligible candidate pairs',fontsize=8)

def gene_panel(ax,d):
    s=d['top_genes'].iloc[::-1];ax.barh(s.Gene,s.max_selection_frequency,color=BLUE)
    horizontal_numbers(ax,s.max_selection_frequency);ax.tick_params(axis='y',labelsize=7)
    ax.set_xlabel('Maximum pair selection frequency',fontsize=8)

def protocol_panel(ax,d):
    q=d['protocols'].copy();labels=['No-adapt\nall100','No-adapt\nheldout70','FT30 refit\nheldout70','FT30 refit\nall100']
    assert len(q)==4
    ax.bar(np.arange(4),q.MAE,color=[BLUE,GRAY,LIGHT,'#6b99b9'])
    for i,r in enumerate(q.itertuples()):ax.text(i,r.MAE+.4,f'{r.MAE:.2f}\nR²={r.R2:.3f}\nn={r.n}',ha='center',fontsize=8)
    ax.set_xticks(range(4),labels,fontsize=8);ax.set_ylabel('MAE (years)');ax.set_ylim(0,max(q.MAE.max()*1.38,1))
    ax.set_xlabel('all100 includes fitted30%; heldout70 remains separate',fontsize=8)

def pathway_panel(ax,d,library=None):
    s=d['pathways'];s=s[(s.scope=='global')&(s.overlap_count>0)]
    if library:s=s[s.library==library]
    s=s.sort_values(['FDR','library','term']).head(12).iloc[::-1]
    if s.empty:note(ax,'No eligible pathway overlap was observed');return
    aliases={'GO_Biological_Process_2025':'GO BP','GO_Molecular_Function_2025':'GO MF','GO_Cellular_Component_2025':'GO CC',
             'Reactome_Pathways_2024':'Reactome','KEGG_2021_Human':'KEGG','MSigDB_Hallmark_2020':'Hallmark'}
    labels=[textwrap.fill(f'{aliases.get(r.library,r.library)}: {r.term}',width=32,break_long_words=False) for r in s.itertuples()]
    ax.scatter(-np.log10(np.maximum(s.FDR,1e-300)),np.arange(len(s)),s=15+6*s.overlap_count,color=BLUE,alpha=.75)
    ax.set_yticks(range(len(s)),labels,fontsize=6);ax.axvline(-np.log10(.05),color=GRAY,ls='--',lw=.7)
    ax.set_xlabel('−log10(BH FDR); dashed: 0.05',fontsize=8)

def aging_panel(ax,d):
    s=d['aging'];s=s[s.scope=='global'].sort_values('term');ax.barh(s.term,s.overlap_count,color=LIGHT)
    for i,r in enumerate(s.itertuples()):ax.annotate(f' FDR={r.FDR:.2g}',(r.overlap_count,i),fontsize=7,va='center')
    ax.set_xlim(0,max(float(s.overlap_count.max()),1)*1.6);ax.set_xlabel('Stable gene overlap\n(user sets: provided_unverified)',fontsize=8)

def association_panel(ax,d):
    s=d['association'].copy();s['partial_r']=pd.to_numeric(s.partial_r,errors='coerce')
    s=s[s.partial_r.notna()].assign(strength=lambda q:q.partial_r.abs()).nlargest(50,'strength').iloc[::-1]
    if s.empty:note(ax,'No adjusted age association was estimable\nSee per-feature residual-DF/missing-covariate reasons');return
    ax.barh(s.feature,s.partial_r,color=np.where(s.partial_r>=0,BLUE,GRAY));ax.axvline(0,color='.6',lw=.7)
    ax.tick_params(axis='y',labelsize=6);ax.set_xlabel('Age-expression partial r (source covariate-adjusted)',fontsize=8);ax.margins(y=.012)

DRAW={'cv':cv_panel,'bootstrap':bootstrap_panel,'scatter':scatter_panel,'coefficients':coefficient_panel,
      'celltypes':celltype_panel,'genes':gene_panel,'protocol':protocol_panel,'pathways':pathway_panel,
      'aging':aging_panel,'association':association_panel}
TITLES={'cv':'Original vs stable EN: internal validation','bootstrap':'Group-bootstrap stability',
        'scatter':'Final clock: OOF and primary external no-adapt','coefficients':'Top 50 stable features (explanatory)',
        'celltypes':'Top 50 candidate-normalized celltypes','genes':'Top 50 stable genes',
        'protocol':'Final clock: distinct external evaluation sets','pathways':'Measured pathway overlaps',
        'aging':'Eight user aging sets (unverified)','association':'Covariate-adjusted age associations'}

def save(fig,stem,dpi):
    paths=[]
    for ext in ['png','pdf']:
        p=Path(str(stem)+'.'+ext);fig.savefig(p,dpi=dpi,facecolor='white');paths.append(p)
    plt.close(fig);return paths

def run(ctx,smoke=False):
    style();d=read_data(ctx);outputs=[];records=[]
    source_names=['models/model_versions.json','stability/final_clock_internal_cv_metrics.csv','stability/final_clock_internal_oof_predictions.csv',
        'stability/bootstrap_feature_selection_frequency.csv','stability/stable_feature_counts_summary.csv','stability/stable_celltype_gene_features.csv',
        'stability/stable_Celltype.csv','stability/stable_genes_unique.csv','external/final_clock_external_combined_predictions.csv',
        'external/final_clock_external_protocol_comparison.csv','enrichment/final_clock_pathway_enrichment.csv','enrichment/aging_gene_set_overlap.csv',
        'enrichment/final_clock_stable_gene_age_association.csv','enrichment/shared_vs_specific_stable_genes.csv']
    sources=[dict(path=str(ctx.path(p).resolve()),sha256=file_hash(ctx.path(p))) for p in source_names]
    for kind,stem in [('top_features','features'),('top_cells','celltypes'),('top_genes','genes')]:
        p=ctx.path(f'figures/final_clock_top50_{stem}.csv');write_csv(p,d[kind]);outputs.append(p)
    specs=[('cv','internal_cv',(9,5)),('bootstrap','bootstrap_selection_frequency',(9,5)),('scatter','external_combined_noadapt',(9,7)),
           ('coefficients','stable_features',(11,9)),('celltypes','stable_celltypes',(11,8)),('genes','stable_genes',(9,9)),
           ('protocol','external_ft30_comparison',(11,6)),('pathways','pathway_overlaps',(12,8)),('aging','aging_sets_overlap',(10,6)),
           ('association','gene_age_associations',(12,10))]
    for kind,name,size in specs:
        fig,ax=plt.subplots(figsize=size,layout='constrained');DRAW[kind](ax,d);ax.set_title(f'{ctx.lineage}: {TITLES[kind]}',loc='left')
        paths=save(fig,ctx.path('figures/final_clock_'+name),ctx.options['dpi']);outputs+=paths
        records += [dict(figure=str(p.resolve()),panels=kind,model_id=d['decision']['versions'][d['decision']['selected_version']]['model_id'],
                    source_files_json=json.dumps(sources),scope='primary final full-development EN; OOF conditional internal CV; FT30 separately identified') for p in paths]
    # Per-GSE representations are never pooled as independent validation cohorts.
    gse=pd.read_csv(ctx.path('external/final_clock_per_gse_metrics.csv'))
    fig,axes=plt.subplots(1,2,figsize=(16,6),layout='constrained')
    for ax,metric in zip(axes,['MAE','R2']):
        q=gse[gse.method==ctx.cfg['external']['primary_method']].sort_values('dataset')
        ax.bar(range(len(q)),q[metric],color=LIGHT);ax.set_xticks(range(len(q)),q.dataset,rotation=65,ha='right',fontsize=6)
        ax.set_ylabel(metric);ax.axhline(0,color='.7',lw=.5)
    fig.suptitle(f'{ctx.lineage}: final clock per-GSE, no adaptation; combined overlap flagged in CSV')
    paths=save(fig,ctx.path('figures/final_clock_per_gse'),ctx.options['dpi']);outputs+=paths
    # All requested functional libraries have explicit detailed panels.
    libraries=sorted(d['pathways'].library.unique());fig,axes=plt.subplots(3,2,figsize=(20,19),layout='constrained')
    for ax,lib in zip(axes.ravel(),libraries):pathway_panel(ax,d,lib);ax.set_title(lib,loc='left',fontsize=10)
    for ax in axes.ravel()[len(libraries):]:ax.remove()
    paths=save(fig,ctx.path('figures/final_clock_GO_Reactome_KEGG_enrichment'),ctx.options['dpi']);outputs+=paths
    # Celltype × pathway heatmap, using actual tests at each celltype universe.
    ct=d['stable'].Celltype.unique();q=d['pathways'][d['pathways'].scope.isin(ct)].copy()
    q['term_id']=q.library+': '+q.term
    chosen=q.sort_values('FDR').term_id.drop_duplicates().head(12)
    heat=q[q.term_id.isin(chosen)].pivot_table(index='scope',columns='term_id',values='FDR',aggfunc='first').reindex(columns=chosen)
    fig,ax=plt.subplots(figsize=(16,5),layout='constrained')
    if len(heat):
        image=ax.pcolormesh(-np.log10(np.maximum(heat.to_numpy(float),1e-300)),cmap='Blues',vmin=0,shading='flat')
        ax.invert_yaxis()
        ax.set_yticks(np.arange(len(heat))+.5,heat.index,fontsize=8)
        ax.set_xticks(np.arange(len(heat.columns))+.5,[textwrap.shorten(x,width=55,placeholder='…') for x in heat.columns],rotation=65,ha='right',fontsize=7)
        colorbar=fig.colorbar(image,ax=ax,label='−log10(BH FDR)')
        if colorbar.solids is not None:colorbar.solids.set_rasterized(False)
    else:note(ax,'No celltype-specific tests available')
    paths=save(fig,ctx.path('figures/final_clock_subtype_pathways'),ctx.options['dpi']);outputs+=paths
    # A detailed lineage Fig3, all panels grounded in the same manifest.
    keys=list(DRAW);fig,axes=plt.subplots(5,2,figsize=(22,34),layout='constrained',gridspec_kw={'height_ratios':[3.5,8,8,5.5,7]})
    for i,(ax,kind) in enumerate(zip(axes.ravel(),keys)):
        DRAW[kind](ax,d);ax.set_title(f'{string.ascii_uppercase[i]}. {TITLES[kind]}',loc='left',fontsize=11)
    fig.suptitle(f'{ctx.lineage}: Final ElasticNet chronological aging clock\nFull-development predictive clock; stable-feature interpretation is separately identified',fontsize=19)
    paths=save(fig,ctx.path(f'figures/Fig3_TabulaStromal_{ctx.lineage}_FINAL'),ctx.options['dpi']);outputs+=paths
    # Ensure every additional detailed/composite image has source/model records.
    recorded={r['figure'] for r in records}
    for p in outputs:
        if p.suffix in ['.png','.pdf'] and str(p.resolve()) not in recorded:
            records.append(dict(figure=str(p.resolve()),panels=p.stem,model_id=d['decision']['versions'][d['decision']['selected_version']]['model_id'],
                                source_files_json=json.dumps(sources),scope='final-clock analysis; exact evaluation scope in source CSV'))
    p=ctx.path('figures/final_clock_figure_data_manifest.csv');write_csv(p,pd.DataFrame(records));outputs.append(p)
    p=ctx.path('figures/final_clock_figure_manifest.json');write_json(p,dict(status='complete',sources=sources,
        figures=[str(p.resolve()) for p in outputs if p.suffix in ['.png','.pdf']],final_model_id=d['decision']['versions'][d['decision']['selected_version']]['model_id'],
        model_path=d['decision']['final_clock_path'],external_primary='no adaptation',internal_scope='out-of-fold models; conditional selected K',
        interpretation_model='full-development bootstrap stable EN; distinct from final predictive model when original wins'))
    outputs.append(p)
    return outputs+[Path(x['path']) for x in sources]

def combine():
    style();contexts={l:Context(l) for l in ['Fibro','Endo']}
    out=ROOT/'Stromal_Combined_Figure/outputs'
    sources=[dict(path=str(c.path('figures/final_clock_figure_manifest.json')),sha256=file_hash(c.path('figures/final_clock_figure_manifest.json'))) for c in contexts.values()]
    fingerprint=digest(dict(sources=sources,renderer=file_hash(__file__),dpi=350))
    marker=out/'FinalClock_FINAL_v2.completion.json'
    if marker.exists():
        previous=json.loads(marker.read_text())
        if previous.get('fingerprint')==fingerprint and all(Path(x['path']).exists() and file_hash(x['path'])==x['sha256'] for x in previous['artifacts']):
            LOG.info('Validated existing final combined figure; no redraw needed');return
    data={l:read_data(c) for l,c in contexts.items()}
    groups=[['cv'],['bootstrap'],['scatter'],['coefficients'],['celltypes','genes'],['protocol'],['pathways','aging']]
    heights=[3.5,3.5,4.2,8,8,4.7,6];fig=plt.figure(figsize=(22,sum(heights)+1),layout='constrained')
    grid=fig.add_gridspec(7,1,height_ratios=heights);records=[]
    for row,kinds in enumerate(groups):
        panels=[(k,l) for k in kinds for l in contexts];sub=grid[row].subgridspec(1,len(panels))
        for col,(kind,l) in enumerate(panels):
            ax=fig.add_subplot(sub[0,col]);DRAW[kind](ax,data[l]);letter=string.ascii_uppercase[len(records)]
            title=f'{letter}. {"Fibroblast" if l=="Fibro" else "Endothelial"} — {TITLES[kind]}'
            if kind in ['coefficients','genes']:title+=f' (n={len(data[l]["top_features" if kind=="coefficients" else "top_genes"])})'
            ax.set_title(textwrap.fill(title,width=50 if len(panels)==2 else 28),loc='left',fontsize=10 if len(panels)==2 else 9,pad=10)
            records.append(dict(panel=letter,lineage=l,kind=kind,row=row+1,column=col+1,
                model_id=data[l]['decision']['versions'][data[l]['decision']['selected_version']]['model_id']))
    fig.suptitle('Stromal final ElasticNet chronological aging clocks\nPrimary external: full-development model, NO adaptation | FT30 heldout70 and all100 separately labeled\nBootstrap-stable EN features are explanatory when the original EN predicts better',fontsize=18)
    paths=save(fig,out/'Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2',350)
    write_csv(out/'FinalClock_v2_panel_sources.csv',pd.DataFrame(records))
    shared=data['Fibro']['shared'];counts=shared.pattern.value_counts().reindex(['Fibroblast_specific','shared','Endothelial_specific'],fill_value=0)
    fig,axes=plt.subplots(2,1,figsize=(10,7),layout='constrained',gridspec_kw={'height_ratios':[3,1]})
    axes[0].bar(range(3),counts,color=[BLUE,LIGHT,GRAY]);axes[0].set_ylabel('Unique mapped stable genes')
    for i,n in enumerate(counts):axes[0].text(i,n+.5,str(n),ha='center')
    axes[0].set_xticks(range(3),['Fibroblast only','Shared','Endothelial only'])
    membership=np.array([[1,1,0],[0,1,1]])
    for i in range(2):
        for j in range(3):axes[1].scatter(j,i,s=90,color=BLUE if membership[i,j] else '#dedede')
    axes[1].set_yticks([0,1],['Fibroblast','Endothelial']);axes[1].set_xticks([]);axes[1].invert_yaxis()
    fig.suptitle('Stable-gene overlap: measured opportunity/power differs by lineage')
    paths+=save(fig,out/'Stromal_FinalClock_shared_specific_genes',350)
    write_csv(out/'FinalClock_shared_vs_specific_stable_genes.csv',shared)
    write_json(marker,dict(status='complete',fingerprint=fingerprint,created_at=now(),renderer_sha256=file_hash(__file__),sources=sources,
        panel_sources=records,artifacts=[dict(path=str(p.resolve()),sha256=file_hash(p),size=p.stat().st_size) for p in paths]))
    LOG.info('Final combined clocks saved: %s',out)

def write_final_report():
    rows=[];text=['# FINAL CLOCK LOCATION — actual completed Phase3/4','Updated: '+now(),
        'All results use the existing lineage directories. Historical 13-model benchmarks and best-fold FT30 A–F retain their original identities. Final-v2 is the new full-development EN analysis.']
    for l in ['Fibro','Endo']:
        c=Context(l);v=versions(c);s=c.csv('stability/stable_celltype_gene_features.csv');cells=c.csv('stability/stable_celltypes_ranked.csv');genes=c.csv('stability/stable_genes_unique.csv')
        no=c.csv('external/final_clock_external_combined_metrics.csv').iloc[0];ft=c.csv('external/external_ft30_all_100_metrics.csv').iloc[0]
        held=c.csv('external/external_ft30_heldout70_metrics.csv');h=held[held.model_version=='final_clock_FT30_refit'].iloc[0]
        enr=c.json('reports/final_clock_enrichment_validation.json');params=v['full_train_hyperparameters'];chosen=v['versions'][v['selected_version']]
        row=dict(lineage=l,K=c.K,alpha=params['alpha'],l1_ratio=params['l1_ratio'],n_bootstrap=300,threshold=.7,
            n_stable_Celltype_Gene=len(s),n_stable_Celltypes=len(cells),n_unique_stable_Genes=len(genes),
            n_development_samples=chosen['training_samples'],n_development_groups=chosen['training_groups'],
            internal_original_R2=v['original_mean_R2'],internal_stable_R2=v['stable_mean_R2'],
            internal_original_MAE=v['original_mean_MAE'],internal_stable_MAE=v['stable_mean_MAE'],selected_version=v['selected_version'],
            final_model_path=v['final_clock_path'],n_external=int(no['n']),noadapt_R2=no.R2,noadapt_MAE=no.MAE,noadapt_MedAE=no.MedAE,
            FT30_all100_R2=ft.R2,FT30_all100_MAE=ft.MAE,FT30_all100_MedAE=ft.MedAE,
            FT30_heldout70_R2=h.R2,FT30_heldout70_MAE=h.MAE,FT30_heldout70_MedAE=h.MedAE,
            n_significant_global_pathways=enr['global_significant_pathways'],n_significant_global_aging_sets=enr['global_significant_aging_sets'])
        rows.append(row);write_csv(c.path('reports/final_clock_results_summary.csv'),pd.DataFrame([row]))
        text += [f'## {l}',f'Final predictive version: **{v["selected_version"]}**; K={c.K}, alpha={params["alpha"]}, l1_ratio={params["l1_ratio"]}.',
            'Final model: `'+v['final_clock_path']+'`',
            'Original full-development model: `'+v['versions']['original_EN']['path']+'`',
            'Stable full-development model: `'+v['versions']['stable_EN']['path']+'`',
            'Best-fold reference (not final): `'+str(c.path('models/EN_best_fold_reference.joblib'))+'`',
            'Schema: `'+str(c.path('models/final_feature_schema.json'))+'`; training manifest: `'+str(c.path('models/final_training_manifest.json'))+'`.',
            f'Development={chosen["training_samples"]} samples/{chosen["training_groups"]} GSE/GSM groups, no verified donor IDs. Bootstrap: 300 valid fits per scope (5 outer TRAIN + full development), 0.70 primary, 0.75 sensitivity.',
            f'Counts: stable Celltype_Gene={len(s)}, Celltypes={len(cells)}, unique Genes={len(genes)}.',
            f'Conditional fivefold R2 original/stable={v["original_mean_R2"]:.6f}/{v["stable_mean_R2"]:.6f}; MAE={v["original_mean_MAE"]:.6f}/{v["stable_mean_MAE"]:.6f}. Stable features are explanatory because predictive accuracy worsened; external did not choose the model.',
            f'Primary no-adapt n={int(no["n"])}, R2={no.R2:.6f}, MAE={no.MAE:.6f}, MedAE={no.MedAE:.6f}.',
            f'FT30 all100 R2={ft.R2:.6f}, MAE={ft.MAE:.6f}, MedAE={ft.MedAE:.6f}; includes fitted30%, mixed apparent/unseen.',
            f'FT30 heldout70 R2={h.R2:.6f}, MAE={h.MAE:.6f}, MedAE={h.MedAE:.6f}. Paired no-adapt70 is in external_ft30_heldout70_metrics.csv.',
            'Complete CSVs: '+', '.join('`'+str(c.path('stability/'+p))+'`' for p in ['stable_celltype_gene_features.csv','stable_celltypes_ranked.csv','stable_genes_unique.csv','stable_genes_by_celltype.csv']),
            'Per-sample predictions: `'+str(c.path('external/final_clock_external_combined_predictions.csv'))+'`; per GSE: `'+str(c.path('external/final_clock_per_gse_predictions.csv'))+'`.',
            'Enrichment: `'+str(c.path('enrichment/enrichment_GO.csv'))+'`, enrichment_Reactome.csv, enrichment_KEGG.csv, aging_gene_set_overlap.csv, shared_vs_specific_stable_genes.csv.',
            f'Global significant pathway terms={enr["global_significant_pathways"]}; significant provided_unverified aging sets={enr["global_significant_aging_sets"]}.',
            'Detailed Fig3: `'+str(c.path(f'figures/Fig3_TabulaStromal_{l}_FINAL.png'))+'` and same stem.pdf.',
            '### Complete stable Celltype list','; '.join(cells.Celltype),
            '### Complete unique stable Gene list','; '.join(genes.Gene),
            '### Complete stable Celltype_Gene list','; '.join(s.feature)]
    table=pd.DataFrame(rows);write_csv(ROOT/'final_clock_results_summary.csv',table)
    write_csv(ROOT/'Stromal_Combined_Figure/outputs/FinalClock_results_summary.csv',table)
    text += ['## Actual executable loading/prediction',
        '```python\nfrom final_clock_common import load_final_clock, predict_age\nclock = load_final_clock("'+str(ROOTS['Fibro']/'models/final_clock.joblib')+'")\npredictions = predict_age("'+str(ROOTS['Fibro']/'features/combined')+'", clock)\nprint(predictions[["sample_id", "age_true", "age_pred", "age_acceleration_raw"]])\n# A new h5ad: predict_age("/absolute/new_data.h5ad", clock, expression="X")\n# GSE bbknn uses expression="raw" exactly as the original audited pipeline.\n```',
        'Missing subtype/gene features are zero-filled BEFORE the frozen training scaler, as explicitly recorded. Named raw pseudobulk DataFrames are reordered by the saved schema. Age units are years; raw residuals are not validated biological age acceleration.',
        'Primary external no-adapt is the final full-development clock. FT30 adapts a deepcopy of THIS final source using the reference fixed-HP adaptation-only regularized refit, not incremental continuation; preprocessing/features stay frozen. 70% and 100% have separate files and no mixed-metric scoring.',
        'Final combined figure: `'+str(ROOT/'Stromal_Combined_Figure/outputs/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.png')+'` and .pdf.',
        'Methods, commands and checks: README_FINAL_CLOCK.md; final_clock_pipeline.ipynb; reports/final_clock_validation.json; logs/final_clock_pipeline.log within each lineage root.']
    write_text(ROOT/'FINAL_CLOCK_LOCATION.md','\n\n'.join(text)+'\n')
    for l in ['Fibro','Endo']:write_text(ROOTS[l]/'reports/FINAL_CLOCK_LOCATION.md','\n\n'.join(text)+'\n')
    print(table.to_string(index=False),flush=True)
