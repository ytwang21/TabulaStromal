"""Cross-lineage A–L figure with panel source data and evidence-aware reports."""
from pathlib import Path
import string
import textwrap
import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import load_cells, read_json, save_table, stage_status, write_json
from .checkpointing import cached_frame, cached_action, save_figure

BLUE, TEAL, WARM = '#326b9e', '#318c8d', '#c66a43'
STYLE_VERSION = 'large_text_v1'
LABEL_SIZE, PATHWAY_SIZE = 12, 11
FIGURE_SIZE = (24, 50)
ROW_HEIGHTS = [.8, 1.5, 1.05, 1.5, 1., 2.]


def single_panel_size(theme, data):
    """Keep dense subtype/term labels legible in the standalone exports too."""
    if theme == 'pathways':
        return (14, 11)
    if theme == 'sid':
        n = data.celltype.nunique() if 'celltype' in data else 0
        return (14, max(7, 2 + .28*n))
    if theme in ['ranking', 'age', 'clock']:
        n = int(data.status.eq('ok').sum()) if 'status' in data else len(data)
        return (14, max(7, 2 + .28*n))
    return (14, 8)


def table(run, lineage, name):
    p = Path(run)/'tables'/f'{lineage}_{name}.csv'
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def unavailable(ax, reason):
    ax.text(.5, .5, reason, ha='center', va='center', wrap=True, transform=ax.transAxes, color='#66717b')
    ax.set_axis_off()


def embedding_source(cfg, run, lineage):
    cells = load_cells(run, lineage)
    n = min(cfg['embedding_max_cells'], len(cells))
    rng = np.random.default_rng(cfg['seed'])
    # Proportional sampling within GSE × donor × subtype; cap plotted points.
    sampled = []
    for _, ix in cells.groupby(['gse', 'sample', 'celltype'], observed=True, dropna=False).indices.items():
        k = min(len(ix), max(1, int(round(len(ix)*n/len(cells)))))
        sampled.extend(rng.choice(ix, k, replace=False))
    if len(sampled) > n:
        sampled = rng.choice(sampled, n, replace=False)
    x = cells.iloc[np.asarray(sampled, dtype=int)].copy()
    positions = x.atlas_row.to_numpy()
    order = np.argsort(positions)
    with h5py.File(cfg['lineages'][lineage]['atlas'], 'r') as f:
        if 'obsm/X_umap' not in f or not isinstance(f['obsm/X_umap'], h5py.Dataset):
            return pd.DataFrame()
        arr = np.asarray(f['obsm/X_umap'][positions[order], :2])
    x = x.iloc[order]
    source = x[['lineage', 'cell_key', 'donor_key', 'gse', 'celltype', 'celltype_general']].reset_index(drop=True)
    source['UMAP1'], source['UMAP2'] = arr[:, 0], arr[:, 1]
    source['sampling'] = 'GSE × donor × subtype stratified; capped for visualization; not abundance estimator'
    return source


def plot_panel(ax, theme, data, lineage):
    color = BLUE if lineage == 'fibro' else TEAL
    if data.empty:
        unavailable(ax, 'Analysis unavailable\nSee audit and stage diagnosis')
        return
    if theme == 'embedding':
        groups = sorted(data.celltype_general.dropna().unique())
        palette = plt.get_cmap('tab10')
        for i, c in enumerate(groups):
            d = data.loc[data.celltype_general.eq(c)]
            ax.scatter(d.UMAP1, d.UMAP2, s=.6, color=palette(i), label=c, rasterized=True, linewidths=0, alpha=.65)
        ax.set(xlabel='UMAP1', ylabel='UMAP2', title=f'Atlas subtype map (stratified n={len(data):,})')
        ax.legend(loc='best', markerscale=6, frameon=False, fontsize=LABEL_SIZE)
        ax.set_xticks([]); ax.set_yticks([])
    elif theme == 'ranking':
        d = data.sort_values('mean')
        y = np.arange(len(d))
        ax.errorbar(d['mean'], y, xerr=np.maximum(np.vstack([d['mean']-d.ci_low, d.ci_high-d['mean']]), 0), fmt='o', ms=4, color=color, elinewidth=1)
        ax.set_yticks(y, [f'{c} (n={int(n)})' for c, n in zip(d.celltype, d.n_donors)], fontsize=LABEL_SIZE)
        ax.set(xlabel='hUSI (native scale)\nEqual donor mean ± 95% bootstrap CI', title='Subtype hUSI; donor unit')
    elif theme in ['age', 'clock']:
        d = data.loc[data.status.eq('ok')].sort_values('effect')
        if d.empty:
            unavailable(ax, 'Adjusted effects not identifiable\nSee model status / support tables')
            return
        y = np.arange(len(d))
        for i, (_, r) in enumerate(d.iterrows()):
            ax.plot([r.ci_low, r.ci_high], [i, i], color='#9ba8b2', lw=1)
            ax.scatter(r.effect, i, s=16, color=WARM if r.q < .05 else color)
        ax.axvline(0, color='#aaa', lw=.6, ls='--')
        ax.set_yticks(y, [f'{c} (n={int(n)}, q={q:.2g})' for c, n, q in zip(d.celltype, d.n_donors, d.q)], fontsize=LABEL_SIZE)
        ax.set(xlabel='hUSI slope per year of chronological age' if theme == 'age' else 'hUSI effect per year of delta age, adjusted for age',
               title='Within-subtype adjusted age effects' if theme == 'age' else 'Held-out GeneModel age acceleration\nAge-adjusted association')
    elif theme == 'sid':
        mat = data.pivot(index='celltype', columns='recommendation', values='mean_recommendation')
        # Columnwise standardization is descriptive and does not sum SID values.
        sd = mat.std(axis=0).replace(0, np.nan)
        z = (mat-mat.mean(axis=0))/sd
        im = ax.imshow(np.ma.masked_invalid(z.to_numpy()), aspect='auto', cmap='BrBG', vmin=-2, vmax=2, interpolation='nearest')
        ax.set_yticks(np.arange(len(z)), z.index, fontsize=LABEL_SIZE)
        ax.set_xticks(np.arange(len(z.columns)), [c.replace('SenCID_rec_', '') for c in z.columns], fontsize=LABEL_SIZE)
        ax.set(title='SID recommendations\nNot a senescence score', xlabel='Each SID standardized across subtypes; donor means')
        bar = plt.colorbar(im, ax=ax, fraction=.03, pad=.02, label='Column z score')
        bar.ax.tick_params(labelsize=LABEL_SIZE)
    elif theme == 'pathways':
        # Display strongest terms per library, all tests remain in source table.
        d = data.loc[data.branch.eq('all_features') & data.direction.eq('all') & (data.overlap_n >= 2)].copy()
        if d.empty:
            unavailable(ax, 'No eligible pathway overlaps\nAll tested terms retained in source table')
            return
        d = d.sort_values('q').drop_duplicates(['library', 'term']).head(14).iloc[::-1]
        val = -np.log10(d.q.clip(lower=1e-300))
        labels = [textwrap.fill(f'{r.library}: {r.term}', width=43,
                               break_long_words=False, break_on_hyphens=False)
                  + f'\n{r.celltype}, {int(r.overlap_n)} genes' for _, r in d.iterrows()]
        # Allocate space by the actual wrapped line count, so long terms cannot
        # collide with the next label when fonts are enlarged.
        lines = np.asarray([label.count('\n') + 1 for label in labels])
        y = np.r_[0., np.cumsum((lines[:-1]+lines[1:])/2 + 1.4)]
        ax.scatter(val, y, s=20+7*d.overlap_n.clip(upper=20), color=color, alpha=.8)
        ax.set_yticks(y, labels, fontsize=PATHWAY_SIZE)
        ax.set_ylim(-lines[0]/2-1, y[-1]+lines[-1]/2+1)
        ax.axvline(-np.log10(.05), lw=.6, ls='--', color='#999')
        ax.set(xlabel='−log10(BH q); full lineage testing family', title='Stable-feature enrichment\nExploratory functional candidates')
    elif theme == 'feature_integration':
        d = data.dropna(subset=['effect', 'signed_coef_mean']) if 'effect' in data else pd.DataFrame()
        if d.empty:
            unavailable(ax, 'No identifiable subtype age slopes\nStable-feature table remains available')
            return
        ax.scatter(d.effect, d.signed_coef_mean, color=color, s=20)
        ax.axhline(0, color='#aaa', lw=.6); ax.axvline(0, color='#aaa', lw=.6)
        # Separate labels in data space and retain leader lines to their points.
        ordered = d.sort_values('signed_coef_mean')
        span = max(float(d.signed_coef_mean.max()-d.signed_coef_mean.min()), .01)
        gap = span * .075
        previous_y = -np.inf
        x_mid = (float(d.effect.min()) + float(d.effect.max()))/2
        for _, r in ordered.iterrows():
            label_y = max(float(r.signed_coef_mean), previous_y + gap)
            previous_y = label_y
            right = r.effect > x_mid
            ax.annotate(r.celltype, (r.effect, r.signed_coef_mean), xytext=(-8 if right else 8, label_y),
                        fontsize=LABEL_SIZE, textcoords=('offset points', 'data'), ha='right' if right else 'left', va='bottom',
                        arrowprops={'arrowstyle':'-', 'color':'#9ba8b2', 'lw':.5})
        ax.set_ylim(float(d.signed_coef_mean.min())-.08*span, previous_y+.18*span)
        ax.set(xlabel='Adjusted hUSI chronological-age slope', ylabel='Mean signed stable EN coefficient', title='Subtype continuous evidence comparison\nExploratory')
    ax.spines[['top', 'right']].set_visible(False)


def make_reports(cfg, run, lineages, source_run=None):
    run = Path(run)
    read_run = Path(source_run) if source_run else run
    findings = ['# Biological findings and evidence strength',
        'These are cross-sectional associations. hUSI is a senescence-associated index, not age in years; recSID is an identity recommendation. No lineage or subtype was preselected for a positive finding.']
    for lin in lineages:
        cohort = table(read_run, lin, 'cohort')
        if not cohort.empty:
            r = cohort.iloc[0]
            findings.append(f'## {lin}\n\nCoverage: {int(r.n_cells):,} Atlas cells, {int(r.donors):,} verified donors, {int(r.hUSI_cells):,} hUSI cells, {int(r.recSID_cells):,} recSID cells; common coverage {int(r.hUSI_recSID_common_cells):,}. Coverage differences are technical support, not biological senescence differences.')
        rank = table(read_run, lin, 'subtype_ranking')
        if not rank.empty:
            r = rank.loc[rank.score.eq('hUSI')].sort_values('mean', ascending=False).iloc[0]
            findings.append(f'Descriptive evidence: highest observed donor-mean hUSI is {r.celltype}: {r["mean"]:.4f}, 95% donor-bootstrap CI [{r.ci_low:.4f}, {r.ci_high:.4f}], n_donors={int(r.n_donors)}, n_GSE={int(r.n_gse)}. This ranking is not an age trend and may reflect cell identity, tissue or disease mixture.')
        age = table(read_run, lin, 'age_associations')
        if not age.empty:
            primary = age.loc[age.sensitivity.eq(f'min_cells_{cfg["min_cells"]}') & age.adjusted.eq(True) & age.status.eq('ok')]
            positives = primary.loc[(primary.effect > 0) & (primary.q < .05)]
            negatives = primary.loc[(primary.effect < 0) & (primary.q < .05)]
            findings.append(f'Adjusted association evidence: {len(primary)} identifiable primary subtype models; {len(positives)} positive and {len(negatives)} negative age slopes at BH q<0.05. Positive discovery does not by itself establish stable pro-aging status or pan-tissue replication.')
            for _, r in primary.sort_values('q').head(3).iterrows():
                healthy = age.loc[age.sensitivity.eq('healthy_only') & age.celltype.eq(r.celltype) & age.score.eq(r.score) & age.status.eq('ok')]
                robust = len(healthy) and np.sign(healthy.iloc[0].effect) == np.sign(r.effect) and healthy.iloc[0].q < .05
                findings.append(f'{r.celltype}: hUSI slope {r.effect:.6f}/year, 95% CI [{r.ci_low:.6f}, {r.ci_high:.6f}], q={r.q:.3g}; donors={int(r.n_donors)}, GSE={int(r.n_gse)}. Healthy-only adjusted direction/significance supported={bool(robust)}. Evidence remains association-level; alternatives include residual batch, age support imbalance and transcriptional identity.')
        clock = table(read_run, lin, 'clock_associations')
        if not clock.empty:
            ok = clock.loc[clock.status.eq('ok')]
            sig = ok.loc[ok.q < .05]
            findings.append(f'Clock integration: {len(ok)} identifiable age-adjusted delta-age models, {len(sig)} at q<0.05, using only completed held-out folds. Stable features derive from a separate earlier EN model. Partial folds limit precision and do not constitute independent external replication.')
        overlap = table(read_run, lin, 'score_feature_gene_overlap')
        if not overlap.empty:
            r = overlap.iloc[0]
            findings.append(f'Circularity check: {int(r.n_shared)}/{int(r.n_stable_genes)} stable genes overlap the {r.method} scoring gene list. Overlap reduces the independence of biological validation; non-overlap ORA is a limited sensitivity check, not a re-estimated independent score.')
    findings.append('## Cross-lineage interpretation\n\nCompare signs, uncertainty and support, not native hUSI absolute values: local hUSI computation minmax-normalizes each input dataset. SID landscape is not evidence for cellular senescence intensity. No causal mechanism, true cellular aging trajectory, or strong cross-lineage mechanistic evidence is established by these analyses. Pathway enrichments remain exploratory stable-feature hypotheses.')
    (run/'results/biological_findings.md').write_text('\n\n'.join(findings))
    (run/'results/statistical_methods.md').write_text('''# Statistical methods

Biological donor is sample, explicitly confirmed by the user; conflicting donor ages/sex are excluded. Cell ID joins are exact and checked against sample/subtype metadata, minimum score-row match rate 95%. Atlas coverage has a separate denominator.

Age is measured in years using Age; age-subset annotations are explicitly selected, lowercase age is audited, numerical cross-file conflicts excluded. Age models use verified members of the age subset with numerical age. No in-sample clock predictions enter association models.

Scores are summarized per donor × tissue × subtype with mean/median/IQR; primary minimum 20 scored cells, sensitivities 10/50. Ranking averages tissues within donor and bootstraps independent donors with fixed seed. Paired subtype contrasts compare donor/tissue matched subtype means to equal-weight means of other present subtypes; these are exploratory identity contrasts.

OLS adjusts for available tissue, sex, GSE, disease, treatment, platform and QC proxies. Complete cases are used; missing metadata are not imputed. Exactly nested nuisance factors are replaced by a full-rank basis of the same covariate span, preserving adjustment. Age/target collinearity or insufficient residual degrees of freedom blocks the model. Confidence intervals and p values use donor-cluster robust covariance with t reference. Unadjusted and adjusted estimates are separate. GSE is a fixed effect; inference still assumes residual independence across donors after adjustment.

Healthy-only adjusted fits and tissue-specific fits are reported when identifiable. Leave-one-GSE-out gives descriptive unadjusted directions, not a claim of independently replicated adjusted effects. Observed donor age quantiles determine descriptive age groups. Subtype composition includes genuine observed zero counts and uses all eligible cells as denominator, independently of score coverage.

BH-FDR families are all subtypes per lineage/score/model-sensitivity/adjustment, all subtype composition tests per lineage, and all tested pathway terms/subtype/direction/library/overlap branches within a lineage. Gene associations, if enabled, use one lineage-wide family. No feature/pathway is selected for reporting before adjustment.

GeneModel uses prespecified K8000 ElasticNet completed outer held-out folds. Train-only preprocessing is documented in the inspected existing training source; donor fold disjointness is verified from saved Train/Valid rows. Delta age = prediction − observed chronological age; regression adjusts for chronological age. Test-set age-bias correction is not fitted. Stable-feature interpretation is explicitly a separate older model branch.

Local GMT enrichment uses actual tested candidate Celltype_Gene features as subtype-specific background, including zero-hit terms in the testing family. The feature parser matches the longest verified actual subtype name prefix. Feature/signature overlap and feature-only non-overlap sensitivity are reported. Stable-feature coefficients are model-dependent and do not identify causal pro-aging genes.
''')
    (run/'results/limitations.md').write_text('''# Limitations and skipped analyses

Existing SenCID outputs are recSID identity recommendations. SID_Score/probability/binary predictions were not saved. Overall senescence intensity, high/high discordance categories and native binary inter-method agreement are blocked, rather than constructed from recommendation values.

hUSI covers an age-named subset, not the whole Atlas; numerical Age is missing for some subset cells. Full Atlas cell and composition context must not be mistaken for full Atlas score coverage. Scores are minmax scaled within the input dataset. SSE assumes two states; local GMM labels are multi-state. Neither supplies a transferable absolute threshold across cohorts/lineages.

The raw-named matrices contain sampled fractional values. Count-based DE is not run without independently verified count provenance; no rounding or pseudo-count reconstruction is performed. Optional pseudobulk marker/score associations are exploratory and do not substitute for a controlled count-based DE study. Canonical marker positivity/negativity alone cannot establish or exclude senescence. SASP/SenMayo score analysis is not asserted without verified signature inputs.

Sparse strata and nested GSE/tissue/age support can make adjusted effects unidentifiable. Complete-case covariate adjustment can change the study population; inspect model-specific n_donors and support. Primary contrasts use the scored cohort and may be biased by selective score availability. Donor age inconsistencies are excluded, which can remove pooled or otherwise ambiguous samples.

GeneModel benchmark files may still be growing: exact prediction snapshots and hashes are saved. Current completed folds may not cover all donors. Earlier stable-feature EN models differ from the rebenchmark clock; the association and interpretation branches must not be conflated. One library representation is used; alternate batch-correction files are not independent cohorts. Circularity from scoring-signature overlap and cross-sectional confounding limits mechanistic claims.
''')


def run(cfg, run, lineages, source_run=None):
    run = Path(run)
    read_run = Path(source_run) if source_run else run
    for folder in ['figures/panels', 'figures/source_data', 'results', 'checkpoints']:
        (run/folder).mkdir(parents=True, exist_ok=True)
    for lin in lineages:
        if not (read_run/'checkpoints'/f'{lin}_audit.json').exists():
            raise RuntimeError('Phase 0 gate not passed')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12, 'axes.titlesize': 16, 'axes.labelsize': 14,
                        'xtick.labelsize': 12, 'ytick.labelsize': LABEL_SIZE, 'axes.titlepad': 12,
                        'axes.labelpad': 9, 'ytick.major.pad': 6,
                        'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white'})
    # The six paired themes share definitions; unavailable clock replaces with
    # actual stable-feature integration, never fake zeros or empty data panels.
    sources = {}
    for lin in lineages:
        saved_embedding = list((read_run/'figures/source_data').glob(f'*_{lin}_embedding.csv'))
        sources[(lin, 'embedding')] = cached_frame(['embedding',lin],lambda:pd.read_csv(saved_embedding[0]) if saved_embedding else embedding_source(cfg,read_run,lin))
        rank = table(read_run, lin, 'subtype_ranking')
        sources[(lin, 'ranking')] = rank.loc[rank.score.eq('hUSI')] if not rank.empty else rank
        age = table(read_run, lin, 'age_associations')
        sources[(lin, 'age')] = age.loc[age.sensitivity.eq(f'min_cells_{cfg["min_cells"]}') & age.adjusted.eq(True) & age.score.eq('hUSI')] if not age.empty else age
        sources[(lin, 'sid')] = table(read_run, lin, 'SID_recommendation_landscape')
        sources[(lin, 'clock')] = table(read_run, lin, 'clock_associations')
        sources[(lin, 'feature_integration')] = table(read_run, lin, 'feature_subtype_integration')
        sources[(lin, 'pathways')] = table(read_run, lin, 'pathway_enrichment')
    fifth = 'clock' if all(not sources[(lin, 'clock')].empty and sources[(lin, 'clock')].status.eq('ok').any() for lin in lineages) else 'feature_integration'
    themes = ['embedding', 'ranking', 'age', 'sid', fifth, 'pathways']
    fig, axes = plt.subplots(6, len(lineages), figsize=(FIGURE_SIZE if len(lineages)==2 else (14, 50)),
                             gridspec_kw={'height_ratios': ROW_HEIGHTS}, squeeze=False, layout='constrained')
    fig.get_layout_engine().set(w_pad=.15, h_pad=.18, wspace=.12, hspace=.035)
    panel_manifest = []
    for row, theme in enumerate(themes):
        for col, lin in enumerate(lineages):
            ax = axes[row, col]
            letter = string.ascii_uppercase[row*len(lineages)+col]
            data = sources[(lin, theme)]
            plot_panel(ax, theme, data, lin)
            ax.text(-.14, 1.04, letter, fontsize=24, fontweight='bold', transform=ax.transAxes)
            if row == 0:
                ax.set_title('')
                ax.set_title(('Fibroblast' if lin=='fibro' else 'Endothelial') + f'\nAtlas map: {len(data):,} stratified cells', loc='left')
            name = f'{letter}_{lin}_{theme}'
            save_table(data, run/'figures/source_data'/f'{name}.csv')
            def draw_single():
                single, sax = plt.subplots(figsize=single_panel_size(theme,data),layout='constrained')
                try:
                    plot_panel(sax,theme,data,lin)
                    for ext in ['png','pdf']:
                        save_figure(single,run/'figures/panels'/f'{name}.{ext}',dpi=cfg['figure_dpi'],bbox_inches='tight')
                finally:
                    plt.close(single)
            cached_action(['panel',STYLE_VERSION,name],[run/'figures/panels'/f'{name}.{ext}' for ext in ['png','pdf']],draw_single)
            panel_manifest.append({'panel': letter, 'lineage': lin, 'theme': theme, 'n_source_rows': len(data), 'source': str(run/'figures/source_data'/f'{name}.csv')})
    # Paired forest axes use common limits for the same statistical definition.
    for row in [1, 2, 4]:
        active = [ax for ax in axes[row] if ax.axison]
        if len(active) > 1 and themes[row] != 'feature_integration':
            limits = [ax.get_xlim() for ax in active]
            for ax in active:
                ax.set_xlim(min(v[0] for v in limits), max(v[1] for v in limits))
    for ext in ['png', 'pdf']:
        target = run/'figures'/f'Fig_CellAge_Stromal_Integration.{ext}'
        cached_action(['integrated_figure',STYLE_VERSION,ext,lineages],[target],lambda:save_figure(fig,target,dpi=cfg['figure_dpi']))
    plt.close(fig)
    save_table(pd.DataFrame(panel_manifest), run/'figures/source_data/panel_manifest.csv')
    caption = '''# TabulaStromal CellAge Exploration: senescence-associated states and chronological-age clocks

Paired panels show Fibroblast (left) and Endothelial (right). A–B: existing Atlas UMAP, GSE × donor × subtype stratified visualization only (point cap recorded in source); color denotes broad annotation. C–D: mean donor hUSI by original Celltype, equal donor weighting, 95% donor-bootstrap CI. Native hUSI values are dataset-minmax scaled; shared axis limits do not imply cross-lineage calibration. E–F: within-subtype chronological-age associations adjusted for available tissue, sex, GSE, disease, treatment, platform and QC. Donor-cluster 95% CI; warm points BH q<0.05. Sample sizes and full model status are retained in source tables. G–H: six SenCID identity recommendation dimensions, each standardized across subtype donor means. Recommendation values are not senescence probabilities or intensity. I–J: '''
    caption += ('age-adjusted associations with GeneModel delta age using completed held-out K8000 ElasticNet folds only. ' if fifth == 'clock' else 'exploratory continuous comparison of hUSI age slopes and signed earlier-model stable coefficients, substituted because clock association was unavailable or non-identifiable. ')
    caption += '''K–L: local versioned pathway enrichment of earlier EN stable features using subtype-specific tested candidate-gene backgrounds, with lineage-wide BH correction; display takes at most 14 representative terms and full tested terms are saved. These are exploratory pathway hypotheses, not SenCID validation or causal mechanisms. SenCID prediction concordance and count-based DE are unavailable for current inputs. Main score statistics use verified biological donors; only visualization uses sampled cells. All associations are cross-sectional. Figure title and this caption are separate from panel content.
'''
    (run/'figures/Fig_CellAge_Stromal_Integration_caption.md').write_text(caption)
    make_reports(cfg, run, lineages, source_run=read_run)
    write_json(run/'figures/figure_style.json', {'style_version': STYLE_VERSION,
        'figure_inches': list(FIGURE_SIZE if len(lineages)==2 else (14,50)),
        'title_pt': 16, 'axis_label_pt': 14, 'subtype_label_pt': LABEL_SIZE,
        'pathway_label_pt': PATHWAY_SIZE, 'legend_pt': LABEL_SIZE, 'panel_letter_pt': 24,
        'row_height_ratios': ROW_HEIGHTS})
    stage_status(run, 'figure', 'both', 'complete', fifth_theme=fifth, n_panels=len(panel_manifest), style_version=STYLE_VERSION)
