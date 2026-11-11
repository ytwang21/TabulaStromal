"""Chronological-age slopes, composition and prespecified robustness checks."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from .common import (adjusted_fit, aggregate_donors, fdr, load_cells,
                     save_table, score_columns, stage_status)


def fit_groups(agg, cfg, sensitivity, level='celltype', adjusted=True):
    rows = []
    if agg.empty:
        return rows
    for (score, subtype), x in agg.groupby(['score', level], observed=True):
        x = x.loc[x.age.notna()].copy()
        rec = adjusted_fit(x, min_donors=cfg['min_donors'], adjusted=adjusted)
        d = x.groupby('donor_key', observed=True).agg(age=('age', 'first'), value=('score_mean', 'mean'))
        rho = stats.spearmanr(d.age, d.value).statistic if len(d) >= cfg['min_donors'] and d.age.nunique() > 3 and d.value.nunique() > 1 else np.nan
        rows.append({'score': score, level: subtype, 'sensitivity': sensitivity, 'adjusted': adjusted,
                     'spearman_rho': rho, 'age_min': x.age.min(), 'age_max': x.age.max(),
                     'scope': 'multi-GSE association' if x.gse.nunique() >= cfg['min_gse'] else 'cohort-specific exploratory', **rec})
    return rows


def run(cfg, run, lineage):
    run = Path(run)
    cells = load_cells(run, lineage)
    # Only verified age-subset members with numerical Age and consistent donor metadata.
    valid = cells.loc[cells.age_atlas_member & cells.age.notna() & cells.donor_qc_pass].copy()
    if valid.donor_key.nunique() < cfg['min_donors']:
        stage_status(run, 'age', lineage, 'blocked', 'Too few verified donors with numerical Age')
        return
    rows = []
    primary = None
    for threshold in cfg['cell_threshold_sensitivity']:
        agg = aggregate_donors(valid, score_columns(valid), threshold)
        rows.extend(fit_groups(agg, cfg, f'min_cells_{threshold}', adjusted=False))
        rows.extend(fit_groups(agg, cfg, f'min_cells_{threshold}', adjusted=True))
        if threshold == cfg['min_cells']:
            primary = agg
    if primary is None or primary.empty:
        stage_status(run, 'age', lineage, 'blocked', 'No donor-subtype groups meet cell minimum')
        return
    healthy = primary.loc[primary.disease_status.astype(str).str.lower().isin(['healthy', 'normal'])]
    rows.extend(fit_groups(healthy, cfg, 'healthy_only', adjusted=True))
    # Leave-one-GSE-out is a directional descriptive robustness test, not 1000s of extra p-values.
    logo = []
    for (score, subtype), x in primary.groupby(['score', 'celltype'], observed=True):
        if x.gse.nunique() < cfg['min_gse']:
            continue
        for gse in x.gse.dropna().unique():
            z = x.loc[x.gse.ne(gse)].groupby('donor_key', observed=True).agg(age=('age', 'first'), value=('score_mean', 'mean')).dropna()
            if len(z) >= cfg['min_donors'] and z.age.nunique() >= 4:
                slope = stats.linregress(z.age, z.value).slope
                logo.append({'lineage': lineage, 'score': score, 'celltype': subtype, 'excluded_gse': gse,
                             'unadjusted_slope': slope, 'n_donors': len(z), 'definition': 'descriptive direction sensitivity; no adjusted-replication claim'})
    for tissue, x in primary.groupby('tissue', observed=True):
        rows.extend(fit_groups(x, cfg, 'tissue:' + str(tissue), adjusted=True))
    sub = aggregate_donors(valid, score_columns(valid), cfg['min_cells'], 'celltype_sub')
    sub_rows = fit_groups(sub, cfg, 'primary', 'celltype_sub', adjusted=True)
    save_table(fdr(pd.DataFrame(sub_rows), ['score', 'sensitivity', 'adjusted']).assign(lineage=lineage), run/'tables'/f'{lineage}_age_celltype_sub.csv')
    effects = fdr(pd.DataFrame(rows), ['score', 'sensitivity', 'adjusted']).assign(lineage=lineage)
    save_table(effects, run/'tables'/f'{lineage}_age_associations.csv')
    save_table(pd.DataFrame(logo), run/'tables'/f'{lineage}_age_leave_one_GSE_out.csv')
    save_table(primary, run/'tables'/f'{lineage}_age_donor_scores.parquet')
    # Composition denominator is ALL age-eligible Atlas cells, regardless of score coverage.
    counts = valid.groupby(['donor_key', 'tissue', 'celltype'], observed=True, dropna=False).size().rename('n_cells').reset_index()
    totals = counts.groupby(['donor_key', 'tissue'], observed=True, dropna=False).n_cells.sum().rename('total_cells').reset_index()
    donor_meta = valid.groupby(['donor_key', 'tissue'], observed=True, dropna=False).agg(age=('age', 'first'), gse=('gse', 'first'), sex=('sex', 'first'), disease_status=('disease_status', 'first')).reset_index()
    # Explicit zero counts for subtypes absent from a donor; zeros here are genuine observed absence.
    grid = donor_meta.merge(pd.DataFrame({'celltype': valid.celltype.dropna().unique()}), how='cross')
    comp = grid.merge(counts, on=['donor_key', 'tissue', 'celltype'], how='left').merge(totals, on=['donor_key', 'tissue'])
    comp['n_cells'] = comp.n_cells.fillna(0)
    comp['fraction'] = comp.n_cells / comp.total_cells
    save_table(comp.assign(lineage=lineage), run/'tables'/f'{lineage}_composition.csv')
    cr = []
    for subtype, x in comp.groupby('celltype', observed=True):
        cr.append({'lineage': lineage, 'celltype': subtype, 'metric': 'observed subtype fraction, independent of senescence score',
                   **adjusted_fit(x, outcome='fraction', min_donors=cfg['min_donors'])})
    save_table(fdr(pd.DataFrame(cr), []), run/'tables'/f'{lineage}_composition_age_effects.csv')
    # Age bins follow observed donor quantiles, retain boundaries; plots only.
    unique = primary[['donor_key', 'age']].drop_duplicates('donor_key').dropna()
    if unique.age.nunique() >= 4:
        bins, edges = pd.qcut(unique.age, q=min(4, unique.age.nunique()), duplicates='drop', retbins=True)
        unique['age_bin'] = bins.astype(str)
        save_table(unique, run/'tables'/f'{lineage}_age_bins.csv')
        save_table(pd.DataFrame({'edge_years': edges}), run/'tables'/f'{lineage}_age_bin_edges.csv')
    stage_status(run, 'age', lineage, 'complete', n_models=len(effects), n_identifiable=int(effects.status.eq('ok').sum()))
