"""Donor-level subtype score ranking, paired contrasts and confounder QC."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from .common import bootstrap_mean, fdr, load_cells, save_table, stage_status


def run(cfg, run, lineage):
    run = Path(run)
    p = run/'tables'/f'{lineage}_donor_celltype_scores.parquet'
    if not p.exists():
        stage_status(run, 'subtypes', lineage, 'blocked', 'No verified donor score summaries')
        return
    agg = pd.read_parquet(p)
    rng = np.random.default_rng(cfg['seed'])
    rows, contrasts = [], []
    for (score, subtype), x in agg.groupby(['score', 'celltype'], observed=True):
        # Equal weighting of donors; average tissues within a donor for the ranking.
        d = x.groupby('donor_key', observed=True).score_mean.mean()
        lo, hi = bootstrap_mean(d, rng, cfg['bootstrap_n']) if len(d) > 1 else (np.nan, np.nan)
        rows.append({'lineage': lineage, 'score': score, 'celltype': subtype, 'mean': d.mean(),
            'ci_low': lo, 'ci_high': hi, 'n_donors': len(d), 'n_gse': x.gse.nunique(), 'n_cells': x.n_cells.sum(),
            'n_tissues': x.tissue.nunique(), 'inference': 'donor bootstrap; equal donor weighting'})
    # Each donor contributes one paired difference against other subtypes in same tissue.
    for score, x in agg.groupby('score', observed=True):
        for subtype in x.celltype.unique():
            a = x.loc[x.celltype.eq(subtype), ['donor_key', 'tissue', 'score_mean']]
            b = x.loc[~x.celltype.eq(subtype)].groupby(['donor_key', 'tissue'], observed=True, dropna=False).score_mean.mean().reset_index(name='reference_mean')
            z = a.merge(b, on=['donor_key', 'tissue'], validate='one_to_one')
            d = (z.score_mean-z.reference_mean).groupby(z.donor_key).mean()
            if len(d) < cfg['min_donors']:
                continue
            lo, hi = bootstrap_mean(d, rng, cfg['bootstrap_n'])
            pval = stats.ttest_1samp(d, 0).pvalue if d.std() > 0 else np.nan
            contrasts.append({'lineage': lineage, 'score': score, 'celltype': subtype, 'effect': d.mean(),
                'ci_low': lo, 'ci_high': hi, 'p': pval, 'n_donors': len(d),
                'comparison': 'paired within donor/tissue vs equal-weight mean of available other subtypes; exploratory identity contrast'})
    save_table(pd.DataFrame(rows), run/'tables'/f'{lineage}_subtype_ranking.csv')
    save_table(fdr(pd.DataFrame(contrasts), ['score']), run/'tables'/f'{lineage}_subtype_contrasts.csv')
    q = []
    for cov in ['gse', 'tissue', 'sex', 'disease_status', 'platform', 'treatment']:
        for (score, group), x in agg.groupby(['score', cov], observed=True, dropna=False):
            d = x.groupby('donor_key', observed=True).score_mean.mean()
            q.append({'lineage': lineage, 'score': score, 'covariate': cov, 'group': group, 'n_donors': len(d), 'mean': d.mean(), 'median': d.median()})
    save_table(pd.DataFrame(q), run/'tables'/f'{lineage}_confounder_support.csv')
    stage_status(run, 'subtypes', lineage, 'complete', n_comparisons=len(contrasts))
