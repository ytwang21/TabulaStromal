"""SID recommendation landscape; true score concordance only when available."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from .common import bootstrap_mean, load_cells, save_table, stage_status


def run(cfg, run, lineage):
    run = Path(run)
    cells = load_cells(run, lineage)
    recs = [c for c in cells if c.startswith('SenCID_rec_')]
    rows, correlations = [], []
    for c in recs:
        sub = cells.loc[cells[c].notna() & cells.donor_key.notna()]
        d = sub.groupby(['donor_key', 'celltype'], observed=True)[c].mean().reset_index()
        for subtype, x in d.groupby('celltype', observed=True):
            rows.append({'lineage': lineage, 'celltype': subtype, 'recommendation': c, 'mean_recommendation': x[c].mean(),
                         'n_donors': x.donor_key.nunique(), 'definition': 'identity recommendation; NOT senescence intensity'})
        if 'hUSI' in cells:
            common = sub.loc[sub.hUSI.notna()]
            for donor, x in common.groupby('donor_key', observed=True):
                if len(x) >= cfg['min_cells'] and x[c].nunique() > 1 and x.hUSI.nunique() > 1:
                    correlations.append({'lineage': lineage, 'donor_key': donor, 'recommendation': c, 'n_common_cells': len(x),
                         'rho': stats.spearmanr(x.hUSI, x[c]).statistic,
                         'definition': 'within-donor descriptive hUSI vs SID recommendation association; not senescence-method concordance'})
    save_table(pd.DataFrame(rows), run/'tables'/f'{lineage}_SID_recommendation_landscape.csv')
    save_table(pd.DataFrame(correlations), run/'tables'/f'{lineage}_hUSI_recSID_associations.csv')
    if 'SenCID_score' not in cells:
        stage_status(run, 'concordance', lineage, 'conditional', 'SID landscape completed; true SenCID senescence scores and binary labels absent; no high-high/discordance calls')
        return
    # This branch is usable with an explicitly configured true SID_Score export.
    both = cells.loc[cells.hUSI.notna() & cells.SenCID_score.notna() & cells.donor_key.notna()]
    out = []
    for keys, x in both.groupby(['donor_key', 'celltype', 'tissue'], observed=True, dropna=False):
        if len(x) >= cfg['min_cells'] and x.hUSI.nunique() > 1 and x.SenCID_score.nunique() > 1:
            out.append(dict(zip(['donor_key', 'celltype', 'tissue'], keys), rho=stats.spearmanr(x.hUSI, x.SenCID_score).statistic, n_common_cells=len(x)))
    save_table(pd.DataFrame(out), run/'tables'/f'{lineage}_true_score_concordance.csv')
    stage_status(run, 'concordance', lineage, 'complete', n_donor_strata=len(out))
