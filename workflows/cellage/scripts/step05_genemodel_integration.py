"""Exact biological-donor bridge to existing held-out chronological-age predictions."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
from scipy import stats
from .common import (adjusted_fit, fdr, load_cells, read_table, save_table,
                     stage_status, write_json)


def run(cfg, run, lineage):
    run = Path(run)
    g = cfg['lineages'][lineage]['genemodel']
    path = Path(g['predictions'])
    if not path.exists():
        stage_status(run, 'clock', lineage, 'blocked', 'No existing prediction file')
        return
    # Snapshot a possibly growing benchmark table, so one analysis uses one version.
    payload = path.read_bytes()
    snapshot = run/'audit'/f'{lineage}_prediction_snapshot.csv'
    snapshot.write_bytes(payload)
    pred = read_table(snapshot)
    if g.get('model_column') in pred:
        pred = pred.loc[pred[g['model_column']].eq(g['model_value'])].copy()
    if g.get('feature_count_column') in pred:
        pred = pred.loc[pred[g['feature_count_column']].eq(g['feature_count'])].copy()
    split_col, sample_col, fold_col = g['split_column'], g['sample_column'], g['fold_column']
    held = pred.loc[pred[split_col].isin(g['allowed_splits'])].copy()
    if held.empty:
        stage_status(run, 'clock', lineage, 'blocked', 'No eligible held-out predictions for prespecified model/K')
        return
    cells = load_cells(run, lineage)
    bridge = cells.loc[cells.donor_key.notna(), ['sample', 'donor_key', 'gse']].drop_duplicates()
    # A sample reused for different biological donors cannot be resolved by sample alone.
    if bridge.groupby('sample', observed=True).donor_key.nunique().gt(1).any():
        stage_status(run, 'clock', lineage, 'blocked', 'Prediction keys ambiguous; require GSE/sample donor crosswalk')
        return
    bridge = bridge.groupby('sample', observed=True).agg(donor_key=('donor_key', 'first'), gse=('gse', 'first')).reset_index()
    pred = pred.merge(bridge, left_on=sample_col, right_on='sample', how='left', validate='many_to_one', suffixes=('', '_atlas'))
    held = pred.loc[pred[split_col].isin(g['allowed_splits'])].copy()
    # Audit disjointness per fold using all archived Train/Valid records.
    disjoint, fold_rows = True, []
    for fold, x in pred.groupby(fold_col):
        train = set(x.loc[x[split_col].astype(str).str.lower().eq('train'), 'donor_key'].dropna())
        valid = set(x.loc[x[split_col].isin(g['allowed_splits']), 'donor_key'].dropna())
        overlap = train & valid
        verified = bool(train and valid and not overlap)
        disjoint &= verified
        fold_rows.append({'lineage': lineage, 'fold': fold, 'n_train_donors': len(train), 'n_heldout_donors': len(valid),
                          'n_overlap': len(overlap), 'donor_disjoint_verified': verified})
    save_table(pd.DataFrame(fold_rows), run/'audit'/f'{lineage}_clock_fold_disjointness.csv')
    # One held-out prediction per donor in the selected CV snapshot, never average folds silently.
    if held.loc[held.donor_key.notna()].duplicated('donor_key').any():
        stage_status(run, 'clock', lineage, 'blocked', 'Repeated held-out predictions per donor; choose explicit CV repeat policy')
        return
    held['predicted_age'] = pd.to_numeric(held[g['prediction_column']], errors='coerce')
    held['clock_age'] = pd.to_numeric(held[g['age_column']], errors='coerce')
    age_bridge = cells.loc[cells.age_atlas_member & cells.donor_key.notna()].groupby('donor_key', observed=True).age.median()
    held['atlas_age'] = held.donor_key.map(age_bridge)
    held['age_consistent'] = held.atlas_age.notna() & np.isclose(held.atlas_age, held.clock_age, atol=.01, rtol=0)
    held['eligible'] = held.donor_key.notna() & held.age_consistent & held.predicted_age.notna()
    held['delta_age'] = held.predicted_age - held.clock_age
    held['source'] = str(path)
    save_table(held.assign(lineage=lineage), run/'audit'/f'{lineage}_genemodel_join_qc.csv')
    write_json(run/'audit'/f'{lineage}_clock_provenance.json', {
        'prediction_sha256': hashlib.sha256(payload).hexdigest(), 'source': str(path),
        'model': g['model_value'], 'K': g.get('feature_count'), 'heldout_n': len(held),
        'matched_n': int(held.donor_key.notna().sum()), 'eligible_n': int(held.eligible.sum()),
        'folds_present': held[fold_col].unique().tolist(), 'complete_five_fold_OOF': held[fold_col].nunique() == 5,
        'donor_disjoint_verified_from_Train_rows': disjoint,
        'train_only_preprocessing_verified': g.get('preprocessing_train_only_verified', False),
        'preprocessing_evidence': g.get('preprocessing_evidence'),
        'age_bias_adjustment': 'Not fitted on held-out test cohort; raw delta_age is adjusted for chronological age in outcome model',
        'stable_feature_source_warning': 'Stable features come from earlier downstream model; interpretation branch is separate from this benchmark model.'})
    if not disjoint or not g.get('preprocessing_train_only_verified', False):
        stage_status(run, 'clock', lineage, 'blocked', 'Donor-disjoint folds or train-only preprocessing provenance unverified; join QC retained')
        return
    a = run/'tables'/f'{lineage}_age_donor_scores.parquet'
    if not a.exists():
        stage_status(run, 'clock', lineage, 'blocked', 'No eligible donor age score summaries')
        return
    agg = pd.read_parquet(a)
    keep = held.loc[held.eligible, ['donor_key', 'delta_age', 'predicted_age', fold_col, split_col]].rename(columns={fold_col:'clock_fold',split_col:'clock_split'})
    joined = agg.merge(keep, on='donor_key', how='inner', validate='many_to_one')
    save_table(joined, run/'tables'/f'{lineage}_clock_donor_scores.parquet')
    rows = []
    for (score, subtype), x in joined.groupby(['score', 'celltype'], observed=True):
        rec = adjusted_fit(x, target='delta_age', extra=['age'], min_donors=cfg['min_donors'])
        rows.append({'lineage': lineage, 'score': score, 'celltype': subtype,
                     'prediction_context': 'partial held-out-fold snapshot, no in-sample predictions', **rec})
    save_table(fdr(pd.DataFrame(rows), ['score']), run/'tables'/f'{lineage}_clock_associations.csv')
    stage_status(run, 'clock', lineage, 'complete', n_eligible=int(held.eligible.sum()), n_fold=held[fold_col].nunique())
