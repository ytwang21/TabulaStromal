"""Phase 1: reusable minimal tables, donor summaries and coverage QC."""
from pathlib import Path
import numpy as np
import pandas as pd
from .common import aggregate_donors, save_table, score_columns, stage_status, write_json
from .checkpointing import cached_action, atomic_partitioned_frame


def run(cfg, run, lineage):
    run = Path(run)
    cells = pd.read_parquet(run/'cache'/f'{lineage}_audited_cells.parquet')
    target = run/'cache'/f'{lineage}_cells.parquet'
    cached_action('minimal_cell_cache',[target],lambda:save_table(cells,target))
    # Partition only the score-covered minimal table; keep the full Atlas cache for denominators.
    cols = score_columns(cells) + [c for c in cells if c.startswith('SenCID_rec_')]
    covered = cells[cols].notna().any(axis=1) if cols else pd.Series(False, index=cells.index)
    if covered.any():
        target = run/'tables'/f'{lineage}_cell_scores'
        cached_action('partitioned_scores',[target],lambda:atomic_partitioned_frame(cells.loc[covered],target,['lineage','gse']))
    write_json(run/'tables'/f'{lineage}_schema.json', {c: str(t) for c, t in cells.dtypes.items()})
    target = run/'tables'/f'{lineage}_key_mapping.csv'
    cached_action('key_mapping',[target],lambda:save_table(cells.loc[covered, ['lineage', 'cell_key', 'gse', 'gsm', 'sample', 'donor_key', 'age_source']],target))
    rows = []
    for c in cols:
        v = cells[c].dropna()
        rows.append({'lineage': lineage, 'score': c, 'n_available': len(v), 'atlas_coverage': len(v)/len(cells),
                     'minimum': v.min(), 'median': v.median(), 'maximum': v.max(), 'missing_fraction': cells[c].isna().mean()})
    save_table(pd.DataFrame(rows), run/'tables'/f'{lineage}_score_qc.csv')
    common = cells.hUSI.notna() & cells.filter(regex='^SenCID_rec_').notna().any(axis=1) if 'hUSI' in cells else pd.Series(False, index=cells.index)
    cohort = {'lineage': lineage, 'n_cells': len(cells), 'age_atlas_cells': int(cells.age_atlas_member.sum()),
              'known_Age_cells': int(cells.age.notna().sum()), 'donors': cells.donor_key.nunique(),
              'tissues': cells.tissue.nunique(), 'celltypes': cells.celltype.nunique(),
              'hUSI_cells': int(cells.hUSI.notna().sum()) if 'hUSI' in cells else 0,
              'recSID_cells': int(cells.filter(regex='^SenCID_rec_').notna().any(axis=1).sum()),
              'hUSI_recSID_common_cells': int(common.sum())}
    save_table(pd.DataFrame([cohort]), run/'tables'/f'{lineage}_cohort.csv')
    scores = score_columns(cells)
    for level in ['celltype', 'celltype_sub']:
        agg = aggregate_donors(cells, scores, cfg['min_cells'], level)
        save_table(agg, run/'tables'/f'{lineage}_donor_{level}_scores.parquet')
    stage_status(run, 'merge', lineage, 'complete', **{k:v for k,v in cohort.items() if k != 'lineage'})
