"""Phase 0: inventory, metadata-only H5AD audit and exact score mapping."""
from pathlib import Path
import gc
import hashlib
import logging

import numpy as np
import pandas as pd

from .common import (ALIASES, attach_donors, clean_strings, h5_metadata,
                     read_table, save_table, stage_status, strict_map, write_json)

LOG = logging.getLogger('cellage')


def inventory(cfg, run):
    root = Path(cfg['project_root'])
    paths = set()
    for lin, conf in cfg['lineages'].items():
        paths.update(Path(conf['atlas']).parent.glob('*.h5ad'))
        for s in conf['scores'].values():
            paths.update(Path(s['path']).parent.rglob('*'))
            paths.add(Path(s['source_script']))
    original = root/'data/scalex/fibro/human/age_prediction/CellAge'
    paths.update(original.rglob('*'))
    # Bounded project search; metadata only. Do not load pickled models.
    for r in [root/'AgePrediction/GeneModel', root/'AgePrediction/GeneModelBench']:
        if r.exists():
            paths.update(r.rglob('*'))
    rows = []
    for p in sorted(paths):
        if not p.is_file() or '__pycache__' in str(p):
            continue
        s = p.stat()
        low = str(p).lower()
        lineage = 'endo' if 'endo' in low else ('fibro' if 'fibro' in low else 'unknown')
        role = 'candidate_GeneModel' if '/genemodel' in low else 'CellAge_artifact'
        if p.suffix == '.h5ad':
            role = 'age_raw_candidate' if '_age_raw' in low else ('age_atlas' if '_age.h5ad' in low else 'atlas_candidate')
        rows.append({'path': str(p), 'size': s.st_size, 'mtime': pd.Timestamp(s.st_mtime, unit='s', tz='UTC').isoformat(),
                     'file_type': p.suffix, 'lineage': lineage, 'role': role,
                     'numeric_candidate': p.suffix in {'.csv', '.tsv', '.parquet', '.h5ad', '.npz', '.pkl'},
                     'selected': any(str(p) == str(v) for c in cfg['lineages'].values() for v in [c['atlas'], c['age_atlas'], c['raw_atlas'], *[s['path'] for s in c['scores'].values()], *[v for v in c['genemodel'].values() if isinstance(v, str)]])})
    save_table(pd.DataFrame(rows), Path(run)/'audit/input_manifest.csv')
    previews = []
    for p in sorted(original.rglob('*')):
        if p.is_file() and p.suffix in {'.csv', '.tsv', '.parquet'}:
            try:
                d = read_table(p, nrows=200)
                previews.append({'path': str(p), 'scope': 'schema preview, at most 200 rows; never concatenated',
                                 'columns': '|'.join(d.columns.astype(str)), 'n_columns': len(d.columns),
                                 'preview_rows': len(d), 'possible_keys': '|'.join(c for c in d if c in {'cellid','cell_barcode','sample','GSE','GSM','Unnamed: 0'}),
                                 'score_columns': '|'.join(c for c in d if any(k in c for k in ['SID','hUSI','husi']))})
            except Exception as e:
                previews.append({'path': str(p), 'scope': 'preview_failed', 'error': repr(e)})
    save_table(pd.DataFrame(previews), Path(run)/'audit/all_CellAge_file_schemas.csv')
    return rows


def run(cfg, run, lineage):
    run = Path(run)
    conf = cfg['lineages'][lineage]
    audit = run/'audit'
    cache = run/'cache'
    cache.mkdir(exist_ok=True, parents=True)
    metadata, diagnostics, age_checks = {}, [], []
    extra = [cfg.get('donor_column')] if cfg.get('donor_column') else []
    full, genes, info = h5_metadata(conf['atlas'], extra)
    if not full.index.is_unique:
        raise ValueError('Full Atlas obs_names must be unique; supply explicit composite identity adapter')
    metadata['atlas'] = info
    save_table(pd.DataFrame({'gene': genes}), audit/f'{lineage}_atlas_genes.csv')
    if Path(conf['age_atlas']).exists():
        age_obs, _, age_info = h5_metadata(conf['age_atlas'], extra)
        metadata['age_atlas'] = age_info
        # Validate containment and sample identity, never assume subset equivalence.
        mapped_age, qc = strict_map(full, age_obs.reset_index(names='cellid'), min_rate=cfg['min_match_rate'])
        qc.update(lineage=lineage, method='age_subset')
        diagnostics.append(qc)
        if mapped_age is None:
            raise ValueError(f'Age subset conflicts with full Atlas: {qc}')
    else:
        age_obs = pd.DataFrame(index=pd.Index([], dtype=str))
    # age_raw is audited only for metadata and matrix sample, no DE by filename.
    if Path(conf['raw_atlas']).exists():
        raw_obs, _, raw_info = h5_metadata(conf['raw_atlas'], extra)
        raw_info['obs_exactly_equal_age_subset'] = bool(raw_obs.index.equals(age_obs.index))
        metadata['raw_atlas'] = raw_info
        del raw_obs
    write_json(audit/f'{lineage}_h5ad_schema.json', metadata)

    # Age is canonical, lowercase age is inspected, never silently substituted.
    for role, x in [('atlas', full), ('age_atlas', age_obs)]:
        if x.empty:
            continue
        upper = pd.to_numeric(x.get('Age', pd.Series(np.nan, index=x.index)), errors='coerce')
        lower = pd.to_numeric(x.get('age', pd.Series(np.nan, index=x.index)), errors='coerce')
        both = upper.notna() & lower.notna()
        age_checks.append({'lineage': lineage, 'role': role, 'n_cells': len(x), 'Age_numeric': int(upper.notna().sum()),
                           'age_numeric': int(lower.notna().sum()), 'both_numeric': int(both.sum()),
                           'different_values': int((both & ~np.isclose(upper, lower, equal_nan=False)).sum()),
                           'Age_missing_age_present': int((upper.isna() & lower.notna()).sum()),
                           'Age_min': upper.min(), 'Age_max': upper.max(), 'units': cfg['age_units']})
    canonical = full.rename(columns=ALIASES).copy()
    for c in ALIASES.values():
        if c not in canonical:
            canonical[c] = pd.NA
    canonical['age_atlas_member'] = canonical.index.isin(age_obs.index)
    canonical['atlas_Age'] = pd.to_numeric(full.get('Age'), errors='coerce')
    canonical['age'] = canonical['atlas_Age']
    canonical['age_source'] = 'full_atlas.Age'
    canonical['age_conflict'] = False
    if 'Age' in age_obs:
        subset_age = pd.to_numeric(age_obs.Age, errors='coerce').reindex(canonical.index)
        conflict = canonical.atlas_Age.notna() & subset_age.notna() & ~np.isclose(canonical.atlas_Age, subset_age, equal_nan=False)
        canonical['age_conflict'] = conflict
        use = canonical.age_atlas_member & subset_age.notna()
        canonical.loc[use, 'age'] = subset_age.loc[use]
        canonical.loc[use, 'age_source'] = 'age_atlas.Age (explicit subset annotation)'
        canonical.loc[conflict, 'age'] = np.nan
        age_checks.append({'lineage': lineage, 'role': 'cross_file_Age', 'different_values': int(conflict.sum()),
                           'subset_fills_full_missing': int((use & canonical.atlas_Age.isna()).sum()),
                           'policy': 'Use documented age-subset Age; exclude numeric conflicts; never fallback to lowercase age'})
    if cfg['age_units'] != 'years':
        raise ValueError('Convert ages explicitly to years in an audited adapter before analysis')
    canonical.loc[(canonical.age < 0) | (canonical.age > 120), 'age'] = np.nan
    for c in ['gse', 'gsm', 'sample', 'tissue', 'sex', 'disease_status', 'treatment', 'platform']:
        canonical[c] = clean_strings(canonical[c])
    canonical['lineage'] = lineage
    canonical['cell_key'] = canonical.index
    canonical['atlas_row'] = np.arange(len(canonical), dtype=np.int32)
    canonical = attach_donors(canonical, cfg, lineage)
    # donor consistency is checked across the entire Atlas, not just score-positive cells.
    good = canonical.loc[canonical.donor_key.notna()]
    donor_qc = good.groupby('donor_key', observed=True).agg(n_cells=('cell_key', 'size'), n_samples=('sample', 'nunique'),
                    n_gse=('gse', 'nunique'), n_tissues=('tissue', 'nunique'), n_ages=('age', 'nunique'), age=('age', 'median'),
                    n_sexes=('sex', 'nunique'))
    donor_qc['age_conflict'] = donor_qc.n_ages > 1
    donor_qc['identity_conflict'] = donor_qc.n_sexes > 1
    bad = donor_qc.index[donor_qc.age_conflict | donor_qc.identity_conflict]
    canonical['donor_qc_pass'] = canonical.donor_key.notna() & ~canonical.donor_key.isin(bad)
    canonical.loc[~canonical.donor_qc_pass, 'donor_key'] = pd.NA
    save_table(donor_qc.reset_index().assign(lineage=lineage), audit/f'{lineage}_donor_qc.csv')
    sq = canonical.groupby(['gse', 'sample'], observed=True, dropna=False).agg(n_cells=('cell_key', 'size'), n_age_cells=('age', 'count'),
        n_ages=('age', 'nunique'), age=('age', 'median'), n_tissues=('tissue', 'nunique'), n_donors=('donor_key', 'nunique'),
        missing_disease=('disease_status', lambda s: s.isna().mean()), missing_sex=('sex', lambda s: s.isna().mean()),
        n_celltypes=('celltype', 'nunique'))
    save_table(sq.reset_index().assign(lineage=lineage), audit/f'{lineage}_sample_qc.csv')
    save_table(canonical.groupby(['gse', 'sample', 'celltype'], observed=True, dropna=False).size().rename('n_cells').reset_index(), audit/f'{lineage}_sample_subtype_counts.csv')
    schemas, usable = [], []
    for method, spec in conf['scores'].items():
        script = Path(spec['source_script'])
        if script.is_file():
            payload = script.read_bytes()
            (audit/f'{lineage}_{method}_source_script.txt').write_bytes(payload)
            write_json(audit/f'{lineage}_{method}_semantics.json', {'source_script': str(script),
                'sha256': hashlib.sha256(payload).hexdigest(), 'semantics': spec['semantics'],
                'official_reference': 'https://github.com/JackieHanLab/SenCID/blob/main/SenCID/api.py' if method == 'SenCID' else 'https://github.com/WJPina/HUSI',
                'checked_on': '2026-10-30', 'source_values': spec['values']})
        if spec['semantics'] != 'senescence_score' and 'SenCID_score' in spec['values'].values():
            raise ValueError('Recommendation outputs cannot be relabeled as a SenCID senescence score')
        if not Path(spec['path']).exists():
            diagnostics.append({'lineage': lineage, 'method': method, 'passed': False, 'reason': 'missing_file'})
            continue
        scores = read_table(spec['path'])
        mapped, qc = strict_map(full, scores, spec['cell_column'], cfg['min_match_rate'])
        qc.update(lineage=lineage, method=method, source=spec['path'])
        diagnostics.append(qc)
        for c in scores:
            numeric = pd.to_numeric(scores[c], errors='coerce')
            schemas.append({'lineage': lineage, 'method': method, 'path': spec['path'], 'n_rows': len(scores), 'n_columns': len(scores.columns),
                            'column': c, 'dtype': str(scores[c].dtype), 'missing_fraction': scores[c].isna().mean(),
                            'numeric_min': numeric.min(), 'numeric_max': numeric.max(), 'n_unique': scores[c].nunique(),
                            'duplicate_cell_keys': int(scores.duplicated(spec['cell_column']).sum()), 'semantics': spec['semantics']})
        if mapped is not None:
            # Direct cell IDs were verified globally unique; map with explicit keys.
            if not mapped.cell_key.is_unique:
                raise ValueError('Composite mapping needs an explicit canonical cell ID adapter')
            d = mapped.set_index('cell_key')
            aligned = d.reindex(canonical.cell_key)
            for src, dst in spec['values'].items():
                if src in d:
                    canonical[dst] = pd.to_numeric(aligned[src], errors='coerce').to_numpy()
            usable.append(method)
        del scores, mapped
        gc.collect()
    canonical['quality_flags'] = np.where(canonical.donor_qc_pass, 'donor_verified', 'donor_missing_or_conflicting')
    canonical['provenance'] = conf['atlas']
    canonical = canonical.drop(columns=[c for c in ['Age', 'donor', 'donor_id'] if c in canonical])
    # Categorical encoding reduces RAM and parquet size without changing identifiers.
    for c in canonical.select_dtypes(include=['object', 'string']).columns:
        if c not in ['cell_key', 'donor_key']:
            canonical[c] = canonical[c].astype('category')
    save_table(canonical.reset_index(drop=True), cache/f'{lineage}_audited_cells.parquet')
    save_table(pd.DataFrame(diagnostics), audit/f'{lineage}_mapping_qc.csv')
    save_table(pd.DataFrame(schemas), audit/f'{lineage}_score_schema.csv')
    save_table(pd.DataFrame(age_checks), audit/f'{lineage}_age_checks.csv')
    # Candidate prediction schema is inspected even if donor integration is blocked.
    gm = conf['genemodel']
    prediction_schema = {}
    if Path(gm['predictions']).exists():
        pred = read_table(gm['predictions'])
        prediction_schema = {'path': gm['predictions'], 'shape': pred.shape, 'columns': list(pred),
                             'splits': pred[gm['split_column']].value_counts().to_dict(), 'split_evidence': gm['split_evidence'],
                             'donor_disjoint_verified': gm['donor_disjoint_verified']}
    write_json(audit/f'{lineage}_genemodel_schema.json', prediction_schema)
    available = {'lineage': lineage, 'score_mapping': 'feasible' if usable else 'blocked',
        'donor_identity': 'feasible' if canonical.donor_qc_pass.any() else 'blocked',
        'age_associations': 'feasible' if canonical.loc[canonical.donor_qc_pass, 'age'].notna().any() else 'blocked',
        'SenCID_senescence_validation': 'blocked: only recSID recommendation scores saved; pred_dict SID_Score absent',
        'SID_landscape': 'feasible' if 'SenCID' in usable else 'blocked',
        'clock_association': 'conditional: verify donor-disjoint prediction provenance, exact IDs and ages',
        'expression_DE': 'conditional: verified counts provenance and donor support required; fractional raw values found in local pilot',
        'pathways': 'conditional: tested candidate feature universe and local versioned GMT required'}
    write_json(audit/f'{lineage}_available_analyses.json', available)
    (audit/f'{lineage}_audit_summary.md').write_text(
        f'# {lineage} Phase 0\n\nAtlas: {len(canonical):,} cells; {len(genes):,} genes.\n'
        f'Age subset: {len(age_obs):,} cells; donor QC pass: {canonical.loc[canonical.donor_qc_pass, "donor_key"].nunique():,} donors.\n'
        f'Donor evidence: {cfg.get("donor_identity_evidence", "Not verified")}\n\n'
        'Cell mapping uses exact cell IDs and annotation consistency; no positional matching or suffix stripping.\n'
        'Age subset membership does not guarantee a known numerical Age; see age_checks.csv.\n'
        'SenCID columns are recommendation values, NOT SID probabilities or an overall senescence score.\n'
        'Existing hUSI is dataset-minmax scaled; absolute values are not calibrated across lineages. SSE is dataset-adaptive; GMM is not binary.\n'
        'The _age_raw matrix has sampled fractional values; its filename alone does not certify count suitability.\n')
    stage_status(run, 'audit', lineage, 'complete', n_cells=len(canonical), n_donors=int(canonical.donor_key.nunique()))
    del canonical, full, age_obs
    gc.collect()


def finalize(cfg, run, lineages):
    run = Path(run)
    for name in ['mapping_qc', 'sample_qc', 'age_checks', 'score_schema']:
        frames = [pd.read_csv(run/'audit'/f'{lin}_{name}.csv') for lin in lineages if (run/'audit'/f'{lin}_{name}.csv').exists()]
        if frames:
            save_table(pd.concat(frames, ignore_index=True), run/'audit'/f'{name}.csv')
    from .common import read_json
    available = {lin: read_json(run/'audit'/f'{lin}_available_analyses.json') for lin in lineages}
    write_json(run/'audit/available_analyses.json', available)
    (run/'audit/data_dictionary.md').write_text('''# Data dictionary

`lineage`: fibro=Fibroblast, endo=Endothelial; physical H5AD source defines lineage.
`cell_key`: exact original obs_names; no stripping of .1/-1 suffixes.
`donor_key`: sample per user confirmation, or explicit documented crosswalk.
`age`: chronological age in years; Age preferred, age subset Age has explicit provenance.
`atlas_Age`: unmodified numeric Age from full Atlas. `age_conflict`: contradictory numeric Age across files, excluded.
`age_atlas_member`: exact verified membership, separate from availability of numeric age.
`hUSI`: existing minmax-scaled rank-correlation score, higher is more senescence-associated; not years.
`husi_SSE`: native dataset-adaptive binary classification under two-state assumption.
`husi_GMM`: native ordered multi-state labels; not a binary senescent call.
`SenCID_rec_SID1..6`: recSID recommendations for selecting identities, not senescence intensity/probabilities.
`celltype`/`celltype_sub`/`celltype_C`: original annotation levels; not interchangeable.
`delta_age`: predicted chronological age minus observed age; only validated held-out predictions are eligible for main inference.
`n_cells`: technical support; `n_donors`: independent biological replication.
Missing values remain missing; no imputed zero scores or fabricated output columns.
''')
    (run/'audit/audit_summary.md').write_text('\n\n'.join((run/'audit'/f'{lin}_audit_summary.md').read_text() for lin in lineages))
