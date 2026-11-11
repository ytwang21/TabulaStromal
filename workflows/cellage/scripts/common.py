"""Shared I/O, identity checks, donor aggregation and statistical models.

No expression matrix is read to obtain metadata. No barcode normalization or
positional joins are permitted. Missing inference is represented by a diagnosis.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import platform
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.linalg import qr
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from .checkpointing import atomic_json, atomic_path, cached_frame, cached_value, ACTIVE

try:
    from anndata.io import read_elem
except ImportError:
    from anndata._io.specs import read_elem

LOG = logging.getLogger("cellage")
ALIASES = {
    'Celltype': 'celltype', 'Celltype_sub': 'celltype_sub',
    'Celltype_C': 'celltype_C', 'Celltype_general': 'celltype_general',
    'GSE': 'gse', 'GSM': 'gsm', 'sample': 'sample', 'tissue': 'tissue',
    'sex': 'sex', 'Disease_Status': 'disease_status',
    'disease_global': 'disease_global', 'treat': 'treatment', 'platform': 'platform',
    'total_counts': 'total_counts', 'pct_counts_mt': 'pct_counts_mt',
    'n_genes_by_counts': 'n_genes_by_counts',
}
MISSING = {'', 'nan', 'none', 'unknown', 'na', 'n/a', 'not available'}


def write_json(path, obj):
    atomic_json(path, obj)


def read_json(path):
    return json.loads(Path(path).read_text())


def prepare_parquet(df):
    df = df.copy(deep=False)
    for c in df.select_dtypes(include='category'):
        df[c] = pd.Categorical.from_codes(df[c].cat.codes, categories=df[c].cat.categories.tolist(), ordered=bool(df[c].cat.ordered))
    return df


def save_table(df, path):
    """Do not write headerless/empty artifacts that could imply stage success."""
    if df is None or df.empty:
        return False
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with atomic_path(p) as tmp:
        if p.suffix == '.parquet':
            prepare_parquet(df).to_parquet(tmp, index=False)
        else:
            df.to_csv(tmp, index=False)
    return True


def read_table(path, nrows=None):
    p = Path(path)
    if p.suffix == '.parquet':
        d = pd.read_parquet(p)
        return d if nrows is None else d.head(nrows)
    if p.suffix not in {'.csv', '.tsv', '.txt'}:
        raise ValueError(f'Unsupported numeric table: {p}; pickle is inventoried, never auto-executed')
    return pd.read_csv(p, sep='\t' if p.suffix in {'.tsv', '.txt'} else ',', nrows=nrows)


def clean_strings(s):
    s = s.astype('string')
    return s.mask(s.str.lower().isin(MISSING))


def h5_metadata(path, extra_columns=()):
    """Read selected obs columns and small matrix samples directly with HDF5."""
    with h5py.File(path, 'r') as f:
        g = f['obs']
        names = read_elem(g[g.attrs['_index']])
        cols = set(ALIASES) | {'Age', 'age', 'donor', 'donor_id', 'donor_key'} | set(extra_columns)
        obs = pd.DataFrame({c: read_elem(g[c]) for c in cols if c in g}, index=names)
        obs.index = obs.index.astype(str)
        var = f['var']
        genes = np.asarray(read_elem(var[var.attrs['_index']])).astype(str)
        x = f['X']
        sample = np.asarray(x['data'][:100000] if isinstance(x, h5py.Group) else x[:min(4, len(obs)), :])
        info = {
            'path': str(path), 'n_obs': len(obs), 'n_vars': len(genes),
            'unique_obs_names': bool(obs.index.is_unique), 'obs_columns': list(g),
            'var_names_head': genes[:10].tolist(), 'var_names_unique': bool(pd.Index(genes).is_unique),
            'matrix_encoding': str(x.attrs.get('encoding-type', 'dense')),
            'layers': list(f.get('layers', {})), 'raw_available': 'raw' in f,
            'embeddings': list(f.get('obsm', {})), 'uns_keys': list(f.get('uns', {})),
            'sampled_matrix_min': float(np.nanmin(sample)) if sample.size else None,
            'sampled_matrix_max': float(np.nanmax(sample)) if sample.size else None,
            'sampled_noninteger_fraction': float(np.mean(np.abs(sample - np.round(sample)) > 1e-6)) if sample.size else None,
            'matrix_sample_note': 'First 100000 stored values; diagnostic only, not a counts provenance certificate',
        }
    return obs, genes, info


def strict_map(atlas, scores, cell_column='cellid', min_rate=.95):
    """Use exact globally unique cell ID; verify sample and annotations if present.

    If IDs repeat, require shared GSE/GSM/sample fields forming a unique composite.
    Match rate denominator is score rows; Atlas coverage is separately reported.
    """
    a = atlas.reset_index(names='cell_key')
    s = scores.rename(columns={cell_column: 'cell_key'}).copy()
    if s['cell_key'].isna().any():
        raise ValueError('Null score cell IDs')
    a['cell_key'] = a.cell_key.astype(str)
    s['cell_key'] = s.cell_key.astype(str)
    keys = ['cell_key']
    if not a.cell_key.is_unique or not s.cell_key.is_unique:
        keys += [c for c in ['GSE', 'GSM', 'sample'] if c in a and c in s]
    if a[keys].isna().any().any() or s[keys].isna().any().any():
        raise ValueError(f'Null composite keys: {keys}')
    adup, sdup = int(a.duplicated(keys).sum()), int(s.duplicated(keys).sum())
    if adup or sdup:
        raise ValueError(f'Mapping is not one-to-one using {keys}: Atlas duplicates={adup}, score duplicates={sdup}')
    check_cols = [c for c in ['sample', 'Celltype', 'Celltype_sub', 'GSE', 'GSM'] if c in a and c in s and c not in keys]
    joined = s.merge(a[keys + check_cols], on=keys, how='left', suffixes=('', '_atlas'), indicator=True, validate='one_to_one')
    matched = joined['_merge'].eq('both')
    conflicts = pd.Series(False, index=joined.index)
    for c in check_cols:
        valid = joined[c].notna() & joined[c + '_atlas'].notna() & matched
        conflicts |= valid & joined[c].astype(str).ne(joined[c + '_atlas'].astype(str))
    rate = float(matched.mean()) if len(s) else 0.
    qc = {'n_scores': len(s), 'n_atlas': len(a), 'n_matched': int(matched.sum()),
          'match_rate': rate, 'atlas_coverage': float(matched.sum() / len(a)),
          'duplicate_atlas_keys': adup, 'duplicate_score_keys': sdup,
          'metadata_conflicts': int(conflicts.sum()), 'keys': '|'.join(keys),
          'unmatched_examples': '|'.join(s.loc[~matched, 'cell_key'].head(5)),
          'passed': bool(rate >= min_rate and not conflicts.any())}
    if not qc['passed']:
        return None, qc
    return joined.loc[matched, list(s.columns)].copy(), qc


def attach_donors(df, cfg, lineage):
    df = df.copy()
    df['donor_key'] = pd.Series(pd.NA, index=df.index, dtype='string')
    cross = cfg.get('donor_crosswalk')
    if cross:
        c = read_table(cross)
        required = {'lineage', 'gse', 'sample', 'donor_key', 'evidence'}
        if not required.issubset(c):
            raise ValueError(f'Donor crosswalk requires {required}')
        c = c.loc[c.lineage.eq(lineage)]
        if c.duplicated(['gse', 'sample']).any() or c.evidence.isna().any():
            raise ValueError('Ambiguous or undocumented donor crosswalk')
        df = df.drop(columns='donor_key').merge(c[['gse', 'sample', 'donor_key']], on=['gse', 'sample'], how='left', validate='many_to_one')
    elif cfg.get('donor_column') and cfg.get('donor_identity_evidence'):
        col = ALIASES.get(cfg['donor_column'], cfg['donor_column'])
        if col not in df:
            raise ValueError(f'Configured donor column missing: {col}')
        # Do not prefix with GSE: one donor can occur in multiple studies/tissues.
        df['donor_key'] = clean_strings(df[col])
    return df


def score_columns(df):
    return [c for c in ['hUSI', 'SenCID_score'] if c in df and df[c].notna().any()]


def aggregate_donors(cells, scores, min_cells=20, group_col='celltype'):
    return cached_frame(['donor_aggregation', scores, min_cells, group_col],
                        lambda: _aggregate_donors_impl(cells,scores,min_cells,group_col))


def _aggregate_donors_impl(cells, scores, min_cells=20, group_col='celltype'):
    """One row per donor × tissue × subtype × method; no cell-level inference."""
    rows = []
    valid = cells.loc[cells.donor_key.notna() & cells[group_col].notna()].copy()
    if valid.empty:
        return pd.DataFrame()
    for keys, x in valid.groupby(['lineage', 'donor_key', 'tissue', group_col], observed=True, dropna=False, sort=False):
        age_values = x.age.dropna().unique()
        base = dict(zip(['lineage', 'donor_key', 'tissue', group_col], keys))
        base['age'] = float(age_values[0]) if len(age_values) == 1 else np.nan
        for c in ['gse', 'sex', 'disease_status', 'treatment', 'platform']:
            vals = x[c].dropna().unique() if c in x else []
            base[c] = str(vals[0]) if len(vals) == 1 else ('mixed' if len(vals) > 1 else pd.NA)
        base['n_gse_within_donor'] = x.gse.nunique()
        for c in ['total_counts', 'pct_counts_mt', 'n_genes_by_counts']:
            if c in x:
                nv = pd.to_numeric(x[c], errors='coerce').dropna()
                base[c] = nv.median() if len(nv) else np.nan
        for score in scores:
            v = pd.to_numeric(x[score], errors='coerce').dropna()
            if len(v) < min_cells:
                continue
            rec = {**base, 'score': score, 'n_cells': len(v), 'score_mean': v.mean(),
                   'score_median': v.median(), 'score_iqr': v.quantile(.75) - v.quantile(.25)}
            if score == 'hUSI' and 'husi_SSE' in x:
                b = pd.to_numeric(x.loc[x[score].notna(), 'husi_SSE'], errors='coerce').dropna()
                if len(b) and set(b.unique()) <= {0, 1}:
                    rec['native_SSE_high_fraction'] = b.mean()
            rows.append(rec)
    return pd.DataFrame(rows)


def fdr(df, family, p='p', q='q'):
    if df.empty or p not in df:
        return df
    df = df.copy()
    df[q] = np.nan
    groups = df.groupby(family, dropna=False).groups.values() if family else [df.index]
    for ix in groups:
        ix = df.loc[ix].index[df.loc[ix, p].notna()]
        if len(ix):
            df.loc[ix, q] = multipletests(df.loc[ix, p], method='fdr_bh')[1]
    return df


def adjusted_fit(frame, outcome='score_mean', target='age', min_donors=12, extra=(), adjusted=True):
    if ACTIVE.get() is None:
        return _adjusted_fit_impl(frame,outcome,target,min_donors,extra,adjusted)
    token = hashlib.sha256(pd.util.hash_pandas_object(frame,index=True).values.tobytes()).hexdigest()
    return cached_value(['regression',token,list(frame.columns),outcome,target,min_donors,list(extra),adjusted],
                        lambda: _adjusted_fit_impl(frame,outcome,target,min_donors,extra,adjusted))


def _adjusted_fit_impl(frame, outcome='score_mean', target='age', min_donors=12, extra=(), adjusted=True):
    """OLS with donor-cluster robust uncertainty; reject aliased designs.

    Fixed GSE/tissue covariates, complete cases, and donor-cluster covariance.
    Return explicit non-identifiability rather than dropping confounders to get p.
    """
    need = list(dict.fromkeys([outcome, target, 'donor_key'] + list(extra)))
    x = frame.copy()
    for c in [outcome, target, *extra]:
        x[c] = pd.to_numeric(x[c], errors='coerce')
    cats, nums = [], []
    if adjusted:
        cats = [c for c in ['tissue', 'sex', 'gse', 'disease_status', 'treatment', 'platform'] if c in x and x[c].nunique() > 1]
        nums = [c for c in ['total_counts', 'pct_counts_mt', 'n_genes_by_counts'] if c in x and x[c].nunique() > 1]
        for c in cats:
            x[c] = clean_strings(x[c])
        need += cats + nums
    x = x.dropna(subset=need)
    rec = {'n_rows': len(x), 'n_donors': x.donor_key.nunique(), 'n_gse': x.gse.nunique() if 'gse' in x else 0,
           'effect': np.nan, 'ci_low': np.nan, 'ci_high': np.nan, 'p': np.nan,
           'adjustment': ','.join(cats + nums + list(extra)), 'status': 'insufficient_support'}
    if rec['n_donors'] < min_donors or x[target].nunique() < 4:
        return rec
    design = x[[target, *extra, *nums]].astype(float)
    # Scale numeric nuisance variables for conditioning, keep target in native units.
    for c in nums:
        sd = design[c].std()
        if sd > 0:
            design[c] = (design[c] - design[c].mean()) / sd
    if cats:
        design = pd.concat([design, pd.get_dummies(x[cats].astype(str), drop_first=True, dtype=float)], axis=1)
    design = sm.add_constant(design, has_constant='add').astype(float)
    # Many metadata factors nest exactly (e.g. platform within GSE). Keep a
    # full-rank basis of their joint nuisance span, preserving all adjustment.
    nuisance = design.drop(columns=target)
    _, r, piv = qr(nuisance.to_numpy(), mode='economic', pivoting=True)
    diag = np.abs(np.diag(r))
    tol = max(nuisance.shape) * np.finfo(float).eps * (diag.max() if len(diag) else 1.)
    rank = int(np.sum(diag > tol))
    reduced = nuisance.iloc[:, sorted(piv[:rank])]
    target_residual = design[target].to_numpy() - reduced.to_numpy() @ np.linalg.lstsq(reduced, design[target], rcond=None)[0]
    if np.linalg.norm(target_residual) <= 1e-8 * max(np.linalg.norm(design[target]), 1):
        rec['status'] = 'not_identifiable: target_fully_confounded'
        return rec
    rec['n_redundant_nuisance_columns'] = nuisance.shape[1] - rank
    design = pd.concat([design[[target]], reduced], axis=1)
    if len(x) - design.shape[1] < 8:
        rec['status'] = 'not_identifiable: aliased_covariates_or_insufficient_residual_df'
        return rec
    fit = sm.OLS(x[outcome].astype(float), design).fit(cov_type='cluster', cov_kwds={'groups': x.donor_key}, use_t=True)
    ci = fit.conf_int().loc[target]
    rec.update(effect=float(fit.params[target]), ci_low=float(ci.iloc[0]), ci_high=float(ci.iloc[1]), p=float(fit.pvalues[target]), status='ok')
    return rec


def bootstrap_mean(values, rng, n=1000):
    v = np.asarray(values, dtype=float)
    b = np.array([rng.choice(v, len(v), replace=True).mean() for _ in range(n)])
    return np.quantile(b, [.025, .975])


def parse_feature(feature, celltypes):
    matches = [c for c in celltypes if str(feature).startswith(str(c) + '_')]
    if not matches:
        return None, None
    c = max(matches, key=len)
    return c, str(feature)[len(c) + 1:]


def fingerprint(cfg):
    return hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()


def environment_record():
    packages = ['numpy', 'pandas', 'scipy', 'statsmodels', 'anndata', 'h5py', 'pyarrow', 'matplotlib', 'PyYAML', 'pytest', 'nbformat']
    versions = {}
    for p in packages:
        try:
            versions[p] = importlib.metadata.version(p)
        except importlib.metadata.PackageNotFoundError:
            versions[p] = 'MISSING'
    return {'python': sys.version, 'executable': sys.executable, 'platform': platform.platform(), 'packages': versions, 'command': sys.argv}


def memory_guard(cfg, estimated_gb=0):
    # Enforce a soft RSS guard at stage boundaries; metadata remains O(n_cells).
    import resource
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
    if rss + estimated_gb > cfg.get('max_memory_gb', 16):
        raise MemoryError(f'RSS high-water mark {rss:.2f} GiB + estimate {estimated_gb:.2f} exceeds configured guard')


def stage_status(run, stage, lineage, status, reason='', **extra):
    LOG.info('%s %s: %s — %s', stage, lineage, status, reason)
    write_json(Path(run)/'checkpoints'/f'{lineage}_{stage}.json', {'stage': stage, 'lineage': lineage, 'status': status, 'reason': reason, **extra})


def load_cells(run, lineage):
    return pd.read_parquet(Path(run)/'cache'/f'{lineage}_cells.parquet')
