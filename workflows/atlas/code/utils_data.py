"""Read-only, bounded-memory H5AD access. Never materialize the atlas or its graph."""
from pathlib import Path
import hashlib
import json
import logging
import re

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
try:
    from anndata.io import read_elem
except ImportError:
    from anndata._io.specs import read_elem

MARKERS = {
    'fibro': {
        'qFibro': ['APOD', 'CFD', 'CLU', 'PI16'],
        'eFibro': ['COL1A1', 'COL1A2', 'COL6A3', 'LUM', 'DCN', 'COL3A1'],
        'iFibro': ['CXCL12', 'CCL2', 'ICAM1', 'IL33', 'CD74'],
        'Myfibro': ['ACTA2', 'MYL9', 'TPM2', 'TAGLN']},
    'endo': {
        'VenEC': ['CLU', 'VCAM1'], 'MCs': ['MYL9', 'TAGLN', 'TPM2'],
        'LEC': ['EFEMP1', 'PDPN'], 'EndoMT': ['VIM', 'HLA-DRA', 'HLA-DQA1'],
        'CapEC': ['TM4SF1', 'RGCC', 'SPARC'],
        'ArtEC': ['CXCL12', 'PTPRM', 'SRGN', 'ARL15']}}
ORDER = {'fibro': ['qFibro', 'eFibro', 'iFibro', 'Myfibro'],
         'endo': ['ArtEC', 'CapEC', 'VenEC', 'EndoMT', 'MCs', 'LEC']}
COLORS = {'qFibro': '#66A88D', 'eFibro': '#619BC2', 'iFibro': '#D7A254',
          'Myfibro': '#BC788E', 'ArtEC': '#AD7FA8', 'CapEC': '#65AC9E',
          'VenEC': '#7497C5', 'EndoMT': '#D19A62', 'MCs': '#BE7F91', 'LEC': '#A5B86F'}
ALIASES = {'celltype': ['Celltype', 'cell_type', 'celltype'],
           'subtype': ['Celltype_sub', 'celltype_sub', 'subtype'],
           'general': ['Celltype_general', 'celltype_general', 'broad_celltype'],
           'tissue': ['tissue', 'Tissue', 'tissue_1'],
           'sample': ['sample', 'sample_id', 'Sample', 'GSM'],
           'donor': ['donor', 'donor_id', 'patient_id', 'individual'],
           'dataset': ['GSE', 'dataset', 'study_id', 'batch'],
           'age_group': ['Age_group', 'age_group', 'age_bin'],
           'age': ['age', 'Age', 'age_years'], 'leiden': ['leiden', 'Leiden']}
AGE_LABELS = dict(zip('ABCDEFGH', ['A: <20*', 'B: 20–29', 'C: 30–39', 'D: 40–49',
                                           'E: 50–59', 'F: 60–69', 'G: 70–89', 'H: ≥90*']))


def dump_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False,
                                   default=lambda x: x.item() if isinstance(x, np.generic) else str(x)))


def read_frame(group, names=None):
    """Only load needed columns; historical objects include large text reasoning fields."""
    idx = read_elem(group[group.attrs['_index']])
    cols = names if names is not None else list(group.attrs['column-order'])
    return pd.DataFrame({k: read_elem(group[k]) for k in cols}, index=pd.Index(idx.astype(str)))


def clean(s):
    s = pd.Series(s, index=getattr(s, 'index', None)).astype('string').str.strip()
    return s.mask(s.str.lower().isin(['', 'nan', 'na', 'none', 'unknown', 'n/a']), 'Unknown').fillna('Unknown').astype('category')


def age_order(values):
    vals = list(pd.unique(values))
    if set(vals) <= set('ABCDEFGH') | {'Unknown'}:
        return [v for v in list('ABCDEFGH') + ['Unknown'] if v in vals]
    stages = ['Embryo', 'Infant', 'Children', 'Adolescent', 'Young Adult', 'Adult', 'Aged', 'Super old']
    if set(vals) <= set(stages) | {'Unknown'}:
        return [v for v in stages + ['Unknown'] if v in vals]
    def key(s):
        m = re.search(r'\d+', str(s))
        return (str(s) == 'Unknown', float(m.group()) if m else float('inf'), str(s))
    return sorted(vals, key=key)


def matrix_shape(node):
    return tuple(node.attrs['shape']) if isinstance(node, h5py.Group) else node.shape


def csr_rows(node, start, stop):
    """CSR reads are contiguous; CSC is rejected rather than secretly loading it whole."""
    if isinstance(node, h5py.Dataset):
        return sparse.csr_matrix(np.asarray(node[start:stop], dtype=np.float32))
    enc = node.attrs.get('encoding-type', '')
    if enc != 'csr_matrix':
        raise ValueError(f'Need CSR or dense row-chunked source, got {enc}; convert a COPY first.')
    ptr = node['indptr'][start:stop + 1]
    lo, hi = int(ptr[0]), int(ptr[-1])
    return sparse.csr_matrix((node['data'][lo:hi].astype(np.float32),
                              node['indices'][lo:hi], ptr - lo),
                             shape=(stop - start, matrix_shape(node)[1]))


def selected_dense_rows(node, rows, dims=None):
    """Avoid multi-million-row HDF5 fancy indexing; gather contiguous row blocks."""
    rows = np.asarray(rows)
    shape = matrix_shape(node)
    out = np.empty((len(rows), dims or shape[1]), np.float32)
    bins = rows // 4096
    for b in np.unique(bins):
        pos = np.flatnonzero(bins == b); lo = int(b * 4096)
        block = node[lo:min(lo + 4096, shape[0]), :out.shape[1]]
        out[pos] = block[rows[pos] - lo]
    if not np.isfinite(out).all():
        raise ValueError('Embedding contains non-finite values')
    return out


def stratified_rows(obs, n, seed, columns=('general', 'tissue')):
    """Reserve one per stratum, then sample remaining cells proportionally."""
    if n >= len(obs):
        return np.arange(len(obs))
    rng = np.random.default_rng(seed)
    keys = pd.MultiIndex.from_frame(obs[list(columns)])
    codes, uniques = pd.factorize(keys)
    perm = rng.permutation(len(obs))
    _, first = np.unique(codes[perm], return_index=True)
    reserve = perm[first]
    if len(reserve) > n:
        raise ValueError(f'{n} cells cannot represent {len(reserve)} strata; increase sampling size.')
    avail = np.ones(len(obs), bool); avail[reserve] = False
    extra = rng.choice(np.flatnonzero(avail), n - len(reserve), replace=False)
    return np.sort(np.r_[reserve, extra])


def load_metadata(path, lineage, config):
    with h5py.File(path, 'r') as f:
        cols = list(f['obs'].attrs['column-order'])
        mapping = {}
        override = config.get('fields', {}).get(lineage, {})
        for role, options in ALIASES.items():
            mapping[role] = override.get(role) or next((x for x in options if x in cols), None)
            if mapping[role] and mapping[role] not in cols:
                raise ValueError(f'{lineage}: configured column missing: {mapping[role]}')
        needed = list(dict.fromkeys(v for v in mapping.values() if v))
        original = read_frame(f['obs'], needed)
        obs = pd.DataFrame(index=original.index)
        for role, col in mapping.items():
            obs[role] = clean(original[col]) if col else pd.Categorical(['Unknown'] * len(original))
        if not mapping['celltype'] or not mapping['subtype']:
            raise ValueError('Celltype and Celltype_sub must be mapped explicitly when absent.')
        if not mapping['general']:
            obs['general'] = clean(obs['subtype'].astype(str).str.split('_').str[0])
            mapping['general'] = 'derived: subtype prefix, validated against expected states'
        invalid = set(obs['general'].unique()) - set(ORDER[lineage])
        if invalid:
            raise ValueError(f'Unrecognized broad states {invalid}; supply correct general field mapping.')
        if not mapping['age_group']:
            # Numeric ages only; ambiguous free text stays Unknown.
            age = pd.to_numeric(original[mapping['age']], errors='coerce') if mapping['age'] else pd.Series(np.nan, index=obs.index)
            obs['age_group'] = pd.cut(age, [0,20,30,40,50,60,70,90,np.inf], labels=list('ABCDEFGH'), right=False).astype('string').fillna('Unknown').astype('category')
            mapping['age_group'] = 'derived numeric age bins; unparsed free text remains Unknown'
        if 'X_umap' not in f['obsm']:
            raise ValueError('No existing X_umap; this atlas overview does not invent a new atlas embedding.')
        source = config['expression_source']
        if source not in f:
            raise ValueError(f'Expression source {source} missing')
        var = read_frame(f['raw/var'] if source.startswith('raw/') else f['var'])
        if not var.index.is_unique or not obs.index.is_unique:
            raise ValueError('Cell/gene names must be unique; resolve ambiguous identifiers explicitly.')
        shape = matrix_shape(f[source])
        if shape != (len(obs), len(var)):
            raise ValueError('Expression/annotation shapes differ')
        probe = csr_rows(f[source], 0, min(128, len(obs))).data
        if not np.isfinite(probe).all() or (probe < 0).any():
            raise ValueError('Expression has negatives/non-finite values; scaled/PCA data are unsuitable.')
        scale = config['expression_scale']
        if scale not in ['log1p', 'counts']:
            raise ValueError('expression_scale must be explicitly log1p or counts')
        if scale == 'counts' and not np.allclose(probe, np.round(probe), atol=1e-5):
            raise ValueError('Configured counts are not integer-like')
        if scale == 'log1p' and 'log1p' not in f['uns']:
            raise ValueError('No log1p provenance; verify input scale before using it as log1p expression.')
        palettes = {}
        for role in ['celltype','subtype','general','tissue']:
            key = f'{mapping[role]}_colors'
            if key in f['uns'] and mapping[role] in original:
                s = original[mapping[role]]
                if isinstance(s.dtype, pd.CategoricalDtype):
                    color = read_elem(f['uns'][key])
                    if len(color) == len(s.cat.categories):
                        palettes[role] = dict(zip(s.cat.categories.astype(str), color.astype(str)))
        details = dict(input=str(Path(path).resolve()), n_cells=len(obs), n_genes=len(var),
                       obs_columns=cols, var_columns=list(var.columns), obsm=list(f['obsm']),
                       layers=list(f['layers']), has_raw='raw' in f, uns_keys=list(f['uns']),
                       fields=mapping, expression_source=source, expression_scale=scale,
                       expression_probe_min=float(probe.min()) if len(probe) else None,
                       expression_probe_max=float(probe.max()) if len(probe) else None,
                       neighbors=read_elem(f['uns/neighbors']) if 'neighbors' in f['uns'] else {},
                       categories={k: obs[k].value_counts().to_dict() for k in ['general','subtype','celltype','tissue','age_group']},
                       n_samples=int(obs['sample'][obs['sample'] != 'Unknown'].nunique()),
                       n_donors=int(obs['donor'][obs['donor'] != 'Unknown'].nunique()) if mapping['donor'] else None,
                       n_tissues=int(obs['tissue'][obs['tissue'] != 'Unknown'].nunique()),
                       n_datasets=int(obs['dataset'][obs['dataset'] != 'Unknown'].nunique()),
                       missing={k: int((obs[k] == 'Unknown').sum()) for k in obs}, palettes=palettes)
    return obs, var, details


def audit(path, lineage, config, root):
    obs, var, info = load_metadata(path, lineage, config)
    shared = root / 'shared'; shared.mkdir(exist_ok=True, parents=True)
    dump_json(shared / f'data_audit_{lineage}.json', info)
    pd.DataFrame([{'item': k, 'value': json.dumps(v, ensure_ascii=False, default=str)} for k,v in info.items()]).to_csv(shared / f'data_audit_{lineage}.csv', index=False)
    out = root / lineage; out.mkdir(exist_ok=True)
    obs.to_pickle(out / 'metadata.pkl')
    var.to_pickle(out / 'var.pkl')
    for role in ['general','subtype','celltype','tissue','age_group']:
        obs[role].value_counts().rename_axis(role).rename('n_cells').to_csv(out / f'counts_{role}.csv')
    # Audit distinct sample annotation conflicts; never call samples donors.
    pairs = obs.loc[obs['sample'] != 'Unknown', ['sample','tissue','age_group','dataset']].drop_duplicates()
    pairs.to_csv(out / 'sample_annotation_combinations.csv', index=False)
    conflicts = pairs.groupby('sample', observed=True).agg({k:'nunique' for k in ['tissue','age_group','dataset']})
    conflicts[(conflicts > 1).any(axis=1)].to_csv(out / 'sample_annotation_conflicts.csv')
    with h5py.File(path, 'r') as f:
        # Compact embedding copy; no expression or full neighbor graph is copied.
        np.save(out / 'umap.npy', f['obsm/X_umap'][:].astype(np.float32))
    logging.info('%s audit: %s cells, %s genes, %s samples', lineage, len(obs), len(var), info['n_samples'])
    return info


def normalize_block(x, scale):
    if not np.isfinite(x.data).all() or (x.data < 0).any():
        raise ValueError('Invalid expression values in streamed block')
    if scale == 'counts':
        sums = np.asarray(x.sum(axis=1)).ravel()
        x = sparse.diags(1e4 / np.maximum(sums, 1)) @ x
        x.data = np.log1p(x.data)
    return x.tocsr()


def marker_signature(path, lineage, config, root):
    """One full sequential expression pass after a small pilot for control bins."""
    out = root / lineage
    obs = pd.read_pickle(out / 'metadata.pkl'); var = pd.read_pickle(out / 'var.pkl')
    n, p = len(obs), len(var); gene_index = {g:i for i,g in enumerate(var.index)}
    groups = ORDER[lineage]; marker_sets = {k:sorted(MARKERS[lineage][k]) for k in groups}
    available = {k:[g for g in v if g in gene_index] for k,v in marker_sets.items()}
    missing = {k:sorted(set(marker_sets[k])-set(available[k])) for k in groups}
    dump_json(out / 'marker_availability.json', dict(available=available, missing=missing))
    for k in groups:
        if not available[k]:
            raise ValueError(f'No markers for {k}; cannot compute its signature')
    marker_genes = [g for k in groups for g in available[k]]
    marker_cols = [gene_index[g] for g in marker_genes]
    traj_rows = stratified_rows(obs, config['trajectory_cells'], config['seed'], ('subtype','tissue'))
    np.save(out / 'trajectory_rows.npy', traj_rows)
    np.save(out / 'trajectory_marker_genes.npy', np.array(marker_genes))
    traj_x = np.empty((len(traj_rows), len(marker_genes)), np.float32)
    score_path = out / 'signature_scores.npy'
    scores = np.lib.format.open_memmap(score_path, mode='w+', dtype=np.float32, shape=(n,len(groups)))
    rng = np.random.default_rng(config['seed'])
    sums = np.zeros((len(groups),len(marker_genes)), np.float64)
    detected = np.zeros_like(sums); counts = np.zeros(len(groups), np.int64)
    score_sums = np.zeros((len(groups),len(groups)), np.float64)
    with h5py.File(path, 'r') as f:
        node = f[config['expression_source']]
        starts = np.unique(np.linspace(0, max(0,n-config['pilot_block_rows']), config['pilot_blocks'], dtype=int))
        pilot_sum = np.zeros(p, np.float64); pilot_n = 0
        for start in starts:
            x = normalize_block(csr_rows(node, int(start), min(n,int(start)+config['pilot_block_rows'])), config['expression_scale'])
            pilot_sum += np.asarray(x.sum(axis=0)).ravel(); pilot_n += x.shape[0]
        means = pilot_sum / pilot_n
        # Rank bins follow the expression-matched control principle of scanpy score_genes.
        bins = pd.qcut(pd.Series(means).rank(method='first'), min(config['signature_bins'],p), labels=False).to_numpy()
        controls = {}; all_marker_cols = set(marker_cols)
        for k in groups:
            occupied = np.unique(bins[[gene_index[g] for g in available[k]]]); chosen=[]
            for b in occupied:
                pool = np.array([i for i in np.flatnonzero(bins == b) if i not in all_marker_cols])
                if len(pool): chosen.extend(rng.choice(pool, min(len(pool),config['controls_per_bin']), replace=False).tolist())
            controls[k] = sorted(set(chosen))
            if not controls[k]: raise ValueError('No expression-matched control genes')
        dump_json(out / 'signature_controls.json', dict(method='mean(log expression of marker genes) minus mean(log expression of fixed matched controls); controls binned using deterministic evenly spaced pilot blocks; all supplied markers excluded from control pool',
                  pilot_starts=starts.tolist(),pilot_n=pilot_n, markers=available,
                  controls={k:var.index[v].tolist() for k,v in controls.items()},seed=config['seed']))
        pd.DataFrame({'gene':var.index,'pilot_mean':means,'expression_bin':bins}).to_csv(out / 'signature_control_bins.csv',index=False)
        selected = sorted(set(marker_cols) | {i for v in controls.values() for i in v})
        remap = {v:i for i,v in enumerate(selected)}
        score_indices = {k:([remap[gene_index[g]] for g in available[k]], [remap[i] for i in controls[k]]) for k in groups}
        code = pd.Categorical(obs['general'], categories=groups).codes
        chunk = config['chunk_rows']
        for start in range(0,n,chunk):
            stop = min(start+chunk,n)
            x = normalize_block(csr_rows(node,start,stop), config['expression_scale'])
            y = x[:,selected]
            for j,k in enumerate(groups):
                m,c = score_indices[k]
                scores[start:stop,j] = np.asarray(y[:,m].mean(axis=1)-y[:,c].mean(axis=1)).ravel()
            mx = x[:,marker_cols].toarray()
            for j,k in enumerate(groups):
                mask = code[start:stop] == j
                counts[j] += int(mask.sum())
                sums[j] += mx[mask].sum(axis=0,dtype=np.float64)
                detected[j] += (mx[mask] > 0).sum(axis=0)
                score_sums[j] += np.asarray(scores[start:stop])[mask].sum(axis=0,dtype=np.float64)
            lo,hi = np.searchsorted(traj_rows,[start,stop])
            traj_x[lo:hi] = mx[traj_rows[lo:hi]-start]
            if start == 0 or (start // chunk) % 50 == 0 or stop == n:
                logging.info('%s expression %s/%s (%.1f%%)',lineage,stop,n,100*stop/n)
        scores.flush()
    np.save(out / 'trajectory_markers.npy',traj_x)
    np.save(out / 'signature_names.npy',np.array(groups))
    rows=[]
    for i,k in enumerate(groups):
        for j,g in enumerate(marker_genes):
            rows.append(dict(subtype=k, gene=g, marker_set=next(s for s in groups if g in available[s]),
                             mean_expression=sums[i,j]/max(counts[i],1), fraction_expressing=detected[i,j]/max(counts[i],1), n_cells=int(counts[i])))
    pd.DataFrame(rows).to_csv(out / 'marker_summary.csv',index=False)
    pd.DataFrame(score_sums / np.maximum(counts[:,None],1), index=groups,columns=groups).to_csv(out / 'signature_mean.csv')
    # Scores are descriptive evidence using the same genes as the annotations; not independent validation.
    pred = np.asarray(scores).argmax(axis=1)
    confusion = pd.crosstab(pd.Categorical(obs['general'],categories=groups),pd.Categorical(np.array(groups)[pred],categories=groups),dropna=False)
    confusion.to_csv(out / 'signature_concordance_counts.csv')
    dump_json(out / 'signature_concordance.json',dict(top_score_agreement=float(np.mean(pred == code)),
              interpretation='Descriptive and marker-dependent; raw module scores with different gene sets are not a calibrated classifier.'))
    sample_rows=[]
    for j,k in enumerate(groups):
        frame=pd.DataFrame({'sample':obs['sample'].to_numpy(),'subtype':obs['general'].to_numpy(),'score':scores[:,j]})
        frame=frame[frame['sample'] != 'Unknown']
        t=frame.groupby(['sample','subtype'],observed=True)['score'].agg(['mean','size']).reset_index()
        t['signature']=k;sample_rows.append(t)
    pd.concat(sample_rows,ignore_index=True).to_csv(out / 'signature_sample_means.csv',index=False)
    logging.info('%s marker/signature complete',lineage)
