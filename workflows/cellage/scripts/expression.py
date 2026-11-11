"""Optional, bounded-memory donor pseudobulk marker/score association analysis.

Counts provenance is an explicit prerequisite. Fractional, negative, or transformed
matrix values block the module. This does not silently round estimated counts.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
from .common import adjusted_fit, fdr, load_cells, save_table, stage_status

MARKERS = ['CDKN1A', 'CDKN2A', 'SERPINE1', 'GLB1', 'LMNB1']


def run_expression(cfg, run, lineage, stable):
    spec = cfg.get('expression', {})
    if not spec.get('enabled') or not spec.get('counts_provenance'):
        stage_status(run, 'expression', lineage, 'blocked', 'Expression analysis disabled: raw filename insufficient; need verified counts provenance. Marker absence does not exclude senescence.')
        return
    path = spec.get('counts_path', {}).get(lineage, cfg['lineages'][lineage]['raw_atlas'])
    a = ad.read_h5ad(path, backed='r')
    try:
        cells = load_cells(run, lineage)
        cells = cells.set_index('cell_key')
        if not a.obs_names.is_unique or not a.var_names.is_unique:
            raise ValueError('Expression IDs must be unique')
        eligible = cells.loc[cells.age_atlas_member & cells.donor_key.notna()]
        positions = a.obs_names.get_indexer(eligible.index)
        good = positions >= 0
        eligible = eligible.loc[good].copy()
        positions = positions[good]
        keys = ['donor_key', 'tissue', 'celltype']
        group_index = pd.MultiIndex.from_frame(eligible[keys])
        codes, groups = pd.factorize(group_index)
        genes = sorted(set(MARKERS) | set(stable.parsed_gene))
        genes = [g for g in genes if g in a.var_names]
        gi = a.var_names.get_indexer(genes)
        estimated = len(groups) * len(genes) * 8 / 1024**3
        if estimated > cfg.get('max_memory_gb', 16)/3:
            raise MemoryError('Pseudobulk array exceeds memory budget; restrict candidate genes explicitly')
        sums = np.zeros((len(groups), len(genes)), dtype=float)
        library = np.zeros(len(groups))
        n_cells = np.zeros(len(groups), dtype=int)
        order = np.argsort(positions)
        chunk = spec.get('chunk_size', 2048)
        for start in range(0, len(order), chunk):
            ix = order[start:start+chunk]
            block = a.X[positions[ix], :]
            values = block.data if sparse.issparse(block) else np.asarray(block).ravel()
            if np.any(~np.isfinite(values)) or np.any(values < 0) or np.any(np.abs(values-np.round(values)) > 1e-6):
                stage_status(run, 'expression', lineage, 'blocked', 'Verified count requirement failed: noninteger/negative/nonfinite values; no DE or marker results emitted')
                return
            depth = np.asarray(block.sum(axis=1)).ravel()
            small = block[:, gi].toarray() if sparse.issparse(block) else np.asarray(block)[:, gi]
            np.add.at(sums, codes[ix], small)
            np.add.at(library, codes[ix], depth)
            np.add.at(n_cells, codes[ix], 1)
        meta = pd.DataFrame(groups.tolist(), columns=keys)
        meta['n_cells'] = n_cells
        meta['library_counts'] = library
        logcpm = np.log2(1 + sums / np.maximum(library[:, None], 1) * 1e6)
        score_path = Path(run)/'tables'/f'{lineage}_age_donor_scores.parquet'
        if not score_path.exists():
            stage_status(run, 'expression', lineage, 'blocked', 'No donor score table')
            return
        scores = pd.read_parquet(score_path)
        scores = scores.loc[scores.score.eq('hUSI')]
        merged = meta.reset_index(names='row').merge(scores, on=keys, suffixes=('_expression', ''), validate='one_to_one')
        rows = []
        for subtype, x in merged.groupby('celltype', observed=True):
            for j, gene in enumerate(genes):
                if np.count_nonzero(sums[x.row, j] > 0) < spec.get('min_detected_donors', 6):
                    continue
                y = x.copy()
                y['logCPM'] = logcpm[y.row, j]
                rec = adjusted_fit(y, outcome='logCPM', target='score_mean', extra=['age'], min_donors=cfg['min_donors'])
                rows.append({'lineage': lineage, 'celltype': subtype, 'gene': gene, 'canonical_marker': gene in MARKERS,
                             'analysis': 'exploratory donor pseudobulk logCPM score association adjusted for age', **rec})
        save_table(fdr(pd.DataFrame(rows), []), Path(run)/'tables'/f'{lineage}_pseudobulk_score_gene_associations.csv')
        save_table(meta, Path(run)/'tables'/f'{lineage}_pseudobulk_coverage.csv')
        stage_status(run, 'expression', lineage, 'complete', n_tests=len(rows), counts_provenance=spec['counts_provenance'])
    finally:
        a.file.close()
