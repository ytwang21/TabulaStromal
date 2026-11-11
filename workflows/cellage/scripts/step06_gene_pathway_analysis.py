"""Audited stable Celltype_Gene features, local-GMT ORA and evidence tables."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
from scipy import stats
import yaml
from .checkpointing import cached_value
from .common import (fdr, load_cells, parse_feature, read_table, save_table,
                     stage_status, write_json)


def feature_table(path, celltypes):
    d = read_table(path)
    if 'Feature' not in d:
        raise ValueError('Feature universe requires actual tested Feature identifiers')
    parsed = [parse_feature(f, celltypes) for f in d.Feature]
    d['parsed_celltype'] = [a for a, b in parsed]
    d['parsed_gene'] = [b for a, b in parsed]
    if 'Celltype' in d:
        conflict = d.parsed_celltype.notna() & d.Celltype.astype(str).ne(d.parsed_celltype.astype(str))
        if conflict.any():
            raise ValueError('Provided Celltype disagrees with longest-prefix feature parser')
    if 'Gene' in d:
        conflict = d.parsed_gene.notna() & d.Gene.astype(str).ne(d.parsed_gene.astype(str))
        if conflict.any():
            raise ValueError('Provided Gene disagrees with feature parser')
    return d


def read_gmt(path):
    terms = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 3:
                terms[fields[0]] = set(fields[2:])
    return terms


def ora(selected, universe, terms, min_size=5, max_size=3000):
    u = set(universe)
    s = set(selected) & u
    out = []
    if len(s) < 2 or len(u) < 10:
        return out
    for term, genes in terms.items():
        tested = genes & u
        if not min_size <= len(tested) <= max_size:
            continue
        hit = tested & s
        # Include zero-overlap terms in the testing family; don't prefilter by p.
        p = stats.hypergeom.sf(len(hit)-1, len(u), len(tested), len(s))
        out.append({'term': term, 'overlap_n': len(hit), 'term_tested_n': len(tested),
                    'selected_n': len(s), 'universe_n': len(u), 'genes': '|'.join(sorted(hit)), 'p': p,
                    'fold_enrichment': (len(hit)/len(s))/(len(tested)/len(u))})
    return out


def run(cfg, run, lineage):
    run = Path(run)
    g = cfg['lineages'][lineage]['genemodel']
    if not Path(g['stable_features']).exists() or not Path(g['feature_universe']).exists():
        stage_status(run, 'genes', lineage, 'blocked', 'Missing stable-feature table or actual tested universe')
        return
    cells = load_cells(run, lineage)
    types = list(cells.celltype.dropna().unique())
    stable = feature_table(g['stable_features'], types)
    tested = feature_table(g['feature_universe'], types)
    save_table(stable.loc[stable.parsed_celltype.isna()], run/'audit'/f'{lineage}_unparsed_stable_features.csv')
    stable = stable.loc[stable.parsed_celltype.notna()].copy()
    tested = tested.loc[tested.parsed_celltype.notna()].copy()
    # Candidate panel membership is part of the actual tested feature universe.
    if 'In_candidate_panel' in tested:
        tested = tested.loc[tested.In_candidate_panel.astype(str).str.lower().eq('true')]
    if not set(stable.Feature) <= set(tested.Feature):
        raise ValueError('Stable features absent from tested candidate universe')
    stable['source'] = g['stable_features']
    stable['evidence_tier'] = 'exploratory; stable model feature does not establish a senescence mechanism'
    save_table(stable, run/'tables'/f'{lineage}_stable_features_parsed.csv')
    save_table(tested[['Feature', 'parsed_celltype', 'parsed_gene']], run/'tables'/f'{lineage}_tested_feature_universe.csv')
    # Transparent sign-weighted descriptive ranking; no causal pro-aging label.
    summary = stable.groupby('parsed_celltype').agg(n_stable=('Feature', 'size'), n_genes=('parsed_gene', 'nunique'),
        signed_coef_mean=('Mean_coefficient', 'mean'), mean_selection_frequency=('Selection_frequency', 'mean')).reset_index().rename(columns={'parsed_celltype':'celltype'})
    summary['definition'] = 'mean signed standardized EN coefficient; model-specific, correlated features; not causal'
    summary['source_model'] = 'earlier downstream EN stable-feature analysis; separate from K8000 rebenchmark clock'
    age_path = run/'tables'/f'{lineage}_age_associations.csv'
    if age_path.exists():
        age = pd.read_csv(age_path)
        age = age.loc[age.sensitivity.eq(f'min_cells_{cfg["min_cells"]}') & age.adjusted.eq(True) & age.score.eq('hUSI')]
        summary = summary.merge(age[['celltype', 'effect', 'q', 'n_donors', 'n_gse']], on='celltype', how='left')
        usable = summary.dropna(subset=['effect', 'signed_coef_mean'])
        rank_rho = stats.spearmanr(usable.effect, usable.signed_coef_mean).statistic if len(usable) >= 4 else None
        write_json(run/'results'/f'{lineage}_subtype_rank_comparison.json', {'n_subtypes': len(usable), 'descriptive_spearman': rank_rho,
                   'definition': 'age slope vs stable coefficient ranking; older model, not independent causal evidence'})
    save_table(summary.assign(lineage=lineage), run/'tables'/f'{lineage}_feature_subtype_integration.csv')
    # Detect overlap with model-based senescence scoring gene lists; association may be circular.
    overlap_rows = []
    for method, path in g.get('score_gene_lists', {}).items():
        if not Path(path).exists():
            continue
        weights = pd.read_csv(path, index_col=0)
        score_genes = set(weights.index.astype(str))
        stable_genes = set(stable.parsed_gene)
        shared = stable_genes & score_genes
        overlap_rows.append({'lineage': lineage, 'method': method, 'n_score_genes': len(score_genes), 'n_stable_genes': len(stable_genes),
                             'n_shared': len(shared), 'jaccard': len(shared)/len(stable_genes | score_genes),
                             'overlap_fraction_of_stable': len(shared)/len(stable_genes) if stable_genes else np.nan,
                             'shared_genes': '|'.join(sorted(shared)), 'source': path,
                             'no_overlap_sensitivity': 'Cannot recompute existing hUSI without changing original signature; non-overlap feature ORA provided separately'})
        stable[f'overlap_{method}'] = stable.parsed_gene.isin(score_genes)
    save_table(pd.DataFrame(overlap_rows), run/'tables'/f'{lineage}_score_feature_gene_overlap.csv')
    save_table(stable, run/'tables'/f'{lineage}_feature_evidence.csv')
    gene_cfg = yaml.safe_load(Path(g['gene_set_config']).read_text()) if Path(g['gene_set_config']).exists() else {}
    out, provenance = [], []
    for library, spec in gene_cfg.get('gene_sets', {}).items():
        path = Path(spec['path'])
        if not path.exists():
            continue
        terms = read_gmt(path)
        provenance.append({'library': library, 'path': str(path), 'source': spec.get('source'), 'version': spec.get('version'),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'local_mtime_utc': pd.Timestamp(path.stat().st_mtime, unit='s', tz='UTC').isoformat(),
            'download_date': 'See original download_manifest.json; file mtime is not asserted as download date'})
        for subtype, x in stable.groupby('parsed_celltype'):
            universe = set(tested.loc[tested.parsed_celltype.eq(subtype), 'parsed_gene'])
            for direction, selected in [('all', x), ('positive', x.loc[x.Mean_coefficient > 0]), ('negative', x.loc[x.Mean_coefficient < 0])]:
                branches = [('all_features', selected, universe)]
                if 'overlap_hUSI' in stable:
                    # This controls only feature overlap in ORA, not the senescence score itself.
                    nonshared = set(tested.loc[tested.parsed_celltype.eq(subtype), 'parsed_gene']) - score_genes
                    branches.append(('exclude_hUSI_signature_genes', selected.loc[~selected.overlap_hUSI], nonshared))
                for branch, selected_rows, u in branches:
                    for rec in cached_value(['ORA',library,subtype,direction,branch],lambda:ora(set(selected_rows.parsed_gene),u,terms)):
                        out.append({'lineage': lineage, 'celltype': subtype, 'direction': direction, 'branch': branch,
                                    'library': library, 'interpretation': 'stable-feature enrichment; not senescence-associated DE', **rec})
    enrichment = fdr(pd.DataFrame(out), ['lineage'])
    save_table(enrichment, run/'tables'/f'{lineage}_pathway_enrichment.csv')
    save_table(pd.DataFrame(provenance), run/'audit'/f'{lineage}_gene_set_provenance.csv')
    from .expression import run_expression
    run_expression(cfg, run, lineage, stable)
    stage_status(run, 'genes', lineage, 'complete' if provenance else 'conditional',
                 'Stable-feature integration; expression and true SenCID evidence remain separately gated', n_stable=len(stable), n_tested_terms=len(enrichment))
