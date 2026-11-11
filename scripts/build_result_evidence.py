#!/usr/bin/env python3
"""Recheck result numerics and publish a scalar/feature evidence ledger.

Reads private prediction/sample files only for verification; exports no rows or
identifiers from them. Uses HDF5 observation metadata without loading expression.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from audit_source_results import SOURCE_ROOTS, write_csv
from publish_result_assets import ATLAS, BENCH, CELL

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, required=True)
    project = parser.parse_args().project_root.resolve()
    evidence, checks, facts = [], [], {}
    source_reads = {}

    def record_source(path, status='content_read_for_evidence; no individual rows exported'):
        path=Path(path)
        original=project/path
        if original.is_file():
            st=original.stat()
            from datetime import datetime,timezone
            source_reads[str(path)]=dict(module=next((k for k,v in SOURCE_ROOTS.items() if str(path).startswith(v+'/')),'shared'),
                source_file=str(path),type=original.suffix,bytes=st.st_size,
                mtime_utc=datetime.fromtimestamp(st.st_mtime,timezone.utc).isoformat(),expected_use='source verification',
                read_status=status,repository_file='',source_sha256='not hashed: metadata-only large input' if original.suffix=='.h5ad' else hashlib.sha256(original.read_bytes()).hexdigest())

    def read(path):
        record_source(path)
        return pd.read_csv(project / path, keep_default_na=False, low_memory=False)

    def check(name, condition, detail=''):
        checks.append(dict(check=name, passed=bool(condition), detail=str(detail)))
        if not condition:
            raise AssertionError(name + ': ' + str(detail))

    def ev(module, section, lin, name, value, source, locator, figure, panel='', **extra):
        record = dict(claim_id=f'R{len(evidence)+1:05d}', section=section, module=module,
                      lineage=lin, dataset='', result_name=name, exact_value=str(value),
                      metric_definition='', unit='', n_donors='not independently verified',
                      n_samples='', n_cells='', split='', K='', model='', fold_or_cv_protocol='',
                      adaptation_status='', stats_method='', p_value='', FDR='',
                      source_file=str(source), source_sheet_or_column=locator,
                      figure_file=figure, figure_panel=panel, evidence_status='verified_from_source',
                      notes='')
        record.update(extra)
        evidence.append(record)
        return record['claim_id']

    def score(y, p):
        y, p = np.asarray(y, float), np.asarray(p, float)
        return dict(n=len(y), R2=1-np.sum((y-p)**2)/np.sum((y-y.mean())**2),
                    MAE=np.mean(np.abs(y-p)), RMSE=np.sqrt(np.mean((y-p)**2)),
                    MedAE=np.median(np.abs(y-p)), MAD=np.mean(np.abs(y-y.mean())),
                    PCC=np.corrcoef(y, p)[0, 1])

    def assert_metrics(name, expected, actual):
        for key in ['n', 'R2', 'MAE', 'RMSE', 'MedAE', 'MAD', 'PCC']:
            check(name + ' ' + key, np.isclose(float(expected[key]), actual[key], rtol=1e-8, atol=1e-8))

    atlas_facts = {}
    for lin, lower in [('Fibro', 'fibro'), ('Endo', 'endo')]:
        audit_source = ATLAS / 'shared' / f'data_audit_{lower}.json'
        audited = json.loads((project / audit_source).read_text())
        record_source(audit_source)
        h5 = Path(audited['input'])
        record_source(h5.relative_to(project),'HDF5 observation metadata only; expression not read')
        counts = {}
        with h5py.File(h5, 'r') as handle:
            group = handle['obs']
            n = len(group[group.attrs['_index']])
            for col in ['sample', 'GSE', 'tissue', 'Celltype', 'Celltype_general']:
                node = group[col]
                if isinstance(node, h5py.Group) and 'codes' in node:
                    valid = set()
                    for start in range(0, n, 250000):
                        valid.update(int(i) for i in np.unique(node['codes'][start:start+250000]) if i >= 0)
                    counts[col] = len(valid)
                else:
                    valid = set()
                    for start in range(0, n, 250000):
                        valid.update(node[start:start+250000].astype(str))
                    counts[col] = len(valid)
            donor_fields = [str(c) for c in group.attrs['column-order']
                            if any(token in str(c).lower() for token in ['donor', 'patient', 'individual', 'subject'])]
        check(lin + ' atlas source cell count', n == audited['n_cells'])
        check(lin + ' atlas sample categories', counts['sample'] == audited['n_samples'])
        check(lin + ' atlas GSE categories', counts['GSE'] == audited['n_datasets'])
        check(lin + ' atlas tissue categories', counts['tissue'] == audited['n_tissues'])
        check(lin + ' no donor identity field', not donor_fields)
        atlas_facts[lin] = dict(n_cells=n, n_sample_labels=counts['sample'], n_GSE_labels=counts['GSE'],
                                n_raw_tissue_labels=counts['tissue'], n_Celltypes=counts['Celltype'],
                                n_broad_states=counts['Celltype_general'])
        for key, value in atlas_facts[lin].items():
            ev('F1_StromalAtlas', 'Results1', lin, key, value, h5.relative_to(project),
               f'HDF5 obs metadata; observed valid codes; {key}', 'Figure 1; Extended Data Figure 1',
               n_cells=n, n_samples=counts['sample'], unit='count')
        suffix_source = Path('F1_StromalAtlas/results/fig1_20261103_v2') / f'{lower}_Tissue_from_Celltype_summary.csv'
        ev('F1_StromalAtlas', 'Results1', lin, 'n_annotation_suffixes', len(read(suffix_source)),
           suffix_source, 'number of rows including pan', 'Figure 1', 'C-E', unit='annotation suffix categories')
        ss = Path('F1_StromalAtlas/extended_figure1_SS/data') / f'marker_combination_train_holdout_{lower}.csv'
        q = read(ss)
        row = q[(q.split == 'holdout') & (q.marker_version == 'selected_global_combination')].iloc[0]
        for key in ['GSE_equal_top1', 'n_GSE', 'GSE_bootstrap95_low', 'GSE_bootstrap95_high']:
            ev('F1_StromalAtlas', 'Results1', lin, 'signature_' + key, row[key], ss,
               f"split=holdout; marker_version=selected_global_combination; column={key}", 'Extended Data Figure 2',
               notes='Retrospective study-held-out agreement, not a newly blinded identity test')
        clusters_source = Path('F1_StromalAtlas/extended_figure1_SS/data/plot_reference_cluster_check.csv')
        cl = read(clusters_source)
        cl = cl[cl.lineage == lower]
        ev('F1_StromalAtlas', 'Results1', lin, 'cluster_signature_agreement',
           f'{int(cl.annotation_matches_highest.astype(str).str.lower().eq("true").sum())}/{len(cl)}',
           clusters_source, f'lineage={lower}; annotation_matches_highest', 'Extended Data Figure 2',
           notes='Retrospective agreement with the marker-informed reference labels')
    facts['atlas'] = atlas_facts

    summary = read(BENCH / 'final_clock_results_summary.csv')
    stable_sets = {}
    for lin, rel in [('Fibro', 'Fibro_Accelerated_20261102/results/Fibro'), ('Endo', 'Endo')]:
        root = BENCH / rel
        module = 'TabulaStromal_Rebenchmark_2026'
        meta = read(root / 'features/train/samples.csv')
        s = summary[summary.lineage == lin].iloc[0]
        check(lin + ' development samples', len(meta) == int(s.n_development_samples))
        check(lin + ' development groups', meta.group_id.nunique() == int(s.n_development_groups))
        splits = read(root / 'benchmark/cv_splits.csv')
        for fold, q in splits.groupby('fold'):
            check(f'{lin} CV fold{fold} group separation',
                  set(q[q.split == 'Train'].group_id).isdisjoint(set(q[q.split == 'Valid'].group_id)))
        for key in ['K', 'alpha', 'l1_ratio', 'n_bootstrap', 'threshold', 'n_stable_Celltype_Gene',
                    'n_stable_Celltypes', 'n_unique_stable_Genes', 'n_development_samples', 'n_development_groups',
                    'internal_original_R2', 'internal_stable_R2', 'internal_original_MAE', 'internal_stable_MAE']:
            ev(module, 'Results3' if 'stable' in key or key in ['alpha','l1_ratio','n_bootstrap','threshold'] else 'Results2',
               lin, key, s[key], BENCH / 'final_clock_results_summary.csv', f'lineage={lin}; column={key}',
               'Figure 3', K=int(s.K), model='ElasticNet', n_samples=len(meta))
        folds = read(root / 'benchmark/cv_fold_metrics.csv')
        means = read(root / 'benchmark/cv_summary.csv')
        check(lin + ' all benchmark tasks', len(folds) == 520 and len(means) == 104)
        for key,value in [('n_models',folds.model.nunique()), ('requested_K_values',';'.join(map(str,sorted(folds.K.unique())))),
                          ('n_outer_folds',folds.fold.nunique())]:
            ev(module,'Results2',lin,key,value,root/'benchmark/cv_fold_metrics.csv',
               f'observed distinct model/K/fold values; {key}','Figure 2; complete benchmark tables')
        bs_source = root / 'stability/bootstrap_runs_manifest.csv'
        record_source(bs_source,'scope/bootstrap/status columns read; group identifiers not exported')
        bs = pd.read_csv(project / bs_source, usecols=['scope','bootstrap','status'])
        check(lin+' actual bootstrap manifest records', len(bs)==1800 and bs.scope.nunique()==6 and
              bs.status.eq('success').all() and not bs.duplicated(['scope','bootstrap']).any() and bs.groupby('scope').size().eq(300).all())
        ev(module,'Results3',lin,'successful_bootstraps_by_scope',
           ';'.join(f'{k}:{v}' for k,v in bs.groupby('scope').size().items()),bs_source,
           'count records by scope; all status=success; unique scope/bootstrap', 'Figure 3','C,D',
           notes='1800 actual bootstrap records across 6 scopes; 300 each')
        for _, row in means.iterrows():
            q = folds[(folds.K == row.K) & (folds.model == row.model) & (folds.split == row.split)]
            for key in ['R2','MAE','MAD','MedAE','RMSE','PCC']:
                check(f'{lin} {row.model} K{row.K} {row.split} {key} mean', np.isclose(q[key].mean(), row[key+'_mean']))
                check(f'{lin} {row.model} K{row.K} {row.split} {key} SD', np.isclose(q[key].std(ddof=1), row[key+'_std']))
        winner = means[means.split == 'Valid'].sort_values(['R2_mean','MAE_mean'], ascending=[False, True]).iloc[0]
        facts[lin] = dict(winner=winner.to_dict(), summary=s.drop('final_model_path').to_dict())
        cvs = read(root / 'stability/final_clock_internal_cv_summary.csv')
        for _, c in cvs.iterrows():
            for metric in ['R2_mean','R2_std','MAE_mean','MAE_std']:
                ev(module,'Results2',lin,'conditional_EN_'+metric,c[metric],root/'stability/final_clock_internal_cv_summary.csv',
                   f'model_version={c.model_version}; column={metric}','Figure 3','A,B',
                   model=c.model_version,K=int(s.K),fold_or_cv_protocol='five outer folds, conditional on selected K')
        for key in ['K','R2_mean','R2_std','MAE_mean','MAE_std','PCC_mean','PCC_std']:
            ev(module, 'Results2', lin, 'global_winner_'+key, winner[key], root / 'benchmark/cv_summary.csv',
               f'model={winner.model}; K={winner.K}; split=Valid; column={key}',
               'Figure 2 provides source-fold example; fivefold table is primary numeric evidence',
               model=winner.model, K=int(winner.K), split='Valid', fold_or_cv_protocol='grouped fivefold mean ± sample SD')
        original = read(root / 'external/final_clock_external_combined_predictions.csv')
        source_metrics = read(root / 'external/final_clock_external_combined_metrics.csv').iloc[0]
        assert_metrics(lin + ' noadapt direct prediction recomputation', source_metrics, score(original.age_true, original.age_pred))
        split = read(root / 'external/external_ft30_split.csv')
        adapt_ids = set(split.loc[split.set == 'adaptation', 'sample_id'])
        hold_ids = set(split.loc[split.set == 'holdout', 'sample_id'])
        check(lin + ' external adaptation separation', not adapt_ids & hold_ids)
        check(lin + ' external group separation',
              set(split[split.set == 'adaptation'].group_id).isdisjoint(set(split[split.set == 'holdout'].group_id)))
        post = read(root / 'external/external_ft30_all_100_predictions.csv')
        protocol = read(root / 'external/final_clock_external_protocol_comparison.csv')
        for idx, row in protocol.iterrows():
            table = post if row.model_version != 'final_clock_noadapt' else original
            q = table[table.sample_id.isin(hold_ids)] if row.evaluation_source == 'external_ft30_heldout70' else table
            assert_metrics(f'{lin} external protocol row{idx}', row, score(q.age_true, q.age_pred))
            for key in ['n', 'n_adapt', 'n_heldout', 'R2', 'MAE', 'RMSE', 'PCC', 'MedAE', 'MAD']:
                ev(module, 'Results2', lin, 'external_'+key, row[key], root / 'external/final_clock_external_protocol_comparison.csv',
                   f'row={idx+1}; model_version={row.model_version}; evaluation_source={row.evaluation_source}; column={key}',
                   'Figure 3', 'M,N', model='ElasticNet', K=int(s.K), split=row.evaluation_source,
                   adaptation_status=row.model_version, n_samples=int(row.n),
                   metric_definition='MAD=mean(abs(age-mean(age))); MedAE=median prediction error' if key in ['MAD','MedAE'] else '')
        st = read(root / 'stability/stable_Celltype_Gene.csv')
        check(lin + ' stable count', len(st) == int(s.n_stable_Celltype_Gene))
        check(lin + ' stable frequency threshold', (st.selection_frequency >= .7).all())
        stable_sets[lin] = set(st.Gene)
        for _, row in st.iterrows():
            ev(module, 'Results3', lin, row.feature, row.selection_frequency, root / 'stability/stable_Celltype_Gene.csv',
               f'feature={row.feature}; column=selection_frequency', 'Figure 3', 'G-L',
               model='ElasticNet', K=int(s.K), metric_definition='selection count / 300 successful group bootstrap fits',
               notes=f'Celltype={row.Celltype}; Gene={row.Gene}; coefficient_mean={row.coefficient_mean}; sign_consistency={row.sign_consistency}')
        sub = read(root / 'stability/stable_celltypes_ranked.csv')
        for _, row in sub.iterrows():
            ev(module, 'Results3', lin, row.Celltype, row.n_stable_features, root / 'stability/stable_celltypes_ranked.csv',
               f'Celltype={row.Celltype}; column=n_stable_features', 'Figure 3', 'I,J',
               n_samples=row.n_samples, n_cells=row.n_cells, notes='Raw count; unequal candidate opportunity and coverage retained')
        pathways = read(root / 'enrichment/final_clock_pathway_enrichment.csv')
        global_rows = pathways[pathways.scope == 'global']
        for _, row in global_rows[pd.to_numeric(global_rows.FDR) < .05].iterrows():
            ev(module, 'Results3', lin, row.term, row.overlap_count, root / 'enrichment/final_clock_pathway_enrichment.csv',
               f'scope=global; library={row.library}; term={row.term}', 'Figure 3', 'O,P',
               p_value=row.p_value, FDR=row.FDR, stats_method='hypergeometric ORA; BH within scope/library incl zero-overlap terms',
               notes=f'query_size={row.query_size}; universe_size={row.universe_size}; set_size={row.set_size}; library={row.library}')
        assoc = read(root / 'enrichment/final_clock_stable_gene_age_association.csv')
        n_assoc = int((pd.to_numeric(assoc.FDR, errors='coerce') < .05).sum())
        ev(module, 'Results3', lin, 'FDR_significant_expression_age_pairs', n_assoc,
           root / 'enrichment/final_clock_stable_gene_age_association.csv', 'count FDR<0.05', 'Figure 3 supporting associations',
           notes='Covariate-adjusted association, not independent replication')
        facts[lin]['n_age_associations'] = n_assoc
        facts[lin]['global_significant_terms'] = global_rows[pd.to_numeric(global_rows.FDR) < .05].groupby('library').size().to_dict()
        sig=global_rows[pd.to_numeric(global_rows.FDR)<.05]
        for key,value in [('n_significant_global_primary_terms',int((sig.library!='MSigDB_Hallmark_2020').sum())),
                          ('n_significant_global_terms_with_Hallmark',len(sig))]:
            ev(module,'Results3',lin,key,value,root/'enrichment/final_clock_pathway_enrichment.csv',
               'scope=global; FDR<0.05; specified library subset','Figure 3','O,P',
               notes='Number of annotated terms, not independent mechanisms')
        facts[lin]['subtype_counts'] = dict(zip(sub.Celltype, sub.n_stable_features.astype(int)))
    shared = sorted(stable_sets['Fibro'] & stable_sets['Endo'])
    check('shared final stable genes', len(shared) == 9)
    facts['shared_stable_genes'] = shared
    ev('TabulaStromal_Rebenchmark_2026','Results3','both','shared_stable_genes',';'.join(shared),
       ';'.join(str(BENCH/p/'stability/stable_Celltype_Gene.csv') for p in ['Fibro_Accelerated_20261102/results/Fibro','Endo']),
       'set intersection of Gene columns','Figure 3 shared-gene supplement')
    sets_source=BENCH/'resources/provided_aging_sets.json'
    sets=json.loads((project/sets_source).read_text())
    record_source(sets_source)
    ev('TabulaStromal_Rebenchmark_2026','Results3','both','n_provided_unverified_aging_lists',len(sets['sets']),
       sets_source,'number of supplied sets; original source remains unverified','Figure 3','Q,R')

    cohort = read(CELL / 'tables/cohort_overview.csv')
    facts['cellage'] = {}
    for lin in ['fibro', 'endo']:
        row = cohort[cohort.lineage == lin].iloc[0]
        for key in ['hUSI_cells', 'recSID_cells', 'hUSI_recSID_common_cells']:
            ev('CellAge_Exploration', 'Results2', lin, key, row[key], CELL / 'tables/cohort_overview.csv',
               f'lineage={lin}; column={key}', 'Figure 4', 'C-H', unit='cells',
               notes='Run 20261030_220917_407544; source donor label is sample-based')
        age = read(CELL / f'tables/{lin}_age_associations.csv')
        q = age[(age.score == 'hUSI') & (age.sensitivity == 'min_cells_20') &
                (age.adjusted.astype(str).str.lower() == 'true') & (age.status == 'ok')]
        n_sig = int((pd.to_numeric(q.q, errors='coerce') < .05).sum())
        check(lin+' no adjusted age discoveries', n_sig == 0)
        facts['cellage'][lin] = dict(n_identifiable=len(q), n_significant=n_sig,
                                     hUSI_cells=int(row.hUSI_cells), recSID_cells=int(row.recSID_cells))
        contrasts = read(CELL / f'tables/{lin}_subtype_contrasts.csv')
        targets = ['iFibro_CD74_pan','eFibro_COL1A1_pan'] if lin == 'fibro' else ['EndoMT_PTN_pan','EndoMT_TAGLN_pan']
        for _, c in contrasts[contrasts.celltype.isin(targets) & contrasts.score.eq('hUSI')].iterrows():
            ev('CellAge_Exploration','Results2',lin,c.celltype+'_paired_hUSI_contrast',c.effect,
               CELL/f'tables/{lin}_subtype_contrasts.csv',f'score=hUSI; celltype={c.celltype}',
               'Figure 4','C,D',p_value=c.p,FDR=c.q,n_samples=c.n_donors,
               stats_method='source paired subtype contrast; bootstrap CI; source BH correction',
               notes=f'95% CI {c.ci_low} to {c.ci_high}; sample-designated units')
        for key, value in [('identifiable_adjusted_hUSI_age_models', len(q)), ('adjusted_hUSI_age_q_lt_0.05', n_sig)]:
            ev('CellAge_Exploration', 'Results2', lin, key, value, CELL / f'tables/{lin}_age_associations.csv',
               'score=hUSI; sensitivity=min_cells_20; adjusted=True; status=ok; BH q<0.05', 'Figure 4', 'E,F',
               stats_method='source adjusted age model; source BH family retained')
        for _, record in q.iterrows():
            ev('CellAge_Exploration', 'Results2', lin, record.celltype+'_adjusted_hUSI_age', record.effect,
               CELL / f'tables/{lin}_age_associations.csv', f'primary adjusted row; Celltype={record.celltype}',
               'Figure 4','E,F', p_value=record.p, FDR=record.q, n_samples=record.n_donors,
               notes=f'95% CI {record.ci_low} to {record.ci_high}; n_donors is a source-designated sample unit')
    write_csv(ROOT / 'evidence/evidence_registry.csv', evidence, list(evidence[0]))
    for module in SOURCE_ROOTS:
        rows = [r for r in evidence if r['module'] == module]
        if module == 'F2_ModelPerformance':
            data = read(Path(module) / 'data/Fig2_scatter_metrics.csv')
            for _, r in data.iterrows():
                ev(module, 'Results2', r.lineage, 'source_fold4_'+r['split'],
                   f'R2={r.R2}; MAE={r.MAE}; n={r.n}', Path(module)/'data/Fig2_scatter_metrics.csv',
                   f'lineage={r.lineage}; split={r["split"]}', 'Figure 2','E,F',
                   model='ElasticNet',K=20000,fold_or_cv_protocol='source fold 4',
                   adaptation_status='FT30 ALL100 for External, including fitted samples')
            rows = [r for r in evidence if r['module'] == module]
        write_csv(ROOT / module / 'evidence/evidence_registry.csv', rows, list(evidence[0]))
    write_csv(ROOT / 'evidence/evidence_registry.csv', evidence, list(evidence[0]))
    write_csv(ROOT / 'audit/numeric_and_statistical_checks.csv', checks, list(checks[0]))
    inv=ROOT/'audit/file_inventory.csv'
    with inv.open() as f:existing=list(csv.DictReader(f))
    by_source={r['source_file']:r for r in source_reads.values()}
    by_source.update({r['source_file']:r for r in existing})
    write_csv(inv,list(by_source.values()),list(existing[0]))
    for module in SOURCE_ROOTS:
        rows=[r for r in by_source.values() if r['module']==module]
        write_csv(ROOT/module/'audit/file_inventory.csv',rows,list(existing[0]))
    (ROOT / 'audit/verified_facts.json').write_text(json.dumps(facts, indent=2, default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
    print(json.dumps(dict(status='passed', numeric_checks=len(checks), evidence_records=len(evidence)),indent=2))


if __name__ == '__main__':
    main()
