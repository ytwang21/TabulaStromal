#!/usr/bin/env python3
"""Generate the three Results sections from released aggregate evidence only."""
import csv
import json
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/results'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    facts = json.loads((ROOT / 'audit/verified_facts.json').read_text())
    with (ROOT / 'evidence/evidence_registry.csv').open() as f:
        ledger = list(csv.DictReader(f))
    claims = []

    def table(path):
        return pd.read_csv(ROOT / path)

    def section(filename, title, blocks):
        lines = ['# ' + title, '']
        for heading, text, names, source_section in blocks:
            lines += ['## ' + heading, '', text.strip(), '']
            relevant = [r for r in ledger if r['section'] == source_section and
                        (not names or any(n in r['result_name'] for n in names))]
            assert relevant, (filename, heading)
            for number, paragraph in enumerate(re.split(r'\n\s*\n', text.strip()), 1):
                if paragraph.startswith('|'):
                    continue
                claims.append(dict(section_file='docs/results/' + filename, heading=heading,
                                   paragraph=number, claim_text=paragraph,
                                   claim_ids=';'.join(r['claim_id'] for r in relevant),
                                   source_files=';'.join(sorted({r['source_file'] for r in relevant})),
                                   status='mapped to source-verified evidence; protocols retained'))
        lines += ['The source-to-claim mapping is available in [the evidence registry](../../evidence/evidence_registry.csv) and [the paragraph evidence report](../../audit/claim_to_evidence_report.csv).', '']
        (OUT / filename).write_text('\n'.join(lines))

    atlas = facts['atlas']
    ss = {}
    for lin in ['fibro','endo']:
        d = table(f'F1_StromalAtlas/tables/signatures/marker_combination_train_holdout_{lin}.csv')
        ss[lin] = d[(d.split == 'holdout') & (d.marker_version == 'selected_global_combination')].iloc[0]
    agreement = ' '.join(
        f'{"Fibroblast" if lin=="fibro" else "Endothelial"} GSE-equal top-signature concordance was '
        f'{100*r.GSE_equal_top1:.2f}% (95% study-bootstrap CI {100*r.GSE_bootstrap95_low:.2f}–{100*r.GSE_bootstrap95_high:.2f}%; '
        f'{int(r.n_GSE)} retained studies).'
        for lin,r in ss.items())
    section('Results_Part1_StromalAtlas.md',
            'Results Part 1 — A human stromal atlas resolves cellular heterogeneity across studies', [
        ('Atlas scale and metadata coverage',
         f'''The integrated collections contained {atlas['Fibro']['n_cells']:,} fibroblast-lineage and {atlas['Endo']['n_cells']:,} endothelial-lineage cells, associated with {atlas['Fibro']['n_sample_labels']:,} and {atlas['Endo']['n_sample_labels']:,} sample labels and {atlas['Fibro']['n_GSE_labels']} and {atlas['Endo']['n_GSE_labels']} GSE study labels, respectively ([Fig. 1A,B](../../F1_StromalAtlas/figures/Fig1_StromalAtlas_main.pdf)). These are lineage-specific counts; overlapping contributions across lineages were not treated as additional independent participants. The source objects did not provide a harmonized donor identifier. Sample labels therefore describe the available sampling units, rather than a verified number of unique individuals.

The full-atlas tissue metadata contained 39 observed labels per lineage, including nonstandard categories such as ascites and carcinoma. In contrast, the Tissue_from_Celltype variable used in Fig. 1C–E comprised 12 fibroblast and 11 endothelial annotation suffixes, including pan. These suffixes describe annotation naming and are not interchangeable with anatomical source metadata. The full atlas also contains heterogeneous disease annotations and cannot be described as uniformly healthy. Age-clock development was performed in separately filtered cohorts.''',
         ['n_cells','n_sample_labels','n_GSE_labels','n_raw_tissue_labels','n_annotation_suffixes'],'Results1'),
        ('Hierarchical annotation and full-source expression summaries',
         '''Fibroblast annotations comprised Myfibro, eFibro, iFibro and qFibro, resolving into 32 detailed Celltype labels. The endothelial-lineage collection retained ArtEC, CapEC, VenEC, EndoMT, MCs and LEC, resolving into 38 detailed labels ([Extended Data Fig. 1](../../F1_StromalAtlas/figures/ExtendFig1_subtype_umap_composition.pdf)). MCs is a source-defined category within that collection; the retained label does not imply that every included cell is a conventional endothelial cell.

The existing UMAP coordinates and broad-state marker patterns described cellular heterogeneity without generating a new embedding. Marker panels summarize mean stored expression and detection fractions; signature panels show full-source state means, with separate quantitative and annotation color scales (Fig. 1G–J). These summaries support descriptive annotation and share the reference data used for labeling; they are not an independent cell-identity validation. Composition summaries retained unknown ages and described pooled-cell distributions. Cell-type proportions were not used as age-clock features.''',
         ['n_Celltypes','n_broad_states'],'Results1'),
        ('Retrospective signature agreement and exploratory resources',
         f'''A separate marker-combination analysis evaluated agreement with the existing broad reference annotations ([Extended Data Fig. 2](../../F1_StromalAtlas/figures/ExtFig1_SignatureScore_Fibro_Endo.pdf)). {agreement} At the annotation-by-Leiden level, the highest signature agreed with the reference label in 51/55 fibroblast and 32/44 endothelial records. Discordant records were retained in the figure.

These measurements are retrospective annotation agreement: reference labels were marker-informed and the retained studies had been inspected during earlier work. They are not a new blinded accuracy test. The stored filtering-strategy resources and metacell trajectories provide additional descriptive context. In particular, assumed trajectory roots and study-dependent transcriptomic structure do not establish a temporal aging sequence or experimentally observed lineage fate.''',
         ['signature_','cluster_signature_agreement'],'Results1'),
    ])

    winners = []
    final = []
    protocol_rows = []
    examples = []
    for lin,label in [('Fibro','fibroblast'),('Endo','endothelial')]:
        w = facts[lin]['winner']
        winners.append(f"For the {label} collection, {w['model']} at K={int(w['K']):,} ranked first by mean validation R², with R²={w['R2_mean']:.3f} ± {w['R2_std']:.3f}, MAE={w['MAE_mean']:.2f} ± {w['MAE_std']:.2f} years and PCC={w['PCC_mean']:.3f} ± {w['PCC_std']:.3f} (mean ± sample SD across folds).")
        cs = table(f'TabulaStromal_Rebenchmark_2026/tables/{lin}/stability/final_clock_internal_cv_summary.csv')
        r = cs[cs.model_version == 'original_EN'].iloc[0]
        k = int(facts[lin]['summary']['K'])
        final.append(f'The {label} EN branch selected K={k:,} and achieved conditional validation R²={r.R2_mean:.3f} ± {r.R2_std:.3f} and MAE={r.MAE_mean:.2f} ± {r.MAE_std:.2f} years.')
        ps = table(f'TabulaStromal_Rebenchmark_2026/tables/{lin}/external/final_clock_external_protocol_comparison.csv')
        labels = ['No adaptation, all external','No adaptation, held-out 70%','After adaptation, held-out 70%','After adaptation, ALL100 (mixed)']
        for label_scope,(_,r) in zip(labels,ps.iterrows()):
            protocol_rows.append(f'| {lin} | {label_scope} | {int(r.n)} | {r.R2:.3f} | {r.MAE:.2f} |')
    external_table = '| Lineage | Evaluation population | n | R² | MAE (years) |\n| --- | --- | ---: | ---: | ---: |\n' + '\n'.join(protocol_rows)
    cohort = facts['cellage']
    for lin,ct in [('fibro','iFibro_CD74_pan'),('endo','EndoMT_PTN_pan')]:
        d = table(f'CellAge_Exploration/tables/{lin}_subtype_contrasts.csv')
        r = d[(d.score=='hUSI') & (d.celltype==ct)].iloc[0]
        examples.append(f'{ct} had a positive paired hUSI contrast against the equal-weight mean of other observed subtypes within the same sample/tissue (effect {r.effect:.4f}, 95% bootstrap CI {r.ci_low:.4f}–{r.ci_high:.4f}; BH q={r.q:.2g}; {int(r.n_donors)} source-designated sampling units).')
    section('Results_Part2_ClockPerformance_CellAge.md',
            'Results Part 2 — Celltype-resolved age clocks distinguish internal prediction, external transfer and cellular senescence', [
        ('Clock overview and complete grouped benchmark',
         '''TabulaStromal predicts recorded chronological age from sample-level Celltype_Gene expression features ([Fig. 2A](../../F2_ModelPerformance/figures/TabulaStromal_Workflow_V6.png)). Expression was aggregated within sample and annotated cell type; explicit cell-type proportions or ratios were excluded. Feature filtering, ranking and scaling were fitted within training partitions. The local DeepSets, Pseudobulk and Ensemble implementations are IDEAL-Age-inspired estimators and are not asserted to reproduce an external method protocol exactly.

The development cohorts contained 640 fibroblast samples in 634 conservative GSE/GSM proxy groups and 437 endothelial samples in 431 groups. Thirteen model families were evaluated at K=8,000, 10,000, 15,000 and 20,000 using grouped five-fold CV. Saved groups were disjoint within every training/validation split. This check cannot resolve unannotated donor reuse across studies.''' ,
         ['n_development','n_models','requested_K_values','n_outer_folds'],'Results2'),
        ('Internal performance and source-figure scope',
         '\n\n'.join(winners) + '''

This is a descriptive ranking across the completed benchmark, not a model-superiority significance test. The supplied [Fig. 2B–F](../../F2_ModelPerformance/figures/Fig2_StromalClock_ModelPerformance.pdf) uses source fold 4, with each model's K chosen from that fold's validation result. It is not a five-fold mean. Its external values use FT30 followed by ALL100 evaluation, including fitted adaptation samples. MAD is the mean absolute deviation of chronological-age labels from their mean; it describes label dispersion rather than prediction error. The figure's biological-age axis wording refers to prediction of the recorded chronological-age target.''',
         ['global_winner','source_fold4_'],'Results2'),
        ('A separate final interpretable ElasticNet branch',
         ' '.join(final) + ''' Both clocks retained the original selected feature representation because stable-feature restriction predicted less well. They were refitted to all development samples ([Fig. 3A,B](../../TabulaStromal_Rebenchmark_2026/figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.pdf)). Selection of the EN family and its K is separate from the global benchmark ranking. Because K was selected using the same internal comparison, these CV estimates are conditional on that selection rather than fully nested estimates of the complete selection procedure.''',
         ['conditional_EN_','K','internal_original','internal_stable'],'Results2'),
        ('External transfer depends on adaptation and the evaluation population',
         '''Without adaptation, the full-development fibroblast EN yielded R²=−2.943 and MAE=29.32 years in 93 merged external samples. The endothelial EN yielded R²=−1.720 and MAE=24.89 years in 77 samples (Fig. 3E,F). Negative R² preserves an important transfer limitation: the predictions had greater squared error than the external cohort-mean reference.

For the final EN adaptation protocol, 28 fibroblast and 23 endothelial samples were used for regularized refitting with frozen features, scaler and source hyperparameters, leaving 65 and 54 group-disjoint samples for evaluation. On the same held-out fibroblast samples, MAE decreased from 29.02 to 13.14 years and R² changed from −3.168 to −0.001. For endothelial samples, MAE decreased from 25.65 to 13.18 years and R² changed from −2.194 to −0.355. Adaptation reduced error, while held-out R² remained near zero or negative (Fig. 3M,N).

''' + external_table + '''

Pooling the fitted and held-out observations gives the more favorable ALL100 metrics shown in the last row for each lineage. Those pooled metrics are mixed apparent/unseen performance and are not independent validation. Primary no-adaptation performance, matched held-out adaptation comparisons and the historical source-fold figure therefore remain distinct results. No independent health-outcome analysis establishes these predictions as biological-age measurements.''',
         ['external_'],'Results2'),
        ('CellAge identifies subtype differences without a confirmed adjusted age association',
         f'''The audited CellAge snapshot included hUSI outputs for {cohort['fibro']['hUSI_cells']:,} fibroblast-lineage and {cohort['endo']['hUSI_cells']:,} endothelial-lineage cells; recSID recommendations covered {cohort['fibro']['recSID_cells']:,} and {cohort['endo']['recSID_cells']:,} cells, respectively ([Fig. 4](../../CellAge_Exploration/figures/Fig_CellAge_Stromal_Integration.pdf)). hUSI was analyzed as a senescence-associated score, not an age in years or a cross-lineage calibrated quantity. The inspected SenCID_rec_SID columns contain identity recommendations and do not provide genuine senescence-intensity measurements.

{' '.join(examples)} These differences describe subtype-associated scores; they are not slopes with chronological age. The original source calls these units donors under its sample=donor assumption, which the available cross-study metadata cannot independently verify.

At the primary threshold of at least 20 scored cells per sample/tissue/subtype, {cohort['fibro']['n_identifiable']} fibroblast and {cohort['endo']['n_identifiable']} endothelial adjusted hUSI–age models were identifiable; none met BH q<0.05 (Fig. 4E,F). The historical K=8,000 EN integration did not yield an identifiable fully adjusted hUSI–delta-age association. That partial historical prediction snapshot and its older stable features differ from the final EN clocks. A new joint analysis using complete final-clock OOF predictions and genuine SenCID intensities remains unavailable; it is not reported as a completed result.''',
         ['hUSI','recSID','paired_hUSI'],'Results2'),
    ])

    subtype_text = '; '.join(f'{k}, {v} '+('pair' if v==1 else 'pairs') for k,v in facts['Fibro']['subtype_counts'].items())
    shared = ', '.join(facts['shared_stable_genes'])
    section('Results_Part3_StableFeatures_Mechanisms.md',
            'Results Part 3 — Stable expression features nominate lineage-associated molecular hypotheses', [
        ('Group-bootstrap stability and subtype support',
         f'''The interpretable branch used 300 successful group bootstrap fits within each of five outer-training partitions and in the full development set, giving 1,800 checked records per lineage. Filtering, ranking and scaling were repeated within resamples. At selection frequency ≥0.70, full-development stability identified 438 fibroblast Celltype_Gene pairs representing 416 genes in five detailed cell types, and 23 endothelial pairs representing 23 genes in one cell type ([Fig. 3C,D,G–L](../../TabulaStromal_Rebenchmark_2026/figures/combined/Stromal_individual_aging_clock_Fibroblast_Endothelial_FINAL_v2.pdf)). Each held-out stable-model comparison used its corresponding outer-training stable set; full-development stability was used for exploratory explanation.

Fibroblast pair counts were {subtype_text}. The endothelial set was confined to CapEC_TM4SF1_pan. Candidate-normalized fractions and source cell/sample coverage accompany the raw counts. Unequal observation and candidate-gene opportunity limit interpretation of these counts as intrinsic subtype importance.''',
         ['n_bootstrap','successful_bootstraps','threshold','n_stable','n_unique','pan'],'Results3'),
        ('Stable-feature restriction did not improve age prediction',
         '''For fibroblasts, restricting EN to stable features reduced conditional mean validation R² from 0.742 to 0.636 and increased MAE from 11.86 to 14.65 years. For endothelial cells, R² decreased from 0.633 to 0.123 and MAE increased from 14.18 to 25.68 years (Fig. 3A,B). Neither stable version met the promotion criterion. The stable signatures therefore provide an interpretation of selected expression features, whereas the original-feature EN remains the final predictor. Selection stability does not itself imply superior prediction.''',
         ['internal_original','internal_stable'],'Results2'),
        ('Shared membership is not matched causal direction',
         f'''Nine genes were shared between the two final stable sets: {shared} ([shared-gene supplement](../../TabulaStromal_Rebenchmark_2026/figures/combined/Stromal_FinalClock_shared_specific_genes.png)). Shared membership means recurrence under the recorded selection procedures; it does not establish equal coefficient direction or independent causal replication. Apparent lineage-restricted membership can also reflect differences in sample support and feature observability.''',
         ['shared_stable_genes'],'Results3'),
        ('Functional annotation separates matrix-associated and vesicle-associated signals',
         '''Global fibroblast enrichment included Collagen-Containing Extracellular Matrix (GO:0062023; 27/416 query genes; 11,337-gene background; BH FDR=6.06×10⁻⁷) and the Reactome term Extracellular Matrix Organization (24/416 query genes; the same background; FDR=6.59×10⁻⁶). Endothelial enrichment included ER to Golgi Transport Vesicle Membrane (GO:0012507; 3/23 query genes; 10,949-gene background; FDR=0.0133; Fig. 3O,P). These are observed overlaps of stable genes with annotated gene sets, not demonstrations of pathway activation or causal aging mechanisms.

Across the specified global GO, Reactome and KEGG libraries, 56 fibroblast and 12 endothelial terms met FDR<0.05. Adding the separately retained Hallmark library increased the fibroblast count to 75. The exported global tables retain zero-overlap tested terms and the original within-scope/library BH values; the multiplicity family was not recomputed after filtering. Overlapping annotations driven by the same genes do not represent independent mechanistic discoveries.

The eight supplied aging lists remain provided_unverified. Their measured overlaps in Fig. 3Q,R are statistical comparisons with those supplied lists, not validation against authenticated official signatures.''',
         ['Matrix','Vesicle','Transport','COLLAGEN','Extracellular','n_significant_global','provided_unverified'],'Results3'),
        ('Expression–age associations define follow-up hypotheses',
         f'''The saved covariate-adjusted expression–age analyses contained FDR<0.05 results for {facts['Fibro']['n_age_associations']} of 438 fibroblast pairs and {facts['Endo']['n_age_associations']} of 23 endothelial pairs. These associations remain distinct from model-coefficient signs. Per-study direction resources can share samples and lack a prespecified universal replication criterion, so they do not independently validate every stable feature. The older CellAge signature branch is also distinct from these final stable sets. No perturbation experiment, authenticated GWAS result or transcription-factor network analysis establishes a molecular driver here. The supported interpretation is a set of predictive features and molecular hypotheses for further testing.''',
         ['FDR_significant_expression_age_pairs'],'Results3'),
    ])
    path = ROOT / 'audit/claim_to_evidence_report.csv'
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(claims[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(claims)
    print(json.dumps({'sections':3,'mapped_paragraphs':len(claims)},indent=2))


if __name__ == '__main__':
    main()
