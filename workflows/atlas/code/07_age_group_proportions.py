#!/usr/bin/env python
"""Sample-level age-group composition boxplots from the two requested age-filtered atlases.

Read obs only; never load X. Main figures: one general state per row, general proportion
at left, up to five Celltype proportions at right. All percentages use the sample's
total cells in the same lineage object as denominator, including zero observations.
"""
import argparse
import hashlib
import json
import logging
from pathlib import Path
import shutil
import zipfile

import stromal_atlas_pipeline as pipeline  # Project-local caches and headless plotting.
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr

from utils_data import read_elem

BASE = pipeline.HERE.parent.parent
INPUTS = {
    'fibro': BASE / 'data/scalex/fibro/human/Fibro_human_FibDEG_count_notAge_0.2_10X_tissueB_bbknn_onlyFib_age.h5ad',
    'endo': BASE / 'data/scalex/endo/human/Endo_human_EndoDEG_count_notAge_0.2_10X_tissueB_bbknn_onlyEndo_age.h5ad',
}
ORDER = {'fibro': ['iFibro', 'eFibro', 'Myfibro', 'qFibro'],
         'endo': ['ArtEC', 'CapEC', 'VenEC', 'EndoMT', 'MCs', 'LEC']}
NAMES = {'fibro': 'Fibroblast', 'endo': 'Endothelial'}
AGES = list('ABCDEFGH')
AGE_COLORS = ['#71B4A3', '#EFA17D', '#95A1BE', '#D48AB5', '#AACA61', '#E2CE57', '#D9BC92', '#8CC4D0']


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    result = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    idx = valid[np.argsort(p[valid], kind='stable')]
    adjusted = p[idx] * len(idx) / np.arange(1, len(idx) + 1)
    if len(idx):
        result[idx] = np.minimum(1., np.minimum.accumulate(adjusted[::-1])[::-1])
    return result


def read_metadata(path):
    fields = ['Celltype_general', 'Celltype', 'Age_group', 'sample', 'GSE', 'tissue']
    with h5py.File(path, 'r') as handle:
        missing = set(fields) - set(handle['obs'])
        if missing:
            raise ValueError(f'{path}: missing fields {missing}')
        obs = pd.DataFrame({field: read_elem(handle['obs'][field]) for field in fields})
        source_cells = len(handle['obs'][handle['obs'].attrs['_index']])
    if len(obs) != source_cells or obs.isna().any().any():
        raise ValueError('Missing or inconsistent obs metadata; no cells silently removed')
    for field in fields:
        values = obs[field].astype(str).str.strip()
        if values.str.lower().isin(['', 'nan', 'unknown', 'none', 'na']).any():
            raise ValueError(f'Unknown {field}; inspect before defining sample proportions')
        obs[field] = pd.Categorical(values)
    if set(obs.Age_group) - set(AGES):
        raise ValueError('Unexpected Age_group labels; do not guess chronological order')
    conflicts = obs.groupby('sample', observed=True)[['Age_group', 'GSE', 'tissue']].nunique()
    if (conflicts > 1).any().any():
        raise ValueError('A sample has conflicting age/study/tissue annotations; resolve before plotting')
    mapping = obs[['Celltype_general', 'Celltype']].drop_duplicates()
    if mapping.groupby('Celltype', observed=True).Celltype_general.nunique().max() != 1:
        raise ValueError('Celltype does not map uniquely to Celltype_general')
    return obs, mapping


def make_proportions(obs, field, samples):
    counts = pd.crosstab(obs['sample'], obs[field], dropna=True).reindex(samples.index, fill_value=0)
    if not (counts.sum(axis=1).to_numpy() == samples.n_cells.to_numpy()).all():
        raise ValueError('Cell counts lost when building sample-by-label table')
    proportions = counts.div(samples.n_cells, axis=0) * 100.
    if not np.allclose(proportions.sum(axis=1), 100., atol=1e-9):
        raise ValueError('Within-sample proportions do not sum to 100%')
    return counts, proportions


def statistics(proportions, samples, level):
    age_rank = samples.Age_group.astype(str).map({a: i for i, a in enumerate(AGES)}).to_numpy()
    groups = samples.Age_group.astype(str).to_numpy()
    records = []
    for label in proportions:
        values = proportions[label].to_numpy()
        if np.ptp(values) == 0:
            rho, p_trend, h_stat, p_kw = 0., 1., 0., 1.
        else:
            rho, p_trend = spearmanr(age_rank, values)
            arrays = [values[groups == age] for age in AGES if np.any(groups == age)]
            h_stat, p_kw = kruskal(*arrays) if len(arrays) > 1 else (np.nan, np.nan)
        records.append(dict(level=level, label=label, n_samples=len(values),
            n_nonzero_samples=int((values > 0).sum()), n_nonzero_age_groups=int(pd.Series(groups[values > 0]).nunique()),
            spearman_rho=float(rho), trend_p=float(p_trend), kruskal_H=float(h_stat), kruskal_p=float(p_kw)))
    return records


def ptext(value):
    if not np.isfinite(value):
        return 'NA'
    return f'{value:.1e}' if value < .001 else f'{value:.3f}'


def stars(qvalue):
    return '***' if qvalue < .001 else '**' if qvalue < .01 else '*' if qvalue < .05 else 'ns'


def panel(ax, label, values, samples, stats, general=False):
    groups = samples.Age_group.astype(str).to_numpy()
    arrays = [values[groups == age] for age in AGES]
    boxes = ax.boxplot(arrays, positions=np.arange(1, 9), widths=.65, patch_artist=True,
        showfliers=False, whis=1.5,
        boxprops=dict(edgecolor='#737373', linewidth=.75),
        whiskerprops=dict(color='#858585', linewidth=.75), capprops=dict(color='#858585', linewidth=.75),
        medianprops=dict(color='#454545', linewidth=1.))
    for box, color in zip(boxes['boxes'], AGE_COLORS):
        box.set_facecolor(color)
        box.set_alpha(.90)
    # Zero-filled tissue-specific types often have zero IQR in every age bin.
    # Display the actual samples and their full range instead of making nonzero
    # observations disappear as suppressed boxplot outliers.
    rng = np.random.default_rng(42)
    for i, array in enumerate(arrays, 1):
        ax.scatter(i + rng.uniform(-.18, .18, len(array)), array, s=3.0,
                   color='#485360', alpha=.24, linewidths=0, rasterized=True, zorder=3)
    top = float(np.max(values)) if len(values) else 0.
    ax.set_ylim(0, 105 if general else max(top * 1.10, .01))
    ax.set_xticks(range(1, 9), AGES, fontsize=7)
    ax.set_xlabel('Age group', fontsize=7, labelpad=2)
    ax.set_ylabel('Cell proportion (%)', fontsize=7)
    ax.tick_params(axis='y', labelsize=7)
    ax.grid(axis='y', color='#dddddd', linewidth=.6)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color('#c2c2c2')
        spine.set_linewidth(.7)
    ax.set_title(f'{label} ({stars(stats.kruskal_q)})', fontsize=8.4, fontweight='bold', pad=24)
    ax.text(.5, 1.035, f'ρ = {stats.spearman_rho:+.3f}; q(trend) = {ptext(stats.trend_q)}\n'
            f'q(KW) = {ptext(stats.kruskal_q)}', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=6.4, linespacing=1.25)


def render(lineage, samples, general_prop, celltype_prop, stats, selected, out, dpi):
    states = [s for s in ORDER[lineage] if s in general_prop]
    fig, axes = plt.subplots(len(states), 6, figsize=(22.8, 3.05 * len(states) + 1.45), squeeze=False)
    fig.subplots_adjust(left=.075, right=.99, bottom=.080 if len(states) == 4 else .065,
                        top=.89 if len(states) == 4 else .92, wspace=.50, hspace=.96)
    for row, general in enumerate(states):
        ax = axes[row, 0]
        record = stats.loc[(stats.level == 'Celltype_general') & (stats.label == general)].iloc[0]
        panel(ax, general, general_prop[general].to_numpy(), samples, record, general=True)
        ax.text(-.38, .5, general, transform=ax.transAxes, rotation=90, ha='center', va='center',
                fontsize=11, fontweight='bold', bbox=dict(facecolor='#fff16b', edgecolor='#d6cb5e', boxstyle='square,pad=.6'))
        names = selected.loc[(selected.Celltype_general == general) & selected.selected_for_plot, 'label'].tolist()
        for col in range(1, 6):
            ax = axes[row, col]
            if col > len(names):
                ax.set_axis_off()
                if col == len(names) + 1:
                    ax.text(.5, .5, f'{len(names)} observed Celltype'+('' if len(names) == 1 else 's')+'\n(all shown)',
                            ha='center', va='center', transform=ax.transAxes, fontsize=9, color='#888888')
                continue
            label = names[col - 1]
            record = stats.loc[(stats.level == 'Celltype') & (stats.label == label)].iloc[0]
            panel(ax, label, celltype_prop[label].to_numpy(), samples, record)
    fig.suptitle(f'{NAMES[lineage]} | cell-type proportions across age groups', fontsize=16, fontweight='bold', y=.985)
    fig.text(.075, .955 if len(states) == 4 else .958, 'Celltype_general', fontsize=12, fontweight='bold')
    fig.text(.36, .955 if len(states) == 4 else .958,
             'Celltype within each general state · top 5 by |Spearman age-trend ρ| when >5 types', fontsize=12, fontweight='bold')
    counts = samples.Age_group.astype(str).value_counts()
    coverage = '; '.join(f'{a}: n={counts.get(a, 0):,}' for a in AGES)
    footer = (f'{len(samples):,} samples; {int(samples.n_cells.sum()):,} cells. Samples per age group: {coverage}.\n'
              'Denominator = all cells in the same lineage object and sample; missing cell types are 0%. Boxes = IQR; line = median; whiskers = 1.5×IQR; dots = all sample values.\n'
              'ρ: sample-level Spearman correlation with ordered A–H bins. Stars: BH-adjusted Kruskal–Wallis q. Both tests adjusted over all general and fine labels within this lineage, before selection.\n'
              'Exploratory pooled associations; not adjusted for tissue/study, repeated donors or other confounders. A–H retain the input Age_group labels.')
    fig.text(.075, .012, footer, fontsize=7.3, linespacing=1.5)
    stem = out / f'{NAMES[lineage]}_AgeGroup_Proportion_Boxplots'
    for suffix in ('.png', '.pdf'):
        fig.savefig(str(stem) + suffix, dpi=dpi, bbox_inches='tight', pad_inches=.18)
    plt.close(fig)
    return stem


def analyse(lineage, path, out, dpi):
    logging.info('Reading %s obs only: %s', lineage, path)
    before = path.stat()
    obs, mapping = read_metadata(path)
    if set(obs.Celltype_general) - set(ORDER[lineage]):
        raise ValueError('Unexpected general states')
    samples = obs.groupby('sample', observed=True).agg(Age_group=('Age_group', 'first'),
        GSE=('GSE', 'first'), tissue=('tissue', 'first'), n_cells=('sample', 'size'))
    gcounts, gprop = make_proportions(obs, 'Celltype_general', samples)
    ccounts, cprop = make_proportions(obs, 'Celltype', samples)
    metrics = pd.DataFrame(statistics(gprop, samples, 'Celltype_general') + statistics(cprop, samples, 'Celltype'))
    # Correct over ALL tested labels, before selecting top-five displayed Celltypes.
    metrics['trend_q'] = bh_adjust(metrics.trend_p)
    metrics['kruskal_q'] = bh_adjust(metrics.kruskal_p)
    selected = metrics.loc[metrics.level == 'Celltype'].merge(mapping.astype(str), left_on='label', right_on='Celltype', validate='one_to_one')
    selected['abs_spearman_rho'] = selected.spearman_rho.abs()
    selected = selected.sort_values(['Celltype_general', 'abs_spearman_rho', 'trend_q', 'n_nonzero_samples', 'label'],
                                    ascending=[True, False, True, False, True], kind='stable')
    selected['rank_within_general'] = selected.groupby('Celltype_general', sort=False).cumcount() + 1
    selected['n_celltypes_in_general'] = selected.groupby('Celltype_general', sort=False).label.transform('size')
    selected['selected_for_plot'] = selected.rank_within_general <= 5
    for general in gprop:
        cols = mapping.loc[mapping.Celltype_general == general, 'Celltype'].astype(str)
        if not np.allclose(cprop[cols].sum(axis=1), gprop[general], atol=1e-9):
            raise ValueError('Child proportions fail to sum to their general parent')
    samples.to_csv(out / f'{lineage}_sample_metadata.csv')
    metrics.to_csv(out / f'{lineage}_all_statistics.csv', index=False)
    selected.to_csv(out / f'{lineage}_celltype_selection.csv', index=False)
    long_tables = []
    age_summary = []
    for level, counts, proportions in [('Celltype_general', gcounts, gprop), ('Celltype', ccounts, cprop)]:
        proportions.to_csv(out / f'{lineage}_{level}_sample_percent.csv')
        for label in proportions:
            table = samples.reset_index().copy()
            table['level'] = level
            table['label'] = label
            table['n_label_cells'] = counts[label].to_numpy()
            table['percent'] = proportions[label].to_numpy()
            long_tables.append(table)
            summary = table.groupby('Age_group', observed=True).agg(n_samples=('sample', 'size'),
                mean_percent=('percent', 'mean'), median_percent=('percent', 'median'),
                q25_percent=('percent', lambda x: x.quantile(.25)), q75_percent=('percent', lambda x: x.quantile(.75)),
                n_label_cells=('n_label_cells', 'sum'), denominator_cells=('n_cells', 'sum')).reset_index()
            summary['pooled_cell_percent'] = 100 * summary.n_label_cells / summary.denominator_cells
            summary['level'], summary['label'] = level, label
            age_summary.append(summary)
    pd.concat(long_tables, ignore_index=True).to_csv(out / f'{lineage}_sample_proportions_long.csv.gz', index=False)
    pd.concat(age_summary, ignore_index=True).to_csv(out / f'{lineage}_age_group_summary.csv', index=False)
    stem = render(lineage, samples, gprop, cprop, metrics, selected, out, dpi)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Source file changed while reading metadata')
    result = dict(input=str(path), n_cells=len(obs), n_samples=len(samples), n_studies=int(samples.GSE.nunique()),
        n_tissues=int(samples.tissue.nunique()), cells_by_age=obs.Age_group.value_counts().to_dict(),
        samples_by_age=samples.Age_group.value_counts().to_dict(),
        source_bytes=before.st_size, source_mtime_ns=before.st_mtime_ns, source_unchanged=True,
        metadata_only=True, sample_age_study_tissue_conflicts=0, denominator='all cells of this lineage within the same sample',
        zero_filled=True, sample_percent_sums_100=True, child_percent_sums_parent=True,
        selection='top min(5, n_Celltype) within Celltype_general by absolute sample-level Spearman rho; ties: trend q, nonzero samples, label',
        p_adjustment='Benjamini–Hochberg separately for trend and KW, all general+fine labels within lineage before selection',
        files=[str(stem)+'.png', str(stem)+'.pdf'])
    logging.info('%s complete: %s cells; %s samples', lineage, len(obs), len(samples))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=pipeline.HERE.parent / 'results/age_group_proportions_20261102')
    parser.add_argument('--dpi', type=int, default=300)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
        handlers=[logging.StreamHandler(), logging.FileHandler(out / 'run.log')], force=True)
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'savefig.facecolor': 'white'})
    results = {lin: analyse(lin, path, out, args.dpi) for lin, path in INPUTS.items()}
    (out / 'analysis_audit.json').write_text(json.dumps(results, indent=2, default=lambda x: x.item()))
    source = Path(__file__).resolve()
    shutil.copy2(source, out / source.name)
    (out / 'code_sha256.txt').write_text(hashlib.sha256(source.read_bytes()).hexdigest() + '\n')
    (out / 'README_说明.md').write_text('''# 年龄分组细胞比例图

每个大图一行对应一个 Celltype_general。左列为该大类占样本中本谱系全部细胞的百分比；右侧为该大类的 Celltype，分母相同。这里“subtype”按参考图解释为每行的大类 Celltype_general，不按 Celltype_sub 再拆出数十行。

每个 sample 是一个箱线图观测值；比例按该 sample 在对应 h5ad 中的总细胞数计算。未检出的类别补零，不能仅保留阳性样本。此图不包含 h5ad 以外的其他细胞类型，分母不是整份组织中所有细胞。样本不是已经确认独立的供体。

每个大类的 Celltype 超过 5 个时，按样本比例与有序 Age_group A–H 的 Spearman |rho| 选前 5 个；不超过 5 个则全部展示。所选趋势可以是上升或下降；排序不保证显著。完整排序及全部未选类别见 *_celltype_selection.csv。

图中 rho 和 q(trend) 对应样本层面 Spearman 趋势；q(KW) 及星号对应跨年龄组 Kruskal–Wallis 差异。星号阈值：ns >=0.05，* <0.05，** <0.01，*** <0.001。每个谱系中，两个检验各自对所有 general 和 Celltype 检验做 BH 校正后再筛图，不对挑出的 5 个单独校正。检验都是未经组织、研究、重复供体或其他混杂调整的探索性汇总关联；筛选后图形不能当作独立验证，不能据此认定衰老因果效应。

箱体为第25–75百分位，线为中位数，须为1.5倍IQR，叠加散点为全部样本值（包含箱须以外的点），纵轴保留全部样本范围。部分组织特异的 Celltype 在超过75%的样本中未检出，箱体会压在0，此时散点展示少数阳性样本，不能删去阴性样本来扩大箱体。各面板单独确定纵轴范围，大类比例统一0–100%。A–H保留原始 Age_group 标签及顺序，没有擅自重分年龄。每个样本均计入，无人为最小细胞数筛选；样本细胞数见 *_sample_metadata.csv，可据此开展后续敏感性检查。

*_sample_proportions_long.csv.gz 包含每个样本、每个类别的原始计数、分母和百分比；*_age_group_summary.csv 包含样本均值、中位数、IQR及按细胞汇总的比例（汇总比例不用于画样本箱线图）；*_all_statistics.csv 保存全部统计检验；analysis_audit.json 保存实际输入、样本覆盖及输入未改变的核验。

程序只读取 obs，不加载表达矩阵，不写回原 h5ad。复现：在 StromalAtlas/code 下运行 `python 07_age_group_proportions.py`，默认结果写入本目录。
''')
    archive = out / 'StromalAtlas_AgeGroup_Proportions_Figures.zip'
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as bundle:
        for name in ['Fibroblast_AgeGroup_Proportion_Boxplots.png', 'Fibroblast_AgeGroup_Proportion_Boxplots.pdf',
                     'Endothelial_AgeGroup_Proportion_Boxplots.png', 'Endothelial_AgeGroup_Proportion_Boxplots.pdf',
                     'fibro_celltype_selection.csv', 'endo_celltype_selection.csv', 'README_说明.md', 'analysis_audit.json']:
            bundle.write(out / name, name)
    print('RESULT_DIR=' + str(out))
    print('DOWNLOAD_ZIP=' + str(archive))


if __name__ == '__main__':
    main()
