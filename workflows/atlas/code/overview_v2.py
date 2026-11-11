"""User-reference Leiden UMAPs, mirrored coverage, compact markers and verified meta-cell trajectories.

Only plot existing coordinates/results. Never read atlas expression or refit a trajectory.
"""
import gc
import hashlib
import json
import logging
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, binary_closing, binary_dilation, binary_fill_holes

from utils_data import AGE_LABELS, age_order, dump_json, read_elem

HERE = Path(__file__).resolve().parent
TITLES = {'fibro': 'Fibroblast', 'endo': 'Endothelial'}


def load_options(overrides=None):
    options = json.loads((HERE / 'overview_replot_config.json').read_text())
    if overrides:
        options.update(overrides)
    if options['outline_group'] != 'subtype':
        raise ValueError('Overview A/B contours must follow Celltype_sub.')
    return options


def stamp(ax, label):
    ax.text(-.08, 1.08, label, transform=ax.transAxes, fontsize=16,
            fontweight='bold', va='bottom')


def save_figure(fig, stem, dpi):
    logging.info('Saving %s', stem)
    for suffix in ('.png', '.pdf'):
        fig.savefig(str(stem) + suffix, dpi=dpi, bbox_inches='tight', pad_inches=.14)
    plt.close(fig)


def verify_trajectory_inputs(config, root, options):
    """Fail on wrong paths, stale inputs, mismatched IDs/annotations/UMAP or wrong model cells."""
    from metacell_trajectory import load_model
    verified, data = {}, {}
    for lineage in ('fibro', 'endo'):
        trajectory = root / lineage / 'trajectory'
        audit = json.loads((trajectory / 'metacell_audit.json').read_text())
        expected = Path(options['expected_metacell_inputs'][lineage]).resolve()
        configured = Path(config['metacell_inputs'][lineage]).resolve()
        if configured != expected or Path(audit['input']).resolve() != expected:
            raise ValueError(f'{lineage}: trajectory input differs from the specified meta-cell object')
        stat = expected.stat()
        if stat.st_size != audit['source_size'] or stat.st_mtime_ns != audit['source_mtime_ns']:
            raise ValueError(f'{lineage}: meta-cell source changed after fitting')
        obs = pd.read_pickle(trajectory / 'prepared/metadata.pkl')
        cohort = pd.read_csv(trajectory / 'prepared/cohort_cells.csv', dtype={'cell_id': str})
        with h5py.File(expected, 'r') as handle:
            source_ids = np.asarray(read_elem(handle['obs'][handle['obs'].attrs['_index']])).astype(str)
            rows = cohort.source_row.to_numpy(dtype=int)
            ids_equal = np.array_equal(source_ids[rows], obs.index.astype(str)) and np.array_equal(cohort.cell_id, obs.index.astype(str))
            if not ids_equal:
                raise ValueError(f'{lineage}: original meta-cell identities/order do not match prepared fit')
            source_umap = handle['obsm/X_umap'][:].astype(np.float32)[rows, :2]
            cached_umap = np.load(trajectory / 'prepared/umap.npy')
            if not np.array_equal(source_umap, cached_umap):
                raise ValueError(f'{lineage}: source and cached meta-cell UMAP coordinates differ')
            general = np.asarray(read_elem(handle['obs/Celltype_general'])).astype(str)[rows]
            if not np.array_equal(general, obs.general.astype(str)):
                raise ValueError(f'{lineage}: meta-cell annotation/identity mismatch')
        choice = json.loads((trajectory / 'trajectory_selection.json').read_text())
        requested = options['trajectory_models'][lineage]
        if requested == 'selected':
            model_path, model_name = trajectory / choice['path'], choice['model']
        else:
            allowed = pd.read_csv(trajectory / 'three_method_comparison.csv').set_index('model')
            if requested not in allowed.index:
                raise ValueError(f'Unknown trajectory model {requested}')
            model_path, model_name = trajectory / allowed.loc[requested, 'path'], requested
        raw_ids = pd.read_csv(model_path / 'cells.csv', usecols=['cell_id'], dtype=str).cell_id
        if not raw_ids.is_unique or not set(raw_ids) <= set(obs.index):
            raise ValueError(f'{lineage}: native model cells are not unique original meta-cells')
        frame = load_model(model_path, obs)
        if len(frame) != len(obs):
            raise ValueError('Model cannot be aligned to the original meta-cell cohort')
        if model_path.name == 'cellrank':
            native = frame[['x', 'y']].to_numpy()
            finite = np.isfinite(native).all(axis=1)
            if not np.allclose(native[finite], source_umap[finite], atol=1e-5):
                raise ValueError('CellRank native coordinates differ from the original meta-cell UMAP')
        verified[lineage] = dict(input=str(expected), requested_input=options['expected_metacell_inputs'][lineage],
            source_size=stat.st_size, source_mtime_ns=stat.st_mtime_ns,
            source_meta_cells=len(source_ids), fitted_meta_cells=len(obs), native_model_meta_cells=len(raw_ids),
            model=model_name, model_path=str(model_path), source_ids_and_order_match=True,
            source_umap_exact_match=True, source_general_annotations_match=True,
            fitted_ids_sha256=hashlib.sha256('\n'.join(obs.index.astype(str)).encode()).hexdigest(),
            pseudotime_finite=int(np.isfinite(frame.pseudotime).sum()),
            native_embedding='UMAP' if model_path.name == 'cellrank' or model_path.name == 'historical_umap' else 'DDRTree' if model_path.parent.name == 'monocle' else 'STREAM PCA10 (displaying the first two components)')
        data[lineage] = dict(obs=obs, xy=source_umap, frame=frame, model_path=model_path,
                             model_name=model_name, source=expected)
        logging.info('%s trajectory verified: %s; %s/%s original meta-cells; model=%s',
                     lineage, expected, len(obs), len(source_ids), model_name)
    return verified, data


def reference_palette(obs, lineage, options):
    """Use the pasted base-to-white Celltype gradient and actual current labels."""
    bases = options['general_colors'][lineage]
    mapping = {}
    for general in obs.general.unique():
        if general not in bases:
            raise ValueError(f'No reference base color for {general}')
        names = obs.loc[obs.general == general, 'celltype'].unique()
        cmap = plt.get_cmap(LinearSegmentedColormap.from_list('', [bases[general], '#ffffff']), len(names))
        for i, name in enumerate(names):
            mapping[str(name)] = to_hex(cmap(i / (len(names) - 1))) if len(names) > 1 else bases[general]
    return bases, mapping


def subtype_contours(ax, obs, xy, bases, options):
    """99%-mass smoothed regions plus a small outward buffer; retain separate islands."""
    grid = options['outline_grid']
    span = np.maximum(np.ptp(xy, axis=0), 1e-6)
    low, high = np.min(xy, axis=0) - .045 * span, np.max(xy, axis=0) + .045 * span
    edges = [np.linspace(low[d], high[d], grid + 1) for d in (0, 1)]
    centers = [(e[:-1] + e[1:]) / 2 for e in edges]
    codes = obs.subtype.cat.codes.to_numpy() if isinstance(obs.subtype.dtype, pd.CategoricalDtype) else pd.Categorical(obs.subtype).codes
    categories = obs.subtype.cat.categories if isinstance(obs.subtype.dtype, pd.CategoricalDtype) else pd.Categorical(obs.subtype).categories
    paths = {}
    for i, subtype in enumerate(categories):
        mask = codes == i
        if mask.sum() < 15:
            continue
        general = obs.loc[mask, 'general'].mode().iloc[0]
        # Outline colors are darker than the last white legend swatch, for visibility.
        color = bases[general]
        histogram, _, _ = np.histogram2d(xy[mask, 0], xy[mask, 1], bins=edges)
        density = gaussian_filter(histogram.T, 1.05)
        ordered = np.sort(density.ravel())[::-1]
        level = ordered[min(np.searchsorted(np.cumsum(ordered), options['outline_mass'] * ordered.sum()), len(ordered) - 1)]
        region = density >= max(level, 1e-12)
        region = binary_fill_holes(binary_closing(region, iterations=1))
        region = binary_dilation(region, iterations=1)
        contour = ax.contour(centers[0], centers[1], region.astype(float), levels=[.5],
                             colors=[color], linewidths=1.1, linestyles='--', alpha=.85, zorder=5)
        segments = [s.tolist() for s in contour.allsegs[0] if len(s) >= 6]
        paths[str(subtype)] = dict(general=str(general), color=color, n_cells=int(mask.sum()), paths=segments)
    return paths


def reference_atlas(fig, spec, root, lineage, options, letter):
    import scanpy as sc
    from anndata import AnnData

    logging.info('Drawing %s: all atlas cells, Scanpy Leiden size=%s alpha=%s', lineage,
                 options['atlas_scatter_size'], options['atlas_scatter_alpha'])
    obs = pd.read_pickle(root / lineage / 'metadata.pkl')
    xy = np.load(root / lineage / 'umap.npy', mmap_mode='r')
    info = json.loads((root / 'shared' / f'data_audit_{lineage}.json').read_text())
    # Reserve space for the pasted two independent legends; preserve the plot's data aspect.
    nested = spec.subgridspec(1, 2, width_ratios=[1, .78], wspace=.015)
    ax = fig.add_subplot(nested[0, 0])
    adata = AnnData(obs=obs[['leiden', 'general', 'celltype']].rename(columns={
        'general': 'Celltype_general', 'celltype': 'Celltype'}).copy(), shape=(len(obs), 0))
    adata.obsm['X_umap'] = np.asarray(xy)
    with h5py.File(info['input'], 'r') as handle:
        if 'leiden_colors' in handle['uns']:
            colors = np.asarray(read_elem(handle['uns/leiden_colors'])).astype(str)
            if len(colors) != len(adata.obs.leiden.cat.categories):
                raise ValueError('Stored Leiden colors do not match categories')
            adata.uns['leiden_colors'] = colors
    # The pasted code: Leiden scatter, labels on data, s=10, alpha=.8.
    sc.pl.umap(adata, color='leiden', ax=ax, show=False, legend_loc='on data',
               title=f'{TITLES[lineage]} Atlas', s=options['atlas_scatter_size'],
               alpha=options['atlas_scatter_alpha'])
    for artist in ax.collections:
        artist.set_rasterized(True)
    ax.set_title(f'{TITLES[lineage]} Atlas', fontsize=16, fontweight='bold', pad=20)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set(xticks=[], yticks=[], xlabel='', ylabel='')
    x_range, y_range = ax.get_xlim(), ax.get_ylim()
    arrow_length_x = (x_range[1] - x_range[0]) * .1
    arrow_length_y = (y_range[1] - y_range[0]) * .1
    ax.annotate('UMAP_1', xy=(0, 0), xytext=(arrow_length_x, 0),
                arrowprops=dict(arrowstyle='<-', color='black', lw=1.5), fontsize=12, ha='left', va='center')
    ax.annotate('UMAP_2', xy=(0, 0), xytext=(0, arrow_length_y),
                arrowprops=dict(arrowstyle='<-', color='black', lw=1.5), fontsize=12,
                ha='center', va='bottom', rotation=90)
    bases, celltype_colors = reference_palette(obs, lineage, options)
    paths = subtype_contours(ax, obs, xy, bases, options)
    columns = [['qFibro', 'Myfibro', 'iFibro'], ['eFibro']] if lineage == 'fibro' else [['MCs', 'VenEC', 'LEC', 'EndoMT'], ['ArtEC', 'CapEC']]
    legends = []
    for column, groups in enumerate(columns):
        handles, labels = [], []
        for general in groups:
            if general not in set(obs.general):
                continue
            handles.append(mpatches.Patch(color=bases[general]))
            labels.append(general)
            for celltype in sorted(obs.loc[obs.general == general, 'celltype'].unique()):
                handles.append(mpatches.Patch(color=celltype_colors[celltype]))
                labels.append('  ' + celltype)
        legend = ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.015 + .40 * column, .96),
                           loc='upper left', borderaxespad=0, frameon=False, fontsize=10,
                           handlelength=.9, handletextpad=.35, labelspacing=.48)
        if column == 0:
            ax.add_artist(legend)
        legends.append(legend)
    # Retain broad annotation labels from the reference, but never replace subtype contours by ellipses.
    label_positions = []
    for general in obs.general.unique():
        points = xy[obs.general.to_numpy() == general]
        center = np.mean(points, axis=0)
        width = np.ptp(points[:, 0]) * 1.2
        proposed_x = center[0] + width / 2 + np.ptp(x_range) * .02
        right_limit = x_range[1] - .04 * np.ptp(x_range)
        label_x = min(proposed_x, right_limit)
        label_positions.append([label_x, float(center[1]), general])
    # Preserve the reference text boxes while separating colliding broad-state names.
    label_positions.sort(key=lambda entry: entry[1])
    separation = .065 * np.ptp(y_range)
    for i in range(1, len(label_positions)):
        label_positions[i][1] = max(label_positions[i][1], label_positions[i-1][1] + separation)
    shift = max(0, label_positions[-1][1] - (y_range[1] - .04 * np.ptp(y_range)))
    for label_x, label_y, general in label_positions:
        ax.text(label_x, label_y-shift, general, fontsize=12, fontweight='bold', color=bases[general],
                ha='right', va='center', zorder=7,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=bases[general], alpha=.7))
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.text(.0, -.05, f'{len(obs):,} cells; scatter and numbers = Leiden; dashed outlines = Celltype_sub',
            transform=ax.transAxes, fontsize=7)
    stamp(ax, letter)
    # Alignment of cached metadata/coordinates to the original source is checked explicitly.
    with h5py.File(info['input'], 'r') as handle:
        source_xy = handle['obsm/X_umap'][:].astype(np.float32)
        source_ids = np.asarray(read_elem(handle['obs'][handle['obs'].attrs['_index']])).astype(str)
        if not np.array_equal(source_xy, xy) or not np.array_equal(source_ids, obs.index.astype(str)):
            raise ValueError('Atlas cache is not aligned to the original full-cell input')
    provenance = dict(source=info['input'], n_source_cells=info['n_cells'], n_drawn_cells=len(obs),
                      sampled=False, scatter='scanpy.pl.umap(color=leiden, legend_loc=on data)',
                      size=options['atlas_scatter_size'], alpha=options['atlas_scatter_alpha'],
                      source_leiden_palette=list(adata.uns['leiden_colors']),
                      outline_group='Celltype_sub', outline_mass=options['outline_mass'],
                      n_subtype_outlines=len(paths), general_colors=bases,
                      celltype_legend_colors=celltype_colors, source_ids_and_umap_match=True)
    del adata
    gc.collect()
    return obs, paths, provenance


def common_coverage(root, metadata, top_tissues):
    counts = pd.concat([pd.read_csv(root / l / 'summary_tissue.csv', index_col=0).n_cells for l in metadata], axis=1).fillna(0).sum(axis=1)
    tissues = counts.drop(index='Unknown', errors='ignore').nlargest(top_tissues).index.tolist()
    tissue_tables, age_tables = {}, {}
    for lineage, obs in metadata.items():
        group = obs.tissue.astype(str).where(obs.tissue.isin(tissues), 'Others')
        frame = pd.DataFrame({'group': group.to_numpy(), 'sample': obs['sample'].to_numpy()})
        cells = frame.groupby('group', observed=True).size().rename('n_cells')
        samples = frame.loc[frame['sample'] != 'Unknown'].groupby('group', observed=True)['sample'].nunique().rename('n_samples')
        tissue_tables[lineage] = pd.concat([cells, samples], axis=1).fillna(0).reindex(tissues + ['Others'], fill_value=0).astype(int)
        age_tables[lineage] = pd.read_csv(root / lineage / 'summary_age_group.csv', index_col=0)
    ages = age_order(pd.Index(age_tables['fibro'].index).union(age_tables['endo'].index))
    age_tables = {l: t.reindex(ages, fill_value=0) for l, t in age_tables.items()}
    return tissue_tables, age_tables


def format_count(value):
    return f'{value / 1e6:.2f}M' if value >= 1e6 else f'{value / 1e3:.1f}k' if value >= 10000 else f'{value:,}'


def mirrored_chart(fig, spec, tables, metric, title):
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    grid = spec.subgridspec(1, 3, width_ratios=[1, .50, 1], wspace=.015)
    fib, middle, endo = [fig.add_subplot(grid[0, i]) for i in range(3)]
    maximum = max(t[metric].max() for t in tables.values()) * 1.22
    maximum = max(maximum, 1)
    labels = list(tables['fibro'].index)
    y = np.arange(len(labels))
    for ax, lineage in ((fib, 'fibro'), (endo, 'endo')):
        vals = tables[lineage][metric].to_numpy()
        color = '#649FB9' if lineage == 'fibro' else '#D39B70'
        ax.barh(y, vals, color=color, height=.66)
        # Fibro grows from the left towards the middle; Endo grows from right to left.
        ax.set_xlim((0, maximum) if lineage == 'fibro' else (maximum, 0))
        ax.set_ylim(len(labels) - .5, -.6)
        ax.set_yticks([])
        ax.spines[['left', 'right', 'top']].set_visible(False)
        ax.xaxis.set_major_locator(MaxNLocator(3, integer=True))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: format_count(int(x))))
        ax.tick_params(axis='x', labelsize=6)
        ax.set_title(TITLES[lineage], color=color, fontsize=8, pad=5)
        for i, value in enumerate(vals):
            ax.text(value + .015 * maximum, i, format_count(value),
                    ha='left' if lineage == 'fibro' else 'right', va='center', fontsize=6)
    middle.set_xlim(0, 1)
    middle.set_ylim(len(labels) - .5, -.6)
    middle.set_axis_off()
    for i, label in enumerate(labels):
        middle.text(.5, i, AGE_LABELS.get(label, label), ha='center', va='center', fontsize=7)
    middle.set_title(title, fontsize=9, fontweight='bold', pad=24)
    return fib


def coverage_panel(fig, spec, tissue, age):
    grid = spec.subgridspec(2, 2, height_ratios=[1.35, 1], wspace=.20, hspace=.68)
    first = None
    for row, tables, name in ((0, tissue, 'Tissue'), (1, age, 'Age group')):
        for col, metric, unit in ((0, 'n_cells', 'cells'), (1, 'n_samples', 'samples')):
            ax = mirrored_chart(fig, grid[row, col], tables, metric, f'{name} · {unit}')
            if first is None:
                first = ax
    stamp(first, 'C')


def compact_marker(ax, root, lineage, fig):
    from utils_plot import marker_panel
    marker_panel(ax, pd.read_csv(root / lineage / 'marker_summary.csv'), lineage, fig)
    ax.set_box_aspect(.56)
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=7)
    ax.title.set_fontsize(10)
    ax.get_legend().set_bbox_to_anchor((.5, -.38))


def compact_signature(ax, root, lineage, fig):
    from utils_plot import signature_panel
    signature_panel(ax, pd.read_csv(root / lineage / 'signature_mean.csv', index_col=0), lineage, fig)
    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)
    ax.title.set_fontsize(10)


def trajectory_pair(fig, spec, lineage, item, letter):
    from metacell_trajectory import draw_model
    inner = spec.subgridspec(1, 2, width_ratios=[1, 1], wspace=.25)
    ax, projected = [fig.add_subplot(inner[0, i]) for i in range(2)]
    draw_model(ax, item['model_path'], item['obs'], lineage)
    # draw_model's standalone footer would collide with the embedding caption here.
    for text in list(ax.texts):
        if 'meta-cells; star = assumed root' in text.get_text():
            text.remove()
    frame, xy = item['frame'], item['xy']
    good = np.isfinite(frame.pseudotime.to_numpy())
    projected.scatter(xy[:, 0], xy[:, 1], c='#dddddd', s=1, lw=0, rasterized=True)
    dots = projected.scatter(xy[good, 0], xy[good, 1], c=frame.pseudotime.to_numpy()[good],
                             cmap='viridis', vmin=0, vmax=1, s=1.7, lw=0, rasterized=True)
    root = np.nanargmin(frame.pseudotime.to_numpy())
    projected.scatter(*xy[root], marker='*', s=75, c='white', edgecolor='#222222', zorder=10)
    projected.set_aspect('equal', adjustable='datalim')
    projected.set(xticks=[], yticks=[])
    for spine in projected.spines.values():
        spine.set_visible(False)
    projected.set_title('Same pseudotime on input meta-cell UMAP', fontsize=9, loc='left')
    fig.colorbar(dots, ax=projected, fraction=.025, pad=.01, label='Root-relative pseudotime')
    embedding = 'UMAP' if item['model_path'].name in ['cellrank', 'historical_umap'] else 'DDRTree components' if item['model_path'].parent.name == 'monocle' else 'PCA components 1 / 2'
    ax.set_xlabel(f'{embedding}; {len(item["obs"]):,} meta-cells; star = assumed root', fontsize=6.5, labelpad=8)
    projected.set_xlabel('Original input UMAP (no recomputation)', fontsize=7)
    projected.text(.0, -.12, item['source'].name + f'\n{len(xy):,} original meta-cells; median constituent cells ≈30',
                   transform=projected.transAxes, fontsize=6.5)
    stamp(ax, letter)


def render_overview(config, root, output_dir=None, options=None):
    root = Path(root)
    options = load_options(options)
    out = Path(output_dir) if output_dir else root / 'integrated_figure'
    out.mkdir(parents=True, exist_ok=True)
    verification, trajectories = verify_trajectory_inputs(config, root, options)
    dump_json(out / 'trajectory_input_verification.json', verification)
    figure = plt.figure(figsize=options['figure_size'])
    outer = figure.add_gridspec(4, 2, height_ratios=options.get('row_height_ratios', [6.7, 5.6, 3.5, 4.2]), hspace=.35, wspace=.16,
                               left=.045, right=.97, bottom=.065, top=.96)
    metadata, atlas_provenance = {}, {}
    for col, lineage in enumerate(('fibro', 'endo')):
        obs, paths, provenance = reference_atlas(figure, outer[0, col], root, lineage, options, 'AB'[col])
        metadata[lineage] = obs
        atlas_provenance[lineage] = provenance
        dump_json(out / f'{lineage}_subtype_contours.json', dict(method='Smoothed 99%-mass density regions with outward grid buffer; separate islands retained', groups=paths))
    tissue, age = common_coverage(root, metadata, options['top_tissues'])
    for kind, tables in [('tissue', tissue), ('age_group', age)]:
        result = pd.concat([t.assign(lineage=l) for l, t in tables.items()]).rename_axis('group').reset_index()
        result.to_csv(out / f'C_{kind}_mirrored_counts.csv', index=False)
    coverage_panel(figure, outer[1, :], tissue, age)
    del metadata
    gc.collect()
    compact = outer[2, :].subgridspec(1, 4, width_ratios=[1.22, 1.22, .86, 1.02], wspace=.65)
    for col, lineage in enumerate(('fibro', 'endo')):
        ax = figure.add_subplot(compact[0, col])
        compact_marker(ax, root, lineage, figure)
        stamp(ax, 'EF'[col])
        ax = figure.add_subplot(compact[0, 2 + col])
        compact_signature(ax, root, lineage, figure)
        stamp(ax, 'GH'[col])
        trajectory_pair(figure, outer[3, col], lineage, trajectories[lineage], 'IJ'[col])
    footer = ('A/B: original atlas Leiden UMAPs, all cells; subtype contours replace ellipses. C: inward mirrored full-metadata coverage; shared positive scales.\n'
              'E/F: marker detection fraction and gene-scaled expression. G/H: raw signature means; colors = column z-scores. I/J: native fits + same pseudotime on original meta-cell UMAP.\n'
              'Age A/H retain historical mixed labels; samples are not donors. Roots and direction are assumptions; meta-cell sizes vary; CellRank macrostate discovery uses 200 landmarks.')
    figure.text(.045, .012, footer, fontsize=7.5, linespacing=1.5)
    save_figure(figure, out / 'Fig_StromalAtlas_Overview', options['dpi'])
    caption = ('A–B: user-reference Scanpy Leiden scatter (s=10, alpha=.8), Leiden labels on data, two independent general/Celltype annotation legends, arrow axes; '
               'all source cells, original coordinates and Leiden palette. Only contour grouping changes to Celltype_sub; no rectangles/ellipses. '
               'Current Endo MCs/LEC/EndoMT labels are retained; legacy PC/SMC, LymEC and EndoMT-I/II names are not fabricated.\n'
               'C: merged old C/D; common tissue categories and ages, Fibro extends from left to right, Endo from right to left. Counts use positive shared scales; '
               'Others unique samples are recomputed after grouping, not summed. D is intentionally merged into C; E–J retain previous panel letters.\n'
               'E–F: compact marker dotplots, all-cell statistics unchanged. G–H: compact square signature matrices, all-cell statistics unchanged.\n'
               'I–J: selected native trajectory plus its root-relative pseudotime on the unchanged original meta-cell UMAP; input, IDs, annotations and coordinates verified.\n'
               + json.dumps(verification, ensure_ascii=False, indent=2) + '\n' + footer + '\n')
    (out / 'figure_caption.txt').write_text(caption)
    dump_json(out / 'overview_plot_provenance.json', dict(layout_version=2, run_dir=str(root.resolve()),
              atlas=atlas_provenance, trajectory=verification, options=options,
              panel_letters=['A', 'B', 'C', 'E', 'F', 'G', 'H', 'I', 'J'],
              panel_D='Merged into C', analysis_refitted=False))
    # Export individual A/B at the reference 12x6 size, and compact E–H / paired I/J.
    for lineage, letter in [('fibro', 'A'), ('endo', 'B')]:
        fig = plt.figure(figsize=(12, 6))
        reference_atlas(fig, fig.add_gridspec(1, 1, left=.06, right=.99, top=.90, bottom=.10)[0], root, lineage, options, letter)
        save_figure(fig, out / f'Panel_{letter}_{lineage}_Leiden_subtype', options['dpi'])
    fig = plt.figure(figsize=(16, 8))
    coverage_panel(fig, fig.add_gridspec(1, 1, left=.06, right=.97, top=.90, bottom=.06)[0], tissue, age)
    save_figure(fig, out / 'Panel_C_mirrored_coverage', options['dpi'])
    for lineage, marker_letter, signature_letter, trajectory_letter in [('fibro', 'E', 'G', 'I'), ('endo', 'F', 'H', 'J')]:
        fig, ax = plt.subplots(figsize=(7, 4.8))
        compact_marker(ax, root, lineage, fig)
        stamp(ax, marker_letter)
        save_figure(fig, out / f'Panel_{marker_letter}_{lineage}_markers', options['dpi'])
        fig, ax = plt.subplots(figsize=(5.2, 5.0))
        compact_signature(ax, root, lineage, fig)
        stamp(ax, signature_letter)
        save_figure(fig, out / f'Panel_{signature_letter}_{lineage}_signatures', options['dpi'])
        fig = plt.figure(figsize=(12, 5.5))
        trajectory_pair(fig, fig.add_gridspec(1, 1, left=.05, right=.95, top=.88, bottom=.15)[0], lineage, trajectories[lineage], trajectory_letter)
        save_figure(fig, out / f'Panel_{trajectory_letter}_{lineage}_verified_metacells', options['dpi'])
    return out
