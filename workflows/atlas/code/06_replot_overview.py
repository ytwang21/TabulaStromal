#!/usr/bin/env python
"""Replot the completed full-meta-cell run without refitting any analysis."""
import argparse
import datetime as dt
import hashlib
import json
import logging
from pathlib import Path
import shutil

# Set project-local plotting/Numba caches before importing Scanpy or Matplotlib.
import stromal_atlas_pipeline as pipeline


def write_notebook(run, revision):
    import nbformat as nbf
    md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
    notebook = nbf.v4.new_notebook(metadata={'kernelspec': {
        'display_name': 'Python 3', 'language': 'python', 'name': 'python3'}})
    notebook.cells = [
        md('# StromalAtlas：按用户参考代码重新排版\n\n默认只查看已生成的结果，不重跑分析。A/B 是原始全细胞 Leiden UMAP；轮廓按 Celltype_sub。C 合并原 C/D，左右相向柱图；E–H 紧凑显示；I/J 同时显示原生轨迹和原始 meta-cell UMAP 拟时序。'),
        code('from pathlib import Path\nimport json, subprocess, sys\nimport pandas as pd\nfrom IPython.display import display, Image\n'
             f'CODE = Path({str(pipeline.HERE)!r})\nRUN_DIR = Path({str(run)!r})\nRESULT = Path({str(revision)!r})\n'
             'def show(name, width=1400):\n    display(Image(filename=str(RESULT/name), width=width))\n'),
        md('如需重新产图，手动运行下面单元；不会执行 Monocle、STREAM 或 CellRank 的拟合。'),
        code('RENDER_AGAIN = False\nif RENDER_AGAIN:\n    subprocess.run([sys.executable, str(CODE/"06_replot_overview.py"), "--run-dir", str(RUN_DIR)], check=True, cwd=CODE)\n    RESULT = Path(json.loads((RUN_DIR/"integrated_figure/latest_overview_revision.json").read_text())["revision"])\n'),
        code('show("Fig_StromalAtlas_Overview.png", 1800)'),
        code('for name in ["Panel_A_fibro_Leiden_subtype.png", "Panel_B_endo_Leiden_subtype.png", "Panel_C_mirrored_coverage.png"]:\n    show(name)\n'),
        code('for letter, lineage, suffix in [("E","fibro","markers"),("F","endo","markers"),("G","fibro","signatures"),("H","endo","signatures"),("I","fibro","verified_metacells"),("J","endo","verified_metacells")]:\n    show(f"Panel_{letter}_{lineage}_{suffix}.png", 1100)\n'),
        code('display(json.loads((RESULT/"trajectory_input_verification.json").read_text()))\ndisplay(pd.read_csv(RESULT/"C_tissue_mirrored_counts.csv"))\ndisplay(pd.read_csv(RESULT/"C_age_group_mirrored_counts.csv"))\nprint((RESULT/"figure_caption.txt").read_text())\n'),
    ]
    nbf.validate(notebook)
    nbf.write(notebook, revision / 'Replot_Overview.ipynb')
    nbf.write(notebook, pipeline.HERE / 'replot_overview.ipynb')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=pipeline.HERE.parent / 'results/runs/meta_full_20261031_232542_2002510')
    parser.add_argument('--plot-config', type=Path, default=pipeline.HERE / 'overview_replot_config.json')
    parser.add_argument('--dpi', type=int)
    args = parser.parse_args()
    root = args.run_dir.resolve()
    if (root / '.running.lock').exists():
        raise RuntimeError('An analysis process is writing this run; finish it before replotting')
    manifest = json.loads((root / 'run_manifest.json').read_text())
    config = manifest['config']
    for lineage in ['fibro', 'endo']:
        for stage in ['atlas', 'markers', 'trajectory']:
            if manifest['stages'].get(f'{stage}:{lineage}', {}).get('status') != 'complete':
                raise ValueError(f'{stage}:{lineage} is incomplete')
    options = json.loads(args.plot_config.read_text())
    if args.dpi:
        options['dpi'] = args.dpi
    stamp = dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).strftime('%Y%m%d_%H%M%S_%f')
    integrated = root / 'integrated_figure'
    revision = integrated / 'revisions' / ('overview_v2_' + stamp)
    revision.mkdir(parents=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',
        handlers=[logging.StreamHandler(), logging.FileHandler(revision / 'replot.log')], force=True)
    from overview_v2 import render_overview
    logging.info('Rendering revision %s', revision)
    render_overview(config, root, revision, options)
    write_notebook(root, revision)
    # Meaningful checks of biological identity and coverage precede publishing the image.
    import pandas as pd
    coverage = pd.read_csv(revision / 'C_tissue_mirrored_counts.csv')
    for lineage in ['fibro', 'endo']:
        audit = json.loads((root / 'shared' / f'data_audit_{lineage}.json').read_text())
        if int(coverage.loc[coverage.lineage == lineage, 'n_cells'].sum()) != audit['n_cells']:
            raise ValueError('Merged coverage lost cells')
    sources = ['overview_v2.py', '06_replot_overview.py', 'overview_replot_config.json',
               'run_replot_overview.sh', 'utils_plot.py', 'utils_data.py', 'metacell_trajectory.py',
               'README_overview_replot.md', 'reference_atlas_user_code.txt',
               'stromal_atlas_pipeline.py', 'requirements.txt', 'requirements_trajectory.txt']
    snapshot = revision / 'source_code'
    snapshot.mkdir()
    for name in sources:
        source = pipeline.HERE / name
        if source.exists():
            shutil.copy2(source, snapshot / name)
    (revision / 'plot_config_used.json').write_text(json.dumps(options, indent=2))
    # Back up the current original before replacing the exact user-requested filenames.
    previous = integrated / 'history' / ('before_' + stamp)
    previous.mkdir(parents=True)
    for name in ['Fig_StromalAtlas_Overview.png', 'Fig_StromalAtlas_Overview.pdf', 'figure_caption.txt',
                 'trajectory_input_verification.json', 'overview_plot_provenance.json']:
        target = integrated / name
        if target.exists():
            shutil.copy2(target, previous / name)
        if (revision / name).exists():
            temporary = integrated / ('.' + name + '.tmp')
            shutil.copy2(revision / name, temporary)
            temporary.replace(target)
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in snapshot.iterdir()}
    record = dict(status='passed', revision=str(revision), previous_overview=str(previous),
                  analysis_refitted=False, source_code_sha256=files,
                  checks=['Original meta-cell paths and unchanged source file fingerprints',
                          'Meta-cell IDs, order, annotation and original UMAP alignment',
                          'Full-cell atlas IDs and UMAP alignment',
                          'Same actual Leiden palette and all atlas cells drawn',
                          'Merged tissue coverage cell totals preserved'])
    (revision / 'replot_validation.json').write_text(json.dumps(record, indent=2))
    (integrated / 'latest_overview_revision.json').write_text(json.dumps(record, indent=2))
    logging.info('SUCCESS. Updated overview: %s', integrated / 'Fig_StromalAtlas_Overview.png')
    print('RESULT_DIR=' + str(revision))
    print('NOTEBOOK=' + str(revision / 'Replot_Overview.ipynb'))


if __name__ == '__main__':
    main()
