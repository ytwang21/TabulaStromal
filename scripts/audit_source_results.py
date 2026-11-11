#!/usr/bin/env python3
"""Inventory the four original result trees without loading expression matrices.

The complete path-level inventory stays in ignored outputs/ because cache paths
may contain sample labels. Public summaries cover every scanned file; public
figure metadata exclude vendor, runtime, backup and source-snapshot copies.
"""
import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOTS = {
    'F1_StromalAtlas': 'F1_StromalAtlas',
    'F2_ModelPerformance': 'F2_ModelPerformance',
    'TabulaStromal_Rebenchmark_2026': 'AgePrediction/GeneModel/TabulaStromal_Rebenchmark_2026',
    'CellAge_Exploration': 'AgePrediction/CellAge_Exploration',
}
NONRESULT_PARTS = {'vendor', 'runtime', 'tmp', 'backups', 'source_snapshots', 'source_code',
                   '__pycache__', '.ipynb_checkpoints', 'smoke_tests', 'tests', 'Rlib'}


def kind(parts, suffix):
    if any(p in NONRESULT_PARTS for p in parts):
        return 'runtime_vendor_backup_or_test'
    if suffix in {'.joblib', '.pkl', '.pickle', '.npz', '.npy', '.pt', '.pth', '.parquet', '.h5ad'}:
        return 'model_checkpoint_or_private_data'
    if suffix in {'.png', '.pdf', '.svg', '.jpg', '.jpeg'}:
        return 'figure_candidate'
    if suffix in {'.csv', '.tsv', '.xlsx'}:
        return 'table_candidate'
    if suffix in {'.py', '.R', '.ipynb', '.yaml', '.json', '.md', '.txt', '.sh'}:
        return 'code_configuration_or_report'
    return 'other'


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    local = ROOT / 'outputs/source_audit'
    local.mkdir(parents=True, exist_ok=True)
    full = local / 'full_file_inventory.csv.gz'
    totals = defaultdict(lambda: [0, 0])
    figures = []
    n = 0
    fields = ['module', 'source_file', 'type', 'bytes', 'mtime_utc', 'expected_use', 'read_status']
    with gzip.open(full, 'wt', compresslevel=1, newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        for module, rel in SOURCE_ROOTS.items():
            source = project / rel
            if not source.is_dir():
                raise FileNotFoundError(source)
            for directory, dirs, files in os.walk(source, followlinks=False):
                links = [d for d in dirs if (Path(directory) / d).is_symlink()]
                dirs[:] = [d for d in dirs if d not in links]
                for name in files + links:
                    path = Path(directory) / name
                    try:
                        stat = path.lstat()
                    except FileNotFoundError:
                        continue
                    relative = path.relative_to(project)
                    suffix = path.suffix.lower()
                    category = kind(path.relative_to(source).parts, suffix)
                    row = dict(module=module, source_file=str(relative), type=suffix or '(none)',
                               bytes=stat.st_size,
                               mtime_utc=datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                               expected_use=category, read_status='metadata_only; file content not inferred')
                    writer.writerow(row)
                    totals[(module, category, suffix or '(none)')][0] += 1
                    totals[(module, category, suffix or '(none)')][1] += stat.st_size
                    if category == 'figure_candidate' and not path.is_symlink():
                        figures.append(row)
                    n += 1
                    if n % 200000 == 0:
                        print(json.dumps({'scanned_files': n}), flush=True)
            print(json.dumps({'completed_module': module, 'scanned_files': n}), flush=True)
    summary = [dict(module=k[0], category=k[1], extension=k[2], n_files=v[0], bytes=v[1])
               for k, v in sorted(totals.items())]
    write_csv(ROOT / 'audit/source_tree_summary.csv', summary,
              ['module', 'category', 'extension', 'n_files', 'bytes'])
    write_csv(local / 'figure_candidates.csv', figures, fields)
    report = dict(scanned_files=n, figure_candidates=len(figures),
                  complete_inventory='outputs/source_audit/full_file_inventory.csv.gz',
                  complete_inventory_public=False, source_roots=SOURCE_ROOTS,
                  policy='Full metadata scan; no expression matrix densified; no model content read by discovery.')
    (local / 'scan_receipt.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
