#!/usr/bin/env python3
"""Audit the public payload without reading private research inputs."""
import argparse
import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import subprocess
import sys
from zipfile import ZipFile
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
CHECKSUMS = 'provenance/export_checksums.csv'
EXCLUDED_FROM_HASH = {CHECKSUMS, 'provenance/validation.json', 'provenance/publication.json'}
PRIVATE_SUFFIXES = {'.h5ad', '.h5', '.hdf5', '.loom', '.joblib', '.pkl', '.pickle',
                    '.parquet', '.npz', '.npy', '.pt', '.pth', '.pem', '.key'}
SECRET_PATTERNS = [r'gh[pousr]_[A-Za-z0-9]{25,}', r'github_pat_[A-Za-z0-9_]{30,}',
                   r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',
                   r'AKIA[0-9A-Z]{16}', r'sk-[A-Za-z0-9]{35,}']


def public_files():
    result = subprocess.run(['git', 'ls-files', '-z', '--cached', '--others', '--exclude-standard'],
                            cwd=ROOT, text=True, capture_output=True, check=True)
    return sorted({p for p in result.stdout.split('\0') if p and (ROOT / p).is_file()})


def table(path):
    with (ROOT / path).open() as handle:
        return list(csv.DictReader(handle))


def close(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-8, abs_tol=1e-9)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-checksums', action='store_true')
    parser.add_argument('--verify-checksums', action='store_true')
    args = parser.parse_args(argv)
    errors = []
    counts = {'python': 0, 'notebooks': 0, 'structured': 0, 'markdown_links': 0}
    files = public_files()
    released_path=ROOT/'audit/released_assets.csv'
    released={r['repository_file']:r for r in table('audit/released_assets.csv')} if released_path.exists() else {}
    for rel in files:
        path = ROOT / rel
        if path.is_symlink():
            errors.append(f'Symlink in public payload: {rel}')
            continue
        if path.suffix in PRIVATE_SUFFIXES:
            errors.append(f'Excluded research/credential format: {rel}')
        if path.stat().st_size > 2_000_000:
            approved=released.get(rel)
            limit=32_000_000 if path.suffix in {'.png','.pdf'} else 8_000_000
            if not (approved and path.suffix in {'.png','.pdf','.csv'} and
                    path.stat().st_size <= limit and int(approved['bytes'])==path.stat().st_size and
                    hashlib.sha256(path.read_bytes()).hexdigest()==approved['sha256']):
                errors.append(f'Large file requires a matching reviewed scientific-asset record: {rel}')
        raw = path.read_bytes()
        if path.suffix in {'.png', '.jpg', '.jpeg', '.pdf'}:
            continue
        if path.suffix == '.docx':
            try:
                with ZipFile(path) as archive:
                    assert archive.testzip() is None
                    document=ET.fromstring(archive.read('word/document.xml'))
                    assert any(n.text for n in document.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t'))
            except Exception as exc:
                errors.append(f'Invalid editable Word file {rel}: {exc}')
            continue
        try:
            text = raw.decode('utf-8')
        except UnicodeDecodeError:
            errors.append(f'Unexpected non-text artifact: {rel}')
            continue
        for pattern in SECRET_PATTERNS:
            if re.search(pattern, text):
                errors.append(f'Potential credential pattern in {rel}; value not printed')
        if re.search(r'/(?:fs|data)/home/[A-Za-z0-9_]', text):
            errors.append(f'Personal absolute server path remains: {rel}')
        try:
            if path.suffix == '.py':
                ast.parse(text, filename=rel)
                counts['python'] += 1
            elif path.suffix == '.ipynb':
                nb = json.loads(text)
                assert nb['nbformat'] == 4 and isinstance(nb['cells'], list)
                for cell in nb['cells']:
                    if cell['cell_type'] == 'code':
                        assert cell.get('outputs') == [] and cell.get('execution_count') is None
                counts['notebooks'] += 1
            elif path.suffix == '.json':
                json.loads(text)
                counts['structured'] += 1
            elif path.suffix in {'.yml', '.yaml', '.cff'}:
                import yaml
                yaml.safe_load(text)
                counts['structured'] += 1
        except Exception as exc:
            errors.append(f'Invalid source/structured file {rel}: {type(exc).__name__}: {exc}')
        if path.suffix == '.md':
            targets=re.findall(r'\]\(([^\s)]+)\)', text)+re.findall(r'(?:src|href)="([^\"]+)"',text)
            for target in targets:
                if target.startswith(('https://', 'http://', 'mailto:', '#')):
                    continue
                target = target.split('#', 1)[0]
                if target and not (path.parent / target).exists():
                    errors.append(f'Broken Markdown link in {rel}: {target}')
                counts['markdown_links'] += 1

    for row in table('provenance/source_manifest.csv'):
        if not (ROOT / row['repository_path']).is_file():
            errors.append('Missing exported source: ' + row['repository_path'])
        if not re.fullmatch('[0-9a-f]{64}', row['source_sha256']):
            errors.append('Malformed source hash: ' + row['source_path'])

    try:
        summary = table('results/final_clock_results_summary.csv')
        assert {r['lineage'] for r in summary} == {'Fibro', 'Endo'} and len(summary) == 2
        metric_comparisons = 0
        for row in summary:
            lin = row['lineage']
            features = table(f'results/{lin}/stability/stable_Celltype_Gene.csv')
            assert len(features) == int(row['n_stable_Celltype_Gene'])
            assert len({f['Gene'] for f in features}) == int(row['n_unique_stable_Genes'])
            assert len({f['Celltype'] for f in features}) == int(row['n_stable_Celltypes'])
            assert all(float(f['selection_frequency']) >= float(row['threshold']) for f in features)
            assert row['selected_version'] == 'original_EN'
            assert float(row['internal_original_MAE']) < float(row['internal_stable_MAE'])
            assert float(row['noadapt_R2']) < 0
            assert float(row['FT30_all100_MAE']) < float(row['FT30_heldout70_MAE'])
            fold_rows = table(f'results/{lin}/benchmark/cv_fold_metrics.csv')
            aggregates = table(f'results/{lin}/benchmark/cv_summary.csv')
            assert len(fold_rows) == 13 * 4 * 5 * 2
            assert len(aggregates) == 13 * 4 * 2
            for agg in aggregates:
                matching = [f for f in fold_rows if all(f[k] == agg[k] for k in ['lineage', 'K', 'model', 'split'])]
                assert len(matching) == 5
                for metric in ['R2', 'MAE', 'MAD', 'MedAE', 'RMSE', 'PCC']:
                    values = [float(f[metric]) for f in matching]
                    assert close(statistics.mean(values), agg[metric + '_mean'])
                    assert close(statistics.stdev(values), agg[metric + '_std'])
                    metric_comparisons += 2
        counts['aggregate_metric_comparisons'] = metric_comparisons
        atlas = table('results/atlas_summary.csv')
        assert {int(r['n_cells']) for r in atlas} == {5916760, 2239277}
        assert all(r['independent_donor_count'] == 'not_verified' for r in atlas)
    except Exception as exc:
        errors.append(f'Aggregate evidence consistency failed: {type(exc).__name__}: {exc}')

    if args.write_checksums and not errors:
        with (ROOT / CHECKSUMS).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=['path', 'bytes', 'sha256'], lineterminator='\n')
            writer.writeheader()
            for rel in files:
                if rel not in EXCLUDED_FROM_HASH:
                    raw = (ROOT / rel).read_bytes()
                    writer.writerow(dict(path=rel, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    if args.verify_checksums and not errors:
        if not (ROOT / CHECKSUMS).exists():
            errors.append('No exported-file checksum manifest. Generate with --write-checksums.')
        else:
            records = table(CHECKSUMS)
            expected = set(files) - EXCLUDED_FROM_HASH
            if {r['path'] for r in records} != expected:
                errors.append('Checksum manifest does not match the public file set.')
            for record in records:
                path = ROOT / record['path']
                if not path.is_file():
                    errors.append('Missing checksummed file: ' + record['path'])
                    continue
                raw = path.read_bytes()
                if len(raw) != int(record['bytes']) or hashlib.sha256(raw).hexdigest() != record['sha256']:
                    errors.append('Checksum differs: ' + record['path'])
    report = dict(status='passed' if not errors else 'failed', files=len(files),
                  payload_bytes=sum((ROOT / f).stat().st_size for f in files),
                  checks=counts, errors=errors)
    print(json.dumps(report, indent=2))
    return 1 if errors else 0


if __name__ == '__main__':
    raise SystemExit(main())
