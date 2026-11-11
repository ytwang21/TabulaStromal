#!/usr/bin/env python3
"""Create ignored configs pointing to a local copy of the original project."""
import argparse
import copy
import json
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]


def localize(value, project):
    if isinstance(value, dict):
        return {k: localize(v, project) for k, v in value.items()}
    if isinstance(value, list):
        return [localize(v, project) for v in value]
    if isinstance(value, str):
        for prefix in ('./data/', 'data/', './AgePrediction/', 'AgePrediction/'):
            if value.startswith(prefix):
                return str(project / value.removeprefix('./'))
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, required=True)
    parser.add_argument('--force', action='store_true', help='Replace previously generated local configs')
    args = parser.parse_args(argv)
    project = args.project_root.resolve()
    if not project.is_dir():
        parser.error('Project root does not exist.')
    out = ROOT / 'configs.local'
    out.mkdir(exist_ok=True)
    configs = {}
    for name, rel in [('benchmark', 'workflows/benchmark/config.yaml'),
                      ('cellage', 'workflows/cellage/config.yaml')]:
        cfg = localize(yaml.safe_load((ROOT / rel).read_text()), project)
        cfg['output_root'] = str(ROOT / 'outputs' / name)
        if name == 'cellage':
            cfg['project_root'] = str(project)
        configs[out / f'{name}.yaml'] = yaml.safe_dump(cfg, sort_keys=False)
    cfg = json.loads((ROOT / 'workflows/atlas/code/config.json').read_text())
    for field in ('inputs', 'metacell_inputs'):
        cfg[field] = {k: str(project / v.removeprefix('../../../')) for k, v in cfg[field].items()}
    configs[out / 'atlas.json'] = json.dumps(cfg, indent=2) + '\n'
    conflicts = [str(p) for p in configs if p.exists() and not args.force]
    if conflicts:
        parser.error('Configs exist; use --force only to replace them: ' + ', '.join(conflicts))
    for path, text in configs.items():
        path.write_text(text)
        print(path.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
