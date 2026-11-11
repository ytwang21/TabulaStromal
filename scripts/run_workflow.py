#!/usr/bin/env python3
"""Launch a workflow from a consistent repository-root working directory."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
ENTRIES = {
    'atlas': ROOT / 'workflows/atlas/code/stromal_atlas_pipeline.py',
    'benchmark': ROOT / 'workflows/benchmark/run_tabulastromal.py',
    'final-clock': ROOT / 'workflows/benchmark/run_final_clock.py',
    'cellage': ROOT / 'workflows/cellage/run_pipeline.py',
    'signatures': ROOT / 'workflows/atlas/extended_figure1_SS/scripts/run_signature_score.py',
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workflow', choices=ENTRIES)
    if not argv or argv[0] in ('-h', '--help'):
        parser.print_help()
        return 0
    args = parser.parse_args(argv[:1])
    forwarded = argv[1:]
    # Resolve caller-specified configs before changing cwd. Other relative paths
    # deliberately follow the documented repository-root convention.
    for i, value in enumerate(forwarded):
        if value == '--config' and i + 1 < len(forwarded):
            forwarded[i + 1] = str(Path(forwarded[i + 1]).resolve())
        elif value.startswith('--config='):
            forwarded[i] = '--config=' + str(Path(value.split('=', 1)[1]).resolve())
    env = dict(os.environ)
    env.setdefault('MPLBACKEND', 'Agg')
    env.setdefault('MPLCONFIGDIR', str(ROOT / 'outputs/runtime/matplotlib'))
    env.setdefault('NUMBA_CACHE_DIR', str(ROOT / 'outputs/runtime/numba'))
    try:
        return subprocess.call([sys.executable, str(ENTRIES[args.workflow]), *forwarded],
                               cwd=ROOT, env=env)
    except KeyboardInterrupt:
        return 130


if __name__ == '__main__':
    raise SystemExit(main())
