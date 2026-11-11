#!/usr/bin/env python3
"""Apply a trusted final-EN bundle to H5AD or named pseudobulk CSV inputs."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', required=True, type=Path,
                        help='Trusted final_clock.joblib bundle (weights are not included)')
    parser.add_argument('--input', required=True, type=Path,
                        help='H5AD, pseudobulk artifact directory, or CSV with sample index first')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--expression', choices=['X', 'raw'], default='X')
    args = parser.parse_args(argv)
    if args.output.resolve() in {args.input.resolve(), args.model.resolve()}:
        parser.error('Output must differ from the input and model.')
    sys.path.insert(0, str(ROOT / 'workflows/benchmark'))
    from final_clock_common import load_final_clock, predict_age
    import pandas as pd
    source = pd.read_csv(args.input, index_col=0) if args.input.suffix.lower() == '.csv' else args.input
    result = predict_age(source, load_final_clock(args.model), expression=args.expression)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f'Wrote {len(result)} sample predictions to {args.output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
