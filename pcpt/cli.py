"""Command-line entrypoint for auditable PCPT research operations."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from .field import evaluate, write_result, DataQualityError
from .science import Scenario, calculate


def main(argv=None):
    parser = argparse.ArgumentParser(prog='pcpt', description='PCPT v3.0 — GRM research (no claims of natural validation)')
    sub = parser.add_subparsers(dest='command', required=True)
    model = sub.add_parser('model', help='Run analytical SYNTHETIC GRM benchmark')
    model.add_argument('--scenario', help='JSON scenario, defaults to documented synthetic baseline')
    model.add_argument('--out', default='outputs/model_summary.json')
    benchmark = sub.add_parser('benchmark', help='Run original finite-volume benchmark + graphs')
    field = sub.add_parser('field', help='Screen a real or synthetic observation CSV with provenance checks')
    field.add_argument('input_csv'); field.add_argument('--period-s', type=float, default=44712.0)
    field.add_argument('--out', default='outputs/field_screen.json')
    args = parser.parse_args(argv)
    if args.command == 'benchmark':
        import grm_benchmark
        grm_benchmark.main()
        return 0
    if args.command == 'model':
        result = calculate(Scenario.read(args.scenario) if args.scenario else Scenario())
        write_result(result, args.out)
    else:
        try:
            result = evaluate(args.input_csv, period_s=args.period_s)
        except (DataQualityError, FileNotFoundError) as exc:
            parser.error(str(exc))
        write_result(result, args.out)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
