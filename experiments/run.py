"""Rerun saved primary configurations without modifying reference results."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--domain', choices=['protein', 'hurricane'])
    parser.add_argument('--pair', help='Exact pair identifier, for example EP111954__EP021984')
    parser.add_argument('--method', choices=['CPS-2F', 'CPS-3F', 'Independent continuous', 'Independent discrete'])
    parser.add_argument('--limit', type=int, help='Run the first N selected configurations')
    parser.add_argument('--list', action='store_true', help='List selections without running solvers')
    parser.add_argument('--output-dir', type=Path, default=ROOT/'reproduced/experiments')
    args = parser.parse_args()
    rows = json.loads((ROOT/'experiments/configurations.json').read_text())
    rows = [c for c in rows if all(not value or c[key] == value for key, value in
            [('domain', args.domain), ('pair', args.pair), ('method', args.method)])]
    if args.limit is not None:
        if args.limit < 1:
            parser.error('--limit must be positive')
        rows = rows[:args.limit]
    if not rows:
        parser.error('No configurations match the selection')
    if args.list:
        for c in rows:
            print(c['domain'], c['pair'], c['method'], c['reference_status'])
        return 0
    if sys.platform != 'win32':
        parser.error('The process supervisor currently requires Windows; saved-result reporting is cross-platform')
    out = args.output_dir.resolve()
    generated = (ROOT/'reproduced').resolve()
    if not out.is_relative_to(generated) or out == generated or out.exists():
        parser.error('Choose a new subdirectory of reproduced/; existing runs are never overwritten')
    out.mkdir(parents=True)
    os.environ['CPS_EXPERIMENT_OUTPUT'] = str(out)
    from supervisor import attempt
    algorithms = {f'src/{n}': hashlib.sha256((ROOT/'src'/n).read_bytes()).hexdigest()
                  for n in ['curve_algorithms.py', 'cps_paper_algorithms.py']}
    results = []
    for index, original in enumerate(rows, 1):
        c = dict(original)
        c['algorithms'] = algorithms
        c['diagnostic'] = f'primary_{index:03}'
        print(f"{index}/{len(rows)} {c['pair']} {c['method']}", flush=True)
        folder = attempt(c, output_directory=out)
        def read(name):
            p = folder/name
            return json.loads(p.read_text()) if p.exists() else {}
        validation, error, terminal = read('validation.json'), read('worker_error.json'), read('supervisor.json')
        sizes = {key: validation.get(key) for key in ['kA', 'kB', 'k']}
        matched = validation.get('status') == 'validated' and sizes == c['reference_sizes']
        status = error.get('status') or validation.get('status') or terminal.get('status')
        resource_exit = status in {'solver_timeout', 'transition_budget_limit', 'state_budget_limit',
                                  'chain_edge_budget_limit', 'memory_limit'}
        expected_exit = c['reference_status'] != 'validated' and resource_exit
        results.append(dict(pair=c['pair'], method=c['method'], status=status,
                            sizes=sizes, matches_reference_sizes=matched,
                            expected_resource_exit=expected_exit,
                            attempt=folder.relative_to(ROOT).as_posix()))
        (out/'comparison.json').write_text(json.dumps(results, indent=2)+'\n')
        print(f'  {status}; sizes match: {matched}; expected resource exit: {expected_exit}', flush=True)
        if status == 'startup_failure':
            print(error.get('reason', 'Worker startup failed'), file=sys.stderr)
            return 1
    return int(any(not (r['matches_reference_sizes'] or r['expected_resource_exit']) for r in results))


if __name__ == '__main__':
    raise SystemExit(main())
