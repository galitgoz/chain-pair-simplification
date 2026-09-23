"""Recompute matched paper summaries from saved records; never run solvers."""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/primary.csv'
METHODS = {'Independent continuous', 'Independent discrete', 'CPS-2F', 'CPS-3F'}


def summarize():
    with SOURCE.open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    groups = defaultdict(dict)
    for row in rows:
        if row['validated'].lower() != 'true' or float(row['alpha']) != 0.5:
            continue
        key = (row['domain'], row['pair'])
        method = row['method']
        if method not in METHODS or method in groups[key]:
            raise ValueError(f'Unexpected or duplicate method for {key}: {method}')
        groups[key][method] = row
    common = {key: value for key, value in groups.items() if set(value) == METHODS}
    for key, methods in common.items():
        signatures = {
            (r['input_sha256'], *(float(r[field]) for field in
              ['nA', 'nB', 'delta1', 'delta2', 'delta3']))
            for r in methods.values()
        }
        if len(signatures) != 1:
            raise ValueError(f'Methods do not share inputs and thresholds: {key}')
    counts = dict(Counter(domain for domain, pair in common))
    if counts != {'protein': 12, 'hurricane': 31}:
        raise ValueError(f'Unexpected common cohort: {counts}')
    observations = [r for key in sorted(common) for _, r in sorted(common[key].items())]
    for r in observations:
        if float(r['k']) != max(float(r['kA']), float(r['kB'])):
            raise ValueError(f'Inconsistent objective: {r["pair"]}')
        if r['solver_time_kind'] != 'solution_return':
            raise ValueError(f'Unexpected timing kind: {r["pair"]}')
    medians = []
    for domain in ['protein', 'hurricane']:
        for method in ['CPS-2F', 'CPS-3F']:
            selected = [r for r in observations if r['domain'] == domain and r['method'] == method]
            medians.append(dict(domain=domain, method=method, pairs=len(selected),
                                median_solver_seconds=median(float(r['solver_seconds']) for r in selected)))
    quality = Counter()
    for methods in common.values():
        ic, cps, discrete = (methods[m] for m in ['Independent continuous', 'CPS-2F', 'CPS-3F'])
        delta = float(cps['k']) - float(ic['k'])
        violation = float(ic['coupling_distance']) > float(ic['delta3']) + 1e-9
        quality['equal_K_IC_CPS2F'] += delta == 0
        quality['CPS2F_K_plus_one'] += delta == 1
        quality['IC_coupling_violations'] += violation
        quality['IC_violations_repaired_at_equal_K'] += violation and delta == 0
        quality['CPS2F_smaller_K_than_CPS3F'] += float(cps['k']) < float(discrete['k'])
        quality['CPS2F_coupling_satisfied'] += float(cps['coupling_distance']) <= float(cps['delta3']) + 1e-9
    summary = dict(common_pairs=len(common), common_observations=len(observations), domains=counts,
                   alpha=0.5, quality=dict(quality),
                   source=SOURCE.relative_to(ROOT).as_posix(), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                   timing='One validated primary observation per pair and method; solver_seconds; resource exits excluded.')
    return observations, medians, summary


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('reproduced/primary'))
    args = parser.parse_args()
    target = args.output_dir.resolve()
    if target == ROOT or any(target.is_relative_to(ROOT / name) for name in ['data', 'results', '.git']):
        parser.error('Choose a generated-output directory outside data/, results/, and .git/.')
    observations, medians, summary = summarize()
    target.mkdir(parents=True, exist_ok=True)
    write_csv(target / 'common_comparison.csv', observations)
    write_csv(target / 'solver_time_medians.csv', medians)
    (target / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))
    print(f'Summaries written to {target}')


if __name__ == '__main__':
    main()
