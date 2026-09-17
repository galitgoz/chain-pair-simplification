from paper_core import ROOT,OUT,sha
import pandas as pd,json,platform,importlib.metadata
rows=[]
for directory in ['cps_papers_w8','hurricane_cps_expanded','hurricane_auxiliary_expanded']:
    p=ROOT/'output'/directory/'comparison.csv'
    if p.exists():
        d=pd.read_csv(p);status=d['status'].astype(str)
        rows.append(dict(experiment=directory,coverage='partial',rows=len(d),optimal=int(status.str.startswith('optimal').sum()),
            nonoptimal=int((~status.str.startswith('optimal')).sum()),note='Legacy three-method protocol; independent discrete reporting missing; not reused as pilot results'))
rows.extend([dict(experiment='Stage 1 pilot',coverage='planned',rows=32,note='Two pairs per domain, two alphas, four methods; validation gates execution'),
             dict(experiment='04 graph cost and 05 paper summary',coverage='planned',note='Not authorized for execution in Stage 1')])
pd.DataFrame(rows).to_csv(OUT/'prior_experiments.csv',index=False)
snapshot=[]
for p in (ROOT/'output').rglob('*'):
    if p.is_file() and OUT not in p.parents:snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
pd.DataFrame(snapshot).to_csv(OUT/'protected_results_snapshot.csv',index=False)
versions={name:importlib.metadata.version(name) for name in ['numpy','pandas','numba','matplotlib','nbclient','nbformat','psutil']}
(OUT/'environment.json').write_text(json.dumps(dict(python=platform.python_version(),platform=platform.platform(),versions=versions),indent=2))
(OUT/'logs/existing_run_inspection.json').write_text(json.dumps(dict(existing_kernel_pids=[11420,22664],observed_cpu_seconds=[483.703125,483.703125],observation='Two read-only process snapshots approximately a minute apart; no separate experiment runner; existing kernel not attached to or stopped',pilot_isolation='separate kernel; sequential below-normal-priority worker processes'),indent=2))
print('Protected old result artifacts:',len(snapshot))
