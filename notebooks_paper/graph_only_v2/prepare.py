"""Read-only reconciliation, inherited provenance, and frozen diagnostic limits."""
from common import ROOT,OUT,sha,atomic,BASE_GRAPH_LIMITS,BUDGETS
import pandas as pd,json,shutil,psutil
OLD=OUT.parent/'graph_only_v1'
assert not (OUT/'diagnostics.json').exists(),'Already initialized'
OUT.mkdir(parents=True,exist_ok=True)
for directory in ['cache','attempts','regressions','references']:(OUT/directory).mkdir(exist_ok=True)
r=pd.read_csv(OLD/'pilot_results.csv');v=pd.read_csv(OLD/'validation.csv')
assert len(v)==100 and v.outcome.eq('PASS').all()
assert len(r)==32 and r.status.eq('optimal').sum()==22
j=r[r.method.str.startswith('CPS') & r.status.eq('optimal')]
assert len(j)==8 and j.raw_status.eq('optimal_configuration_graph').all() and j.configuration_graph_entered.all()
missing=r[r.status.ne('optimal')].copy()
def cause(row):
    s=str(row.failure_reason)
    if 'product transition' in s:return 'transition_budget_limit'
    if 'single-chain graph edge' in s:return 'single_chain_edge_budget_limit'
    if row.status=='time_limit':return 'solver_or_measurement_timeout_unsplit_v1'
    return str(row.status)
counts=missing.apply(cause,axis=1).value_counts().to_dict()
assert counts=={'transition_budget_limit':5,'solver_or_measurement_timeout_unsplit_v1':3,'single_chain_edge_budget_limit':2},counts
catalog=[]
for p in OUT.parent.rglob('*results.csv'):
    if OUT in p.parents:continue
    catalog.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),mtime=p.stat().st_mtime))
newer=[x for x in catalog if 'graph_only_' in x['path'] and 'graph_only_v1' not in x['path']]
assert not newer,('Newer run requires reconciliation',newer)
active=[]
for p in psutil.process_iter(['pid','ppid','cmdline']):
    command=' '.join(p.info['cmdline'] or [])
    if 'python' in command.lower() and any(x in command for x in ['ipykernel','run_pilot.py','pilot_worker.py']):active.append(p.info)
atomic(OUT/'reconciliation.json',dict(checks_passed=100,verified_results=22,successful_graph_CPS=8,missing=10,missing_causes=counts,newer_graph_results=newer,result_catalog=catalog,active_processes=active,action='read-only; existing kernel untouched',note='v1 timeout budget included solve and post-solver measurements; the files alone do not isolate the timed-out phase'))
inherited=r[r.status.eq('optimal')].copy()
inherited['provenance_file']='output/paper_v1/graph_only_v1/pilot_results.csv'
inherited['provenance_sha256']=sha(OLD/'pilot_results.csv')
inherited['new_execution']=False;inherited['new_timing_repetition']=False
inherited.to_csv(OUT/'inherited_verified_results.csv',index=False)
for domain in ['protein','hurricane']:shutil.copy2(OLD/f'{domain}_paths.json',OUT/f'inherited_{domain}_paths.json')
for name in ['input_manifest.csv','pilot_parameters.csv']:shutil.copy2(OLD/name,OUT/name)
params=pd.read_csv(OUT/'pilot_parameters.csv');manifest=pd.read_csv(OUT/'input_manifest.csv')
configs=[]
for id,pair,alpha,method,change in [('A','AL031854__AL011991',.5,'CPS-2F',{'max_transitions':2000000}),('B','1o7j.a__1hfj.c',1.,'CPS-2F',{'max_chain_edges':800000}),('C','1o7j.a__1d9q.d',.5,'Independent discrete',{})]:
    row=params[params.pair.eq(pair)&params.alpha.eq(alpha)].iloc[0].to_dict()
    inp=manifest[manifest.pair.eq(pair)].iloc[0]
    assert sha(ROOT/row['input_file'])==inp.input_sha256
    row.update(diagnostic=id,method=method,certificate=False,graph_limits={**BASE_GRAPH_LIMITS,**change},phase_budgets=BUDGETS,
               memory_mib=512,input_sha256=inp.input_sha256,parameter_file_sha256=sha(OUT/'pilot_parameters.csv'),
               graph_limits_applicable=method.startswith('CPS'),internal_seconds_applicable=method.startswith('CPS'),
               independent_discrete_chain_edge_default=2000000 if method=='Independent discrete' else None,
               algorithms={x:sha(ROOT/x) for x in ['src/cps_paper_algorithms.py','src/curve_algorithms.py']})
    configs.append(row)
atomic(OUT/'diagnostics.json',configs)
snapshot=[]
for base in [ROOT/'notebooks_paper',ROOT/'output']:
    for p in base.rglob('*'):
        if not p.is_file() or OUT in p.parents or ROOT/'notebooks_paper/graph_only_v2' in p.parents:continue
        if any(x in p.parts for x in ['cache','__pycache__']):continue
        snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
for x in ['src/cps_paper_algorithms.py','src/curve_algorithms.py']:snapshot.append(dict(path=x,sha256=sha(ROOT/x)))
pd.DataFrame(snapshot).to_csv(OUT/'preservation_snapshot.csv',index=False)
(OUT/'status.md').write_text('# Graph-only v2: repair and three diagnostics\n\nReconciled v1: 100 passing checks; 22 verified results, including 8 graph CPS results; 10 missing (5 transition limits, 2 chain-edge limits, 3 unsplit external timeouts). No newer graph results found. Inherited results are not new executions. Regression and reference checks precede the three diagnostic attempts.\n')
print(json.dumps(dict(verified=22,CPS_graph=8,checks=100,missing_causes=counts,diagnostics=3),indent=2))
