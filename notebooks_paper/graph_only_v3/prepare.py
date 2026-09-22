"""Freeze one v3 follow-up and reuse v2 orchestration without touching previous files."""
from pathlib import Path
import sys,json,hashlib,shutil,datetime
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).parent
OUT=ROOT/'output/paper_v1/graph_only_v3';V2=OUT.parent/'graph_only_v2'
def sha(p):
    from repository_paths import resolve_recorded_path
    return hashlib.sha256(resolve_recorded_path(p).read_bytes()).hexdigest()
if (OUT/'launch.json').exists():raise SystemExit('An attempt already exists: inspect it; do not duplicate.')
assert not (OUT/'configuration.json').exists(),'Prepared version exists; inspect before proceeding.'
OUT.mkdir(parents=True,exist_ok=True)
for name in ['cache','attempts']:(OUT/name).mkdir(exist_ok=True)
import pandas as pd,psutil,numpy as np
catalog=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p),mtime=p.stat().st_mtime) for p in OUT.parent.rglob('*results.csv') if OUT not in p.parents]
newer=[x for x in catalog if 'graph_only_v' in x['path'] and not any(v in x['path'] for v in ['graph_only_v1','graph_only_v2'])]
assert not newer,('Reconcile newer results first',newer)
coverage=pd.read_csv(V2/'coverage_with_provenance.csv');assert len(coverage)==32 and coverage.status.eq('optimal').sum()==23
v1=pd.read_csv(OUT.parent/'graph_only_v1/pilot_results.csv')
assert len(v1[v1.method.eq('CPS-2F') & v1.status.eq('optimal') & v1.route.eq('configuration_graph')])==1
prior=pd.read_csv(V2/'diagnostic_summary.csv');old=prior[prior.diagnostic.eq('A')].iloc[0]
assert old.limiting_resource=='transition_budget_limit' and old.limiting_counter==2000001 and old.states_at_limit==5452
assert not old.raw_solution_saved and abs(old.solver_seconds-2.109742799992091)<1e-8
config=next(c for c in json.loads((V2/'diagnostics.json').read_text()) if c['diagnostic']=='A')
before=json.loads(json.dumps(config));config['graph_limits']['max_transitions']=10000000
assert config['method']=='CPS-2F' and config['pair']=='AL031854__AL011991' and config['alpha']==.5
assert config['graph_limits']==dict(seconds=55.,max_states=30000,max_transitions=10000000,max_chain_edges=200000)
assert config['phase_budgets']['solver']==60 and config['memory_mib']==512 and config['certificate'] is False
assert sha(ROOT/config['input_file'])==config['input_sha256']
data=np.load(ROOT/config['input_file']);assert data['A'].shape==data['B'].shape==(25,2)
assert sha(V2/'pilot_parameters.csv')==config['parameter_file_sha256']
for name,digest in config['algorithms'].items():assert sha(ROOT/name)==digest
for key in before:
    if key!='graph_limits':assert config[key]==before[key]
for key in before['graph_limits']:
    if key!='max_transitions':assert config['graph_limits'][key]==before['graph_limits'][key]
(OUT/'configuration.json').write_text(json.dumps(config,indent=2))
shutil.copy2(V2/'coverage_with_provenance.csv',OUT/'inherited_coverage.csv')
inherited=coverage[coverage.status.eq('optimal')].copy()
inherited['source_coverage_file']=str((V2/'coverage_with_provenance.csv').relative_to(ROOT))
inherited['source_coverage_sha256']=sha(V2/'coverage_with_provenance.csv')
inherited['new_execution_in_v3']=False;inherited['new_timing_repetition_in_v3']=False
inherited.to_csv(OUT/'inherited_verified_results.csv',index=False)
shutil.copy2(V2/'references/A.json',OUT/'reference_A.json')  # never read by the solver
for name in ['input_manifest.csv','pilot_parameters.csv']:shutil.copy2(V2/name,OUT/name)
copies=[]
for name in ['common.py','worker.py','supervisor.py']:
    original=ROOT/'notebooks_paper/graph_only_v2'/name;target=HERE/name
    assert not target.exists()
    text=original.read_text(encoding='utf-8')
    if name=='common.py':text=text.replace("OUT=ROOT/'output/paper_v1/graph_only_v2'","OUT=ROOT/'output/paper_v1/graph_only_v3'")
    if name=='worker.py':text=text.replace('checkpoint=dict(schema_version=2,','checkpoint=dict(schema_version=3,certificate_enabled=False,')
    target.write_text(text,encoding='utf-8')
    copies.append(dict(file=name,parent_file=str(original.relative_to(ROOT)),parent_sha256=sha(original)))
(OUT/'orchestration_parent.json').write_text(json.dumps(copies,indent=2))
snapshot=[]
for base in [ROOT/'notebooks_paper',ROOT/'output']:
    for p in base.rglob('*'):
        if not p.is_file() or OUT in p.parents or HERE in p.parents:continue
        if any(part in p.parts for part in ['cache','__pycache__']):continue
        snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
for p in list(ROOT.glob('*.py'))+list(ROOT.glob('*.ipynb')):snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
pd.DataFrame(snapshot).to_csv(OUT/'preservation_snapshot.csv',index=False)
active=[]
for p in psutil.process_iter(['pid','ppid','cmdline']):
    c=' '.join(p.info['cmdline'] or [])
    if 'python' in c.lower() and any(s in c for s in ['ipykernel','run_diagnostics.py','worker.py','run_pilot.py']):active.append(p.info)
(OUT/'reconciliation.json').write_text(json.dumps(dict(verified_coverage=23,total_configurations=32,successful_CPS2F_graph_results=1,prior_A=old.to_dict(),newer_graph_versions=newer,result_catalog=catalog,active_processes=active,action='read-only; no existing run stopped'),indent=2))
(OUT/'status.md').write_text('# Graph-only v3\n\nPrepared one attempt: A, CPS-2F, AL031854__AL011991, alpha=0.5. Only max_transitions changes (2,000,000 to 10,000,000). Inherited verified coverage 23/32; no previous v3 run found. All previous files preserved.\n')
print('Verified 23/32; no v3 or newer graph results. Exact 25+25 input and parameter hashes verified. Only transition cap changed.')
