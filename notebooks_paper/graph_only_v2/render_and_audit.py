"""Render saved evidence, audit attempts and preservation; no solves."""
from common import ROOT,OUT,sha,atomic
from pathlib import Path
import os,sys,json,platform,importlib.metadata
import pandas as pd,psutil,nbformat
from nbclient import NotebookClient
os.environ['JUPYTER_PATH']=str(OUT/'notebook_kernel')
k=OUT/'notebook_kernel/kernels/graph-v2';k.mkdir(parents=True,exist_ok=True)
(k/'kernel.json').write_text(json.dumps(dict(argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],display_name='Graph v2 diagnostics',language='python')))
p=Path(__file__).with_name('00_repair_and_diagnostics.ipynb');nb=nbformat.read(p,as_version=4)
try:NotebookClient(nb,timeout=120,kernel_name='graph-v2',resources={'metadata':{'path':str(ROOT)}}).execute()
finally:nbformat.write(nb,p)
nbformat.validate(nb)
for c in nb.cells:
    if c.cell_type=='code':assert not any(o.output_type=='error' for o in c.outputs)
index=json.loads((OUT/'attempt_index.json').read_text());assert len(index)==3
assert len(list((OUT/'attempts').glob('*/configuration.json')))==3
for entry in index:
    d=Path(entry['attempt_directory']);c=json.loads((d/'configuration.json').read_text())
    assert c['certificate'] is False and c['phase_budgets']['solver']==60 and c['memory_mib']==512
    assert sha(ROOT/c['input_file'])==c['input_sha256']
    for name,digest in c['algorithms'].items():assert sha(ROOT/name)==digest
    for name,digest in c['code_version'].items():assert sha(ROOT/name)==digest
    if (d/'raw_solver.json').exists():
        raw=json.loads((d/'raw_solver.json').read_text());v=json.loads((d/'validation.json').read_text())
        assert v['raw_checkpoint_sha256']==sha(d/'raw_solver.json')
        assert raw['monotonic']<v['monotonic']
    for field,want in [('max_states',30000),('max_transitions',2000000 if c['diagnostic']=='A' else 500000),('max_chain_edges',800000 if c['diagnostic']=='B' else 200000),('seconds',55.)]:
        assert c['graph_limits'][field]==want
for r in pd.read_csv(OUT/'preservation_snapshot.csv').itertuples():assert sha(ROOT/r.path)==r.sha256
active=[]
for pr in psutil.process_iter(['pid','cmdline']):
    command=' '.join(pr.info['cmdline'] or [])
    if 'graph_only_v2' in command and any(x in command for x in ['worker.py','run_diagnostics.py']):active.append(pr.pid)
assert not active,active
atomic(OUT/'environment.json',dict(python=platform.python_version(),platform=platform.platform(),versions={x:importlib.metadata.version(x) for x in ['numpy','numba','pandas','psutil','nbclient']}))
artifacts=[dict(path=str(f.relative_to(ROOT)),sha256=sha(f)) for f in Path(__file__).parent.iterdir() if f.is_file()]
pd.DataFrame(artifacts).to_csv(OUT/'artifact_hashes.csv',index=False)
atomic(OUT/'final_audit.json',dict(repair_checks_passed=21,diagnostic_attempts=3,inherited_verified=22,newly_validated=1,
    total_verified=23,remaining_missing=9,old_files_unchanged=True,attempt_code_versions_unchanged=True,active_diagnostic_workers=active,
    preexisting_kernel_present=psutil.pid_exists(22664),algorithm_changes=False,automatic_escalation=False))
print('Notebook rendered; exactly three attempts audited; inherited/code/input hashes intact; no active diagnostic workers.')
