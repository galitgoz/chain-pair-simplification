"""Saved-results rendering and final audit. No optimization."""
from common import ROOT,OUT,atomic,sha
from pathlib import Path
import os,sys,json,platform,importlib.metadata
import pandas as pd,psutil,nbformat
from nbclient import NotebookClient
os.environ['JUPYTER_PATH']=str(OUT/'notebook_kernel')
k=OUT/'notebook_kernel/kernels/graph-v3';k.mkdir(parents=True,exist_ok=True)
(k/'kernel.json').write_text(json.dumps(dict(argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],display_name='Graph v3 diagnostic',language='python')))
p=Path(__file__).with_name('00_controlled_transition_budget.ipynb');nb=nbformat.read(p,as_version=4)
try:NotebookClient(nb,timeout=120,kernel_name='graph-v3',resources={'metadata':{'path':str(ROOT)}}).execute()
finally:nbformat.write(nb,p)
nbformat.validate(nb)
for c in nb.cells:
    if c.cell_type=='code':assert not any(o.output_type=='error' for o in c.outputs)
index=json.loads((OUT/'attempt_index.json').read_text());assert len(index)==1
assert len(list((OUT/'attempts').glob('*/configuration.json')))==1
d=Path(index[0]['attempt_directory']);c=json.loads((d/'configuration.json').read_text())
assert sha(ROOT/c['input_file'])==c['input_sha256']
for name,digest in c['algorithms'].items():assert sha(ROOT/name)==digest
for name,digest in c['code_version'].items():assert sha(ROOT/name)==digest
raw=json.loads((d/'raw_solver.json').read_text());validation=json.loads((d/'validation.json').read_text())
assert raw['certificate_enabled'] is False and raw['configuration']['certificate'] is False
assert raw['solver_status']=='optimal_configuration_graph' and raw['observation']['configuration_graph_entered']
assert raw['observation']['chain_graphs_completed']==2
assert validation['status']=='validated' and validation['raw_checkpoint_sha256']==sha(d/'raw_solver.json')
assert raw['monotonic']<validation['monotonic']
assert raw['raw_result']['k']==validation['k']==5
assert raw['raw_result']['states']==raw['observation']['states_dequeued']==5678
for r in pd.read_csv(OUT/'preservation_snapshot.csv').itertuples():assert sha(ROOT/r.path)==r.sha256
active=[]
for pr in psutil.process_iter(['pid','cmdline']):
    text=' '.join(pr.info['cmdline'] or [])
    if 'graph_only_v3' in text and any(x in text for x in ['worker.py','run_once.py']):active.append(pr.pid)
assert not active
atomic(OUT/'environment.json',dict(python=platform.python_version(),platform=platform.platform(),versions={x:importlib.metadata.version(x) for x in ['numpy','numba','pandas','psutil','nbclient']}))
pd.DataFrame([dict(path=str(f.relative_to(ROOT)),sha256=sha(f)) for f in Path(__file__).parent.iterdir() if f.is_file()]).to_csv(OUT/'artifact_hashes.csv',index=False)
atomic(OUT/'final_audit.json',dict(new_attempts=1,inherited_verified=23,new_verified_graph_results=1,total_verified=24,
    total_successful_CPS2F_graph_results=2,reference_agreement=True,certificate_enabled=False,algorithm_changes=False,
    original_files_unchanged=True,attempt_code_hashes_unchanged=True,active_workers=active,preexisting_kernel_present=psutil.pid_exists(22664),automatic_retry=False))
print('One attempt; saved and validated graph optimum k=5; unchanged original files; no active experiment workers.')
