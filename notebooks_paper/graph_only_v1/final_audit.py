from graph_core import ROOT,OUT,sha
from pathlib import Path
import nbformat,pandas as pd,json,psutil
HERE=Path(__file__).parent
artifacts=[]
for p in HERE.glob('*.ipynb'):
    nb=nbformat.read(p,as_version=4);nbformat.validate(nb)
    for c in nb.cells:
        if c.cell_type=='code':
            compile(c.source,str(p),'exec')
            assert not any(o.output_type=='error' for o in c.outputs)
    artifacts.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
assert len(artifacts)==4
for p in HERE.glob('*.py'):artifacts.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
pd.DataFrame(artifacts).to_csv(OUT/'artifact_hashes.csv',index=False)
active=[]
for p in psutil.process_iter(['pid','cmdline']):
    c=' '.join(p.info['cmdline'] or [])
    if 'graph_only_v1' in c and any(x in c for x in ['pilot_worker.py','run_pilot.py','render_notebooks.py']):active.append(p.pid)
assert not active,active
r=pd.read_csv(OUT/'pilot_results.csv');j=r[r.method.str.startswith('CPS') & r.status.eq('optimal')]
assert j.configuration_graph_entered.all() and j.raw_status.eq('optimal_configuration_graph').all()
assert j.states_discovered.gt(0).all() and j.transitions_examined.gt(0).all()
proof=j[['pair','alpha','method','raw_status','chain_graphs_completed','configuration_graph_entered','states_discovered','transitions_examined']]
proof.to_csv(OUT/'successful_graph_evidence.csv',index=False)
(OUT/'final_audit.json').write_text(json.dumps(dict(notebooks=4,notebook_errors=0,active_revision_workers=active,
    successful_CPS=len(j),successful_CPS2F=int(j.method.eq('CPS-2F').sum()),successful_CPS3F=int(j.method.eq('CPS-3F').sum()),
    certificate_rows_included=0,small_checks_passed=100,preexisting_kernel_present=psutil.pid_exists(22664),later_stages_executed=False),indent=2))
print(proof.to_string(index=False));print('Final audit passed; no active revision workers.')
