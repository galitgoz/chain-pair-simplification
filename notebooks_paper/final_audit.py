from paper_core import ROOT,OUT,sha
import nbformat,pandas as pd,json,psutil
notebooks=[]
for p in sorted((ROOT/'notebooks_paper').glob('*.ipynb')):
    nb=nbformat.read(p,as_version=4);nbformat.validate(nb)
    for c in nb.cells:
        if c.cell_type=='code':
            compile(c.source,str(p),'exec')
            assert not any(o.output_type=='error' for o in c.outputs),(p,c.outputs)
    notebooks.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),code_cells=sum(c.cell_type=='code' for c in nb.cells)))
assert len(notebooks)==4
active=[]
for p in psutil.process_iter(['pid','cmdline']):
    if any('pilot_worker.py' in x or 'execute_stage1.py' in x for x in (p.info['cmdline'] or [])):active.append(p.pid)
assert not active,active
results=pd.read_csv(OUT/'pilot_results.csv');assert len(results)==32
assert pd.read_csv(OUT/'validation.csv').outcome.eq('PASS').all()
records=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in sorted((ROOT/'notebooks_paper').glob('*.py'))]
pd.DataFrame(records).to_csv(OUT/'stage1_modules.csv',index=False)
(OUT/'final_audit.json').write_text(json.dumps(dict(notebooks=notebooks,active_stage1_workers=active,
    preexisting_kernel_still_present=psutil.pid_exists(22664),optimal=int(results.status.eq('optimal').sum()),
    missing=int(results.status.ne('optimal').sum()),optimization_attempts=32,measured_solver_starts=int(results.measured_execution_count.sum()),
    later_stages_executed=False),indent=2))
print('Four valid notebooks; no notebook errors; no remaining Stage 1 workers; 14 optimal, 18 missing.')
