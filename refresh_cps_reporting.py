"""Recompute reporting while retaining unchanged, verified solver outputs."""
from pathlib import Path
import sys,json,subprocess,hashlib
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

root=Path(__file__).resolve().parent
path=root/'0.analyze_CPS.ipynb';old=nbformat.read(path,as_version=4)
out=root/'output/cps_papers';previous=json.loads((out/'manifest.json').read_text())
hashes={k.replace('\\','/'):v for k,v in previous['sha256'].items()}
for name in ['curve_algorithms.py','cps_paper_algorithms.py','cps_notebook_data.py',
             'cps_notebook_run.py','data/decimated protein backbones/decimated_curves.csv']:
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==hashes[name],name
assert previous['ws']==[16] and previous['alphas']==[1,2,4] and previous['input_scope']=='supplied'
(out/'solver_run_manifest.json').write_text(json.dumps(previous,indent=2))
subprocess.run([sys.executable,str(root/'build_cps_notebook.py')],cwd=root,check=True)
nb=nbformat.read(path,as_version=4)
saved={c.source:c for c in old.cells if c.cell_type=='code'}
for c in nb.cells:
    if c.cell_type=='code' and c.source in saved:
        c.outputs=saved[c.source].outputs;c.execution_count=saved[c.source].execution_count
start=next(i for i,c in enumerate(nb.cells) if c.source.startswith('### 6.1'))
demo=next(c for c in nb.cells if c.source.startswith('window_A'))
bootstrap=nbformat.v4.new_code_cell('''
instances = data['instances']
example = next(iter(instances.values()))
results = pd.read_csv(OUT/'comparison.csv')
paths = json.loads((OUT/'paths.json').read_text())
checks = json.loads((OUT/'algorithm_checks.json').read_text())
from curve_algorithms import discrete_frechet
from cps_paper_algorithms import solve_cps
''')
report=nbformat.v4.new_notebook(cells=[nb.cells[2],bootstrap,demo]+nb.cells[start:])
kernel_root=root/'.analysis_kernels'
km=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(kernel_root)]))
def progress(cell,cell_index,**kwargs):
    print(f'Report cell {cell_index+1}/{len(report.cells)}: '+cell.source.splitlines()[0].encode('ascii','replace').decode(),flush=True)
NotebookClient(report,km=km,timeout=600,resources={'metadata':{'path':str(root)}},on_cell_start=progress).execute()
updated={c.source:c for c in report.cells if c.cell_type=='code'}
count=0
for c in nb.cells:
    if c.cell_type!='code':continue
    if c.source in updated:c.outputs=updated[c.source].outputs
    count+=1;c.execution_count=count
    assert c.outputs or c.source.startswith(('"""','from pathlib'))
nb.cells[0].source+='\n\n*Reporting update: saved solver outputs were reused after verifying algorithm/input hashes. Auxiliary counts and quality metrics were recomputed. Run All executes the full experiment again.*'
nbformat.write(nb,path)
print('Updated notebook reporting and figures.',flush=True)
