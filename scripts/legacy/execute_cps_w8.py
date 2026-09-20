"""Execute the user's current notebook at w=8, preserving the w=16 artifacts."""

# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

from pathlib import Path
import re,json,sys,subprocess
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

root=Path(__file__).resolve().parents[2]
path=root/'0.analyze_CPS.ipynb'
backup=root/'output/cps_papers/0.analyze_CPS_before_w8.ipynb'
if not backup.exists():backup.write_bytes(path.read_bytes())
nb=nbformat.read(path,as_version=4)
for c in nb.cells:
    c.source=c.source.replace("'cps_papers'","'cps_papers_w8'").replace('output/cps_papers/','output/cps_papers_w8/')
    c.source=c.source.replace('dual_bars_w16','dual_bars_w8')
    if c.cell_type=='code':
        c.source=re.sub(r'^WS\s*=.*$', 'WS = (8,)  # Run the supplied stride-8 inputs',c.source,flags=re.M)
        c.source=re.sub(r'^ALPHAS\s*=.*$', 'ALPHAS = (1, 2)',c.source,flags=re.M)
        c.outputs=[];c.execution_count=None
    else:
        c.source=c.source.replace('`w=16`','`w=8`').replace('α ∈ {1,2,4}','α ∈ {1,2}').replace('@16','@8')
        c.source=re.sub(r'\n\n\*Reporting update:.*?\*','',c.source,flags=re.S)
        c.source=c.source.replace('**Saved table:** computed directly from the current input configuration; it does not require successful joint-solver runs. Existing solver outputs are preserved.',
                                 '**Saved table:** computed from the same stride-8 inputs and thresholds used for this run.')
nb.cells[0].source+='\n\n**Executed configuration:** w=8, α=1 and 2; results are written to `output/cps_papers_w8/`. Earlier w=16 artifacts remain in `output/cps_papers/`.'
kernel_root=root/'.analysis_kernels';kernel=kernel_root/'protein-analysis'
kernel.mkdir(parents=True,exist_ok=True)
(kernel/'kernel.json').write_text(json.dumps(dict(argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],display_name='Protein analysis',language='python')))
km=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(kernel_root)]))
def progress(cell,cell_index,**kwargs):
    nbformat.write(nb,path)
    print(f'Cell {cell_index+1}/{len(nb.cells)}: {cell.source.splitlines()[0][:85]}'.encode('ascii','replace').decode(),flush=True)
client=NotebookClient(nb,km=km,timeout=7200,resources={'metadata':{'path':str(root)}},on_cell_start=progress)
try:client.execute()
finally:nbformat.write(nb,path)
subprocess.run([sys.executable,str(root/'validate_cps_deliverable.py'),'output/cps_papers_w8'],check=True,cwd=root)
print('Completed and validated the stride-8 notebook.',flush=True)
