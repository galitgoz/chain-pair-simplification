from pathlib import Path
import sys,json,subprocess
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
root=Path(__file__).resolve().parents[1]
subprocess.run([sys.executable,str(root/'src/build_decimation_notebook.py')],check=True,cwd=root)
path=root/'notebooks/legacy/1.analyze_independent.ipynb';nb=nbformat.read(path,as_version=4)
kernel_root=root/'.analysis_kernels';kernel=kernel_root/'protein-analysis';kernel.mkdir(parents=True,exist_ok=True)
(kernel/'kernel.json').write_text(json.dumps(dict(argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],display_name='Protein analysis',language='python')))
km=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(kernel_root)]))
def progress(cell,cell_index,**kwargs):print(f'Cell {cell_index+1}/{len(nb.cells)}: {cell.source.splitlines()[0][:80]}'.encode('ascii','replace').decode(),flush=True)
client=NotebookClient(nb,km=km,timeout=1800,resources={'metadata':{'path':str(root)}},on_cell_start=progress)
try:client.execute()
finally:nbformat.write(nb,path)
assert not [o for c in nb.cells if c.cell_type=='code' for o in c.outputs if o.output_type=='error']
print(f'Executed notebook successfully: {path.stat().st_size} bytes',flush=True)
