from pathlib import Path
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
root=Path(__file__).resolve().parents[1]
path=root/'notebooks/legacy/1.analyze_CPS.ipynb'
nb=nbformat.read(path,4)
km=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(root/'.analysis_kernels')]))
def progress(cell,cell_index,**kwargs):
    nbformat.write(nb,path)
    print(f'Cell {cell_index+1}/{len(nb.cells)}: {cell.source.splitlines()[0][:70]}',flush=True)
try:NotebookClient(nb,km=km,timeout=7200,resources={'metadata':{'path':str(root)}},on_cell_start=progress).execute()
finally:nbformat.write(nb,path)
assert all(c.execution_count is not None and not any(o.output_type=='error' for o in c.outputs) for c in nb.cells if c.cell_type=='code')
print('Hurricane notebook executed successfully.')
