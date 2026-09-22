"""Refresh track figures from saved solutions without repeating optimization."""

# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path[:0] = [str(_Path(__file__).resolve().parents[2] / 'src'), str(_Path(__file__).resolve().parents[2])]

from pathlib import Path
import copy
import nbformat as nbf
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
root=Path(__file__).resolve().parents[2];path=root/'notebooks/legacy/1.analyze_CPS.ipynb'
nb=nbf.read(path,4)
loading=next(c for c in nb.cells if c.cell_type=='code' and 'tracks, metadata, candidate_pairs=load_hurricanes' in c.source)
index=next(i for i,c in enumerate(nb.cells) if c.cell_type=='code' and c.source.startswith('track_plots('))
bootstrap=nbf.v4.new_code_cell("data,selected=organize(tracks,candidate_pairs,PAIRS_PER_BASIN)\nresults=pd.read_csv(OUT/'comparison.csv')\npaths=json.loads((OUT/'paths.json').read_text())")
mini=nbf.v4.new_notebook(cells=[copy.deepcopy(loading),bootstrap,nbf.v4.new_code_cell(nb.cells[index].source)])
km=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(root/'.analysis_kernels')]))
NotebookClient(mini,km=km,timeout=600,resources={'metadata':{'path':str(root)}}).execute()
nb.cells[index].outputs=mini.cells[-1].outputs
nbf.write(nb,path)
print('Refreshed eight track figures from saved solutions.')
