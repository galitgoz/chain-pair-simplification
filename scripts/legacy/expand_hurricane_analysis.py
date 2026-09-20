"""Expand existing notebooks without rebuilding or discarding their custom cells."""

# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

from pathlib import Path
import json,re,shutil,copy,sys
import nbformat as nbf
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

ROOT=Path(__file__).resolve().parents[2]
per_basin=None if '--all' in sys.argv else 16
suffix='all' if per_basin is None else 'expanded'
base_out=f'output/hurricane_cps_{suffix}'
aux_out=f'output/hurricane_auxiliary_{suffix}'
backup=ROOT/'output/hurricane_before_expansion';backup.mkdir(exist_ok=True)
for filename in ['1.analyze_CPS.ipynb','0.analyze_hurricans.ipynb']:
    if not (backup/filename).exists():shutil.copy2(ROOT/filename,backup/filename)
target=ROOT/base_out;target.mkdir(exist_ok=True)
for filename in ['comparison.csv','paths.json','manifest.json']:
    if not (target/filename).exists():shutil.copy2(ROOT/'output/hurricane_cps'/filename,target/filename)

def configure(nb,old_out,new_out):
    nb=nbf.reads(re.sub(re.escape(old_out)+r'(?:_expanded|_all)?',new_out,nbf.writes(nb)),4)
    for c in nb.cells:
        if c.cell_type=='code':
            c.source=re.sub(r'^PAIRS_PER_BASIN=.*$',f'PAIRS_PER_BASIN={per_basin!r}  # None selects all supplied pairs',c.source,flags=re.M)
            c.source=re.sub(r'^BUDGETS=.*$',"BUDGETS=dict(max_states=300_000,max_transitions=20_000_000,seconds=120.,max_chain_edges=2_000_000)",c.source,flags=re.M)
            c.outputs=[];c.execution_count=None
        else:
            c.source=c.source.replace('four pairs per basin','16 pairs per basin' if per_basin else 'all available pairs in each basin')
            c.source=c.source.replace('(eight pairs total)','(32 pairs total)' if per_basin else '(2,222 pairs total)')
            c.source=c.source.replace('This is a selected comparison, not a census of all 2,222 candidate pairs.',
                                    'This is a selected comparison, not a census of all 2,222 candidate pairs.' if per_basin else 'This run covers all 2,222 supplied screened pairs.')
            c.source=c.source.replace('eight selected pairs','32 selected pairs' if per_basin else 'all supplied pairs')
            c.source=c.source.replace('The three-method results are reused from the verified hurricane run on these identical inputs and settings (45 optimal results, three CPS-2F resource limits). Running all cells recomputes them.',
                'The three-method results are synchronized with the expanded verified hurricane run on identical inputs and thresholds. Verified optimal results are reused; unfinished cases are retried with the configured limits.')
    nb.cells[0].source+='\n\n**Expanded execution:** '+('32 selected pairs' if per_basin else 'all 2,222 supplied pairs')+'; α=1 and 2; 120-second soft time limit and 20-million-transition limit per joint solve. The 300,000-state and two-million single-chain-edge limits remain. Saved optimal results are reused only after input and algorithm hash checks; unfinished cases are retried. Earlier results remain in their original output folders.'
    return nb

def client(nb):
    km=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(ROOT/'.analysis_kernels')]))
    return NotebookClient(nb,km=km,timeout=172800,resources={'metadata':{'path':str(ROOT)}})

path=ROOT/'1.analyze_CPS.ipynb'
nb=configure(nbf.read(path,4),'output/hurricane_cps',base_out)
def progress(cell,cell_index,**kwargs):
    nbf.write(nb,path);print(f'Base cell {cell_index+1}/{len(nb.cells)}',flush=True)
c=client(nb);c.on_cell_start=progress
try:c.execute()
finally:nbf.write(nb,path)

# Preserve the auxiliary notebook's own instructional cells. Shared cells receive
# matching refreshed results; only auxiliary construction/reporting runs again.
path_aux=ROOT/'0.analyze_hurricans.ipynb'
aux=configure(nbf.read(path_aux,4),'output/hurricane_auxiliary',aux_out)
def normalized(source):
    return source.replace(base_out,'OUTPUT').replace(aux_out,'OUTPUT')
base_code={normalized(c.source):c for c in nb.cells if c.cell_type=='code'}
new_indices=[]
for i,c in enumerate(aux.cells):
    if c.cell_type!='code':continue
    match=base_code.get(normalized(c.source))
    if match is not None:
        c.outputs=copy.deepcopy(match.outputs);c.execution_count=match.execution_count
    else:new_indices.append(i)
shutil.copytree(ROOT/base_out,ROOT/aux_out,dirs_exist_ok=True)
loading=next(c for c in aux.cells if c.cell_type=='code' and 'tracks, metadata, candidate_pairs=load_hurricanes' in c.source)
bootstrap=nbf.v4.new_code_cell("data,selected=organize(tracks,candidate_pairs,PAIRS_PER_BASIN)\nparameters=pd.read_csv(OUT/'parameters.csv')\nresults=pd.read_csv(OUT/'comparison.csv')\npaths=json.loads((OUT/'paths.json').read_text())")
mini=nbf.v4.new_notebook(cells=[copy.deepcopy(loading),bootstrap]+[copy.deepcopy(aux.cells[i]) for i in new_indices])
c=client(mini)
c.on_cell_start=lambda cell,cell_index,**kwargs:print(f'Auxiliary cell {cell_index+1}/{len(mini.cells)}',flush=True)
c.execute()
count=max(c.execution_count or 0 for c in aux.cells if c.cell_type=='code')
for i,executed in zip(new_indices,mini.cells[2:]):
    count+=1;aux.cells[i].outputs=executed.outputs;aux.cells[i].execution_count=count
nbf.write(aux,path_aux)
for delivered in [nb,aux]:
    nbf.validate(delivered)
    assert all(c.execution_count is not None and not any(o.output_type=='error' for o in c.outputs) for c in delivered.cells if c.cell_type=='code')
print('Both notebooks expanded and executed successfully.',flush=True)
