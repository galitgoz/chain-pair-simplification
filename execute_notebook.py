"""Build, validate cached sweep, and execute the notebook using this interpreter."""
from pathlib import Path
import sys,json,hashlib,subprocess,importlib.metadata,time
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

root=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(root/'build_notebook.py')],check=True,cwd=root)
path=root/'0.analyze_independent.ipynb'
nb=nbformat.read(path,as_version=4)
namespace={}
exec(next(c.source for c in nb.cells if c.cell_type=='code'),namespace)
signature_cell=next(c.source for c in nb.cells if c.cell_type=='code' and 'def input_signature()' in c.source)
exec(signature_cell[signature_cell.index('def input_signature():'):signature_cell.index('\nsignature = input_signature()')],namespace)
signature=namespace['input_signature']()
# Validate existing outputs against current preparation before blessing the cache.
from protein_data import prepare
from curve_algorithms import continuous_decision
import pandas as pd
import numpy as np
chains,audit,breaks,missing=prepare(root)
fragments={f['id']:f for c in chains.values() for f in c['fragments']}
results=pd.read_csv(root/'output/independent/results.csv')
paths=json.loads((root/'output/independent/paths.json').read_text())
assert len(results)==len(fragments)*3
for row in results.itertuples():
    f=fragments[row.fragment]; ids=paths[f'{row.fragment}@{row.delta_A:g}']['fragment_indices']
    assert len(ids)==row.k and len(f['xyz'])==row.n
    assert continuous_decision(f['xyz'],f['xyz'][ids],row.delta_A)
(root/'output/independent/manifest.json').write_text(json.dumps(dict(signature=signature,
    mode='executed by benchmark_fragments.py; revalidated against current notebook inputs',
    runs=len(results),thresholds=[2,4,8]),indent=2))
requirements=['numpy','pandas','matplotlib','biopython','numba','nbformat','nbclient','ipykernel','pymupdf']
(root/'requirements-analysis.txt').write_text('\n'.join(f'{name}=={importlib.metadata.version(name)}' for name in requirements)+'\n')
kernel_root=root/'.analysis_kernels'
kernel=kernel_root/'protein-analysis'; kernel.mkdir(parents=True,exist_ok=True)
(kernel/'kernel.json').write_text(json.dumps(dict(argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],
    display_name='Protein analysis',language='python')))
manager=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(kernel_root)]))
def progress(cell,cell_index,**kwargs):
    print(f'Executing cell {cell_index+1}/{len(nb.cells)}: {cell.source.splitlines()[0][:100]}',flush=True)
client=NotebookClient(nb,km=manager,timeout=1800,resources={'metadata':{'path':str(root)}},on_cell_start=progress)
try:
    client.execute()
finally:
    nbformat.write(nb,path)
errors=[(i,o) for i,c in enumerate(nb.cells) if c.cell_type=='code' for o in c.get('outputs',[]) if o.output_type=='error']
assert not errors,errors
print(f'Notebook executed successfully: {len(nb.cells)} cells, {path.stat().st_size} bytes',flush=True)
