"""Final notebook/output integrity check (read-only)."""
from pathlib import Path
import json
import nbformat
import pandas as pd
import numpy as np

root=Path(__file__).resolve().parent
nb=nbformat.read(root/'0.analyze_independent.ipynb',as_version=4)
nbformat.validate(nb)
code_cells=[c for c in nb.cells if c.cell_type=='code']
assert all(c.execution_count is not None for c in code_cells)
errors=[o for c in code_cells for o in c.outputs if o.output_type=='error']
assert not errors
images=sum('image/png' in o.get('data',{}) for c in code_cells for o in c.outputs)
assert images==5,images
audit=pd.read_csv(root/'output/identifier_audit.csv')
experiments=pd.read_csv(root/'output/paper_experiments.csv')
results=pd.read_csv(root/'output/independent/results.csv')
assert len(audit)==15 and len(experiments)==21 and len(results)==69
assert audit.fragments.sum()==23
for path in (root/'data/prepared').glob('*.npz'):
    data=np.load(path)
    assert data['xyz'].shape==(len(data['chain_indices']),3)
    assert np.isfinite(data['xyz']).all()
    assert np.all(np.diff(data['chain_indices'])==1)
report=dict(notebook_cells=len(nb.cells),executed_code_cells=len(code_cells),embedded_figures=images,
    chains=len(audit),fragments=int(audit.fragments.sum()),paper_experiments=len(experiments),independent_runs=len(results),errors=len(errors))
print(json.dumps(report,indent=2))
for cell in code_cells:
    if cell.source.startswith('window_A='):
        for out in cell.outputs:
            if 'text/plain' in out.get('data',{}): print(out.data['text/plain'].encode('ascii','replace').decode())
