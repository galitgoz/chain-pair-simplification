"""Add the explicit per-pair parameter table without repeating solver work."""

# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

from pathlib import Path
import json,subprocess,sys,hashlib,copy
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

root=Path(__file__).resolve().parents[2];path=root/'0.analyze_CPS.ipynb'
old=nbformat.read(path,as_version=4)
previous=json.loads((root/'output/cps_papers/manifest.json').read_text())
for name in ('curve_algorithms.py','cps_paper_algorithms.py','cps_notebook_run.py','cps_notebook_data.py'):
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==previous['sha256'][name]
try:
    subprocess.run([sys.executable,str(root/'build_cps_notebook.py')],check=True,cwd=root)
    generated=nbformat.read(path,as_version=4)
    index=next(i for i,c in enumerate(generated.cells) if c.source.startswith('from cps_parameter_table import configured_parameter_table'))
    source=generated.cells[index]
    nb=copy.deepcopy(old)
    existing=next((i for i,c in enumerate(nb.cells) if c.source.startswith('### 6.2 Per-pair parameter table')),None)
    if existing is None:
        insertion=next(i for i,c in enumerate(nb.cells) if c.source.startswith('## 7.'))
        nb.cells[insertion:insertion]=[generated.cells[index-1],source]
    else:
        insertion=existing
        nb.cells[insertion:insertion+2]=[generated.cells[index-1],source]
    bootstrap=nbformat.v4.new_code_cell('''from pathlib import Path
import json, hashlib
import pandas as pd
from IPython.display import display, HTML, Markdown
ROOT=Path.cwd(); OUT=ROOT/'output/cps_papers'
results=pd.read_csv(OUT/'comparison.csv')
previous=json.loads((OUT/'manifest.json').read_text())
WS=tuple(previous['ws']); ALPHAS=tuple(previous['alphas'])
INCLUDE_TABLE3=previous['include_table3']; INPUT_SCOPE=previous['input_scope']
BUDGETS=previous['budgets']; checks=previous['checks']
CONTACT_CUTOFF=previous['contact_cutoff']
CONTACT_SEQUENCE_SEPARATION=previous['contact_sequence_separation']
''')
    # Use the notebook's current input settings, rather than the old result CSV.
    loading=copy.deepcopy(next(c for c in old.cells if c.cell_type=='code' and 'data = load_inputs(' in c.source))
    report=nbformat.v4.new_notebook(cells=[loading,source])
    km=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(root/'.analysis_kernels')]))
    NotebookClient(report,km=km,timeout=600,resources={'metadata':{'path':str(root)}}).execute()
    outputs={c.source:c.outputs for c in old.cells if c.cell_type=='code'}
    outputs[source.source]=report.cells[-1].outputs
    count=0
    for c in nb.cells:
        if c.cell_type=='code':
            count+=1;c.outputs=outputs[c.source];c.execution_count=count
    nb.cells[insertion].source+='\n\n**Saved table:** computed directly from the current input configuration; it does not require successful joint-solver runs. Existing solver outputs are preserved.'
    nbformat.write(nb,path)
    assert any('text/html' in o.get('data',{}) for o in source.outputs)
    nbformat.validate(nb)
    assert not [o for c in report.cells for o in c.get('outputs',[]) if o.output_type=='error']
    print('Added per-pair parameter table and auxiliary counts to the notebook, HTML and CSV.')
    import csv
    with (root/'output/cps_papers/pair_parameters.csv').open(encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
        assert len(rows)>0 and len({r['Protein pair (A / B)'] for r in rows})==14
        for row in rows:
            for side in ('A','B'):
                n=int(row[f'|{side}|']);bound=int(row[f'Theoretical aux upper bound {side}'])
                assert bound==2*n*(n-1)
                assert int(row[f'Auxiliary points {side}'])<=bound
        print(f'Checked theoretical bounds and measured counts for {len(rows)} configurations.')
except Exception:
    nbformat.write(old,path)
    raise
