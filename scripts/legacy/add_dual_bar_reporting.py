"""Execute only the new plot cells; preserve existing verified notebook outputs."""

# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))

from pathlib import Path
import json,subprocess,sys,hashlib
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
    nb=nbformat.read(path,as_version=4)
    source=next(c for c in nb.cells if c.source.startswith('from cps_quality_plots import academic_dual_bars'))
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
    report=nbformat.v4.new_notebook(cells=[bootstrap,source,nb.cells[-1]])
    km=KernelManager(kernel_name='protein-analysis',kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(root/'.analysis_kernels')]))
    NotebookClient(report,km=km,timeout=600,resources={'metadata':{'path':str(root)}}).execute()
    outputs={c.source:c.outputs for c in old.cells if c.cell_type=='code'}
    outputs.update({c.source:c.outputs for c in report.cells if c.cell_type=='code'})
    count=0
    for c in nb.cells:
        if c.cell_type=='code':
            count+=1;c.outputs=outputs[c.source];c.execution_count=count
    nb.cells[0].source+='\n\n*Reporting update: verified existing solver outputs were retained; the additional grouped-bar figures were newly rendered. Run All recomputes the full experiment.*'
    nbformat.write(nb,path)
    assert sum('image/png' in o.get('data',{}) for o in source.outputs)==len(previous['ws'])*len(previous['alphas'])
    subprocess.run([sys.executable,str(root/'validate_cps_deliverable.py')],check=True,cwd=root)
    print('Added 3 dual-bar figures (PDF, SVG, PNG) to the executed notebook.')
except Exception:
    nbformat.write(old,path)
    raise
