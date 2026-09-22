"""Finalize reporting units and verify the hurricane deliverable independently."""
from pathlib import Path
import json,hashlib,sys
import nbformat
import numpy as np
import pandas as pd
from curve_algorithms import discrete_frechet,continuous_decision
root=Path(__file__).resolve().parents[1];out=root/(sys.argv[1] if len(sys.argv)>1 else 'output/hurricane_cps')
path=root/'notebooks/legacy/1.analyze_CPS.ipynb';nb=nbformat.read(path,4)
# Correct the generic graph HTML's protein-unit label without changing counts.
line="(OUT/'graph_sizes.html').write_text((OUT/'graph_sizes.html').read_text(encoding='utf-8').replace('All distances are in angstroms.','All distances are in projected kilometres.'),encoding='utf-8')"
for c in nb.cells:
    if c.cell_type=='code' and c.source.startswith('from cps_graph_sizes import') and line not in c.source:
        c.source=c.source.replace('graph_sizes=graph_size_table(data,results,OUT)','graph_sizes=graph_size_table(data,results,OUT)\n'+line)
    if c.cell_type=='markdown' and c.source.startswith('## 8.'):
        c.source=c.source.replace('m*m n*n (where m*=|A*| and n*=|B*|)',r'$m^*m\,n^*n$ (where $m^*=|A^*|$ and $n^*=|B^*|$)')
html=out/'graph_sizes.html'
html.write_text(html.read_text(encoding='utf-8').replace('All distances are in angstroms.','All distances are in projected kilometres.'),encoding='utf-8')
nbformat.validate(nb);nbformat.write(nb,path)
code=[c for c in nb.cells if c.cell_type=='code']
assert all(c.execution_count is not None for c in code)
assert not [o for c in code for o in c.outputs if o.output_type=='error']
r=pd.read_csv(out/'comparison.csv');paths=json.loads((out/'paths.json').read_text())
manifest=json.loads((out/'manifest.json').read_text())
manifest['sha256']['src/hurricane_cps.py']=hashlib.sha256((root/'src/hurricane_cps.py').read_bytes()).hexdigest()
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
assert len(r)==manifest['selected_pairs']*len(manifest['alphas'])*3
assert not r.duplicated(['case','method']).any()
for row in r[r.status.str.startswith('optimal')].itertuples():
    a,b=row.pair.split('__');A=np.loadtxt(root/f'data/hurricanes/{a}.txt');B=np.loadtxt(root/f'data/hurricanes/{b}.txt')
    p=paths[f'{row.case}@{row.method}'];ia=np.array(p['indices_A']);ib=np.array(p['indices_B'])
    assert ia[0]==ib[0]==0 and ia[-1]==len(A)-1 and ib[-1]==len(B)-1
    assert (np.diff(ia)>0).all() and (np.diff(ib)>0).all()
    assert (row.nA,row.nB,row.kA,row.kB)==(len(A),len(B),len(ia),len(ib))
    assert continuous_decision(A,A[ia],row.delta1) and continuous_decision(B,B[ib],row.delta2)
    if row.method!='GCS':assert discrete_frechet(A[ia],B[ib])<=row.delta3+1e-9
    assert np.isclose(row.budget_ratio,discrete_frechet(A[ia],B[ib])/row.delta3)
g=pd.read_csv(out/'graph_sizes.csv');u=pd.read_csv(out/'auxiliary_utilization.csv')
assert (g.valid_vertices<=g.cartesian_bound).all() and (g.cartesian_bound<=g.theoretical_vertex_bound).all()
assert (g[g.method.eq('CPS-3F')][['auxiliary_A','auxiliary_B']]==0).all().all()
for curve in ['A','B']:
    assert np.allclose(u[f'utilization_pct_{curve}'],100*u[f'auxiliary_{curve}']/u[f'auxiliary_upper_{curve}'])
report=dict(cells=len(nb.cells),executed_code_cells=len(code),pairs=r.pair.nunique(),rows=len(r),
    statuses=r.status.value_counts().to_dict(),png_figures=len(list((out/'figures').rglob('*.png'))))
(out/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
