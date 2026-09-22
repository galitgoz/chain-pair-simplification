"""Verify the executed notebook and its experiment exports agree."""
from pathlib import Path
import json,sys
import nbformat
import pandas as pd
import numpy as np

root=Path(__file__).resolve().parents[1]
nb=nbformat.read(root/'notebooks/legacy/0.analyze_CPS.ipynb',as_version=4)
nbformat.validate(nb)
code=[c for c in nb.cells if c.cell_type=='code']
assert all(c.execution_count is not None for c in code)
assert not [o for c in code for o in c.outputs if o.output_type=='error']
out=root/(sys.argv[1] if len(sys.argv)>1 else 'output/cps_papers')
manifest=json.loads((out/'manifest.json').read_text())
r=pd.read_csv(out/'comparison.csv')
inputs=pd.read_csv(out/'inputs.csv')
assert set(r.w)==set(manifest['ws'])
assert set(r.alpha)==set(manifest['alphas'])
parameters=pd.read_csv(out/'pair_parameters.csv')
assert len(parameters)==len(inputs)*len(manifest['alphas'])
assert set(parameters['Stride w'])==set(manifest['ws'])
for side in ('A','B'):
    n=parameters[f'|{side}|'];bound=parameters[f'Theoretical aux upper bound {side}']
    assert (bound==2*n*(n-1)).all()
    assert (parameters[f'Auxiliary points {side}']<=bound).all()
assert r.pair.nunique()==14
assert len(r)==len(inputs)*len(manifest['alphas'])*3
assert manifest['checks']['checks']==100
aux=pd.read_csv(out/'auxiliary_counts.csv')
paired=pd.read_csv(out/'cps2f_vs_cps3f.csv')
assert len(aux)==len(paired)==len(r)//3
assert aux[['auxiliary_A','auxiliary_B']].notna().all().all()
assert (aux.extended_A==aux.nA+aux.auxiliary_A).all()
assert (aux.extended_B==aux.nB+aux.auxiliary_B).all()
assert (aux.auxiliary_total==aux.auxiliary_A+aux.auxiliary_B).all()
assert paired[['pair','delta1','delta2','delta3']].notna().all().all()
for side in ('A','B'):
    valid=r[r[f'k{side}'].notna()]
    assert np.allclose(valid[f'compression_factor_{side}'],valid[f'n{side}']/valid[f'k{side}'])
    assert np.allclose(valid[f'removed_pct_{side}'],100*(1-valid[f'k{side}']/valid[f'n{side}']))
    for metric in ('contact_precision','contact_recall','contact_F1'):
        assert valid[f'{metric}_{side}'].dropna().between(0,1).all()
assert inputs.full_chain_endpoints_preserved.all()
assert inputs.query("tables == '1,2'").dropped_decimated_vertices_per_chain.eq(0).all()
solved=r[r.status.str.startswith('optimal')]
joint=solved[solved.method!='GCS']
assert joint.pair_bound_verified.all()
for method in ('CPS-3F','CPS-2F'):
    rows=joint[joint.method==method]
    for side,delta in [('A','delta1'),('B','delta2')]:
        metric=f'discrete_fidelity_{side}' if method=='CPS-3F' else f'continuous_fidelity_{side}_upper'
        assert (rows[metric]<=rows[delta]+1e-5).all()
pivot=solved.pivot(index='case',columns='method',values='k').dropna()
assert (pivot.GCS<=pivot['CPS-2F']).all()
assert (pivot['CPS-2F']<=pivot['CPS-3F']).all()
paths=json.loads((out/'paths.json').read_text())
for row in solved.itertuples():
    p=paths[f'{row.case}@{row.method}']
    for side,n in [('A',row.nA),('B',row.nB)]:
        ids=p[f'indices_{side}']
        assert ids[0]==0 and ids[-1]==n-1 and np.all(np.diff(ids)>0)
    c=np.array(p['coupling']);steps=np.diff(c,axis=0)
    assert tuple(c[0])==(0,0)
    assert tuple(c[-1])==(len(p['indices_A'])-1,len(p['indices_B'])-1)
    assert np.all((steps>=0)&(steps<=1)) and np.all(steps.sum(axis=1)>0)
figures=list((out/'figures').glob('*.png'))
assert len(figures)>=6
images=sum('image/png' in o.get('data',{}) for c in code for o in c.outputs)
assert images>=6
print(json.dumps(dict(cells=len(nb.cells),executed_code_cells=len(code),rows=len(r),
    pairs=r.pair.nunique(),embedded_figures=images,statuses=r.status.value_counts().to_dict()),indent=2))
