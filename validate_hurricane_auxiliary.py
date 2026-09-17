from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import nbformat
root=Path(__file__).resolve().parent;out=root/(sys.argv[1] if len(sys.argv)>1 else 'output/hurricane_auxiliary')
nb=nbformat.read(root/'0.analyze_hurricans.ipynb',4)
nbformat.validate(nb)
code=[c for c in nb.cells if c.cell_type=='code']
assert all(c.execution_count is not None for c in code)
assert not [o for c in code for o in c.outputs if o.output_type=='error']
s=pd.read_csv(out/'auxiliary_construction_summary.csv')
p=pd.read_csv(out/'extended_curve_points.csv')
g=pd.read_csv(out/'graph_sizes.csv');g=g[g.method.eq('CPS-2F')]
for r in s.itertuples():
    pts=p[p.pair.eq(r.pair)&p.alpha.eq(r.alpha)&p.curve.eq(r.curve)]
    assert len(pts)==r.extended_size==r.n_original+r.added_auxiliary
    assert pts.point_type.eq('auxiliary').sum()==r.added_auxiliary
    assert (np.diff(pts.curve_parameter)>0).all()
    name=r.pair.split('__')[0 if r.curve=='A' else 1]
    P=np.loadtxt(root/f'data/hurricanes/{name}.txt')
    t=pts.curve_parameter.to_numpy();i=np.minimum(t.astype(int),len(P)-2)
    expected=P[i]+(t-i)[:,None]*(P[i+1]-P[i])
    assert np.allclose(expected,pts[['x_km','y_km']],rtol=0,atol=1e-7)
    assert np.allclose(pts[pts.point_type.eq('original')][['x_km','y_km']],P)
    prior=g[g.pair.eq(r.pair)&g.alpha.eq(r.alpha)].iloc[0]
    assert r.added_auxiliary==prior[f'auxiliary_{r.curve}']
    assert r.auxiliary_upper==2*len(P)*(len(P)-1)
    assert r.accepted_root_events-r.duplicate_events_removed==r.added_auxiliary
report=dict(cells=len(nb.cells),executed_code_cells=len(code),curve_parameter_cases=len(s),
    original_and_auxiliary_catalog_rows=len(p),verified_against_solver_counts=True,
    construction_figures=len(list((out/'figures').glob('auxiliary_on_curves*.png'))))
(out/'auxiliary_validation.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
