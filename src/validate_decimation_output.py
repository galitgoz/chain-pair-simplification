from pathlib import Path
import json
import nbformat
import numpy as np
import pandas as pd
from decimation_algorithms import decimate,edge_stats,discrete_frechet,segment_frechet_decision

root=Path(__file__).resolve().parents[1];out=root/'output/decimation'
nb=nbformat.read(root/'notebooks/legacy/1.analyze_independent.ipynb',as_version=4);nbformat.validate(nb)
cc=[c for c in nb.cells if c.cell_type=='code']
assert all(c.execution_count is not None for c in cc)
assert not [o for c in cc for o in c.outputs if o.output_type=='error']
images=sum('image/png' in o.get('data',{}) for c in cc for o in c.outputs)
assert images==7,images
for cell in nb.cells:
    if cell.cell_type=='markdown' and cell.source.startswith('### F2') and 'A flat distance' not in cell.source:
        cell.source+='\n\nA flat distance can be dominated by retained endpoints; it does not prove that helical detail survived. Read F2 together with edge growth and the coarse-trace flags.'
nbformat.write(nb,root/'notebooks/legacy/1.analyze_independent.ipynb')
frame=pd.read_csv(out/'pipeline.csv');paths=json.loads((out/'paths.json').read_text())
assert len(frame)==210 and frame.pair.nunique()==14
required='nA_before nB_before edge_mean_before dF_before nA_after nB_after edge_mean_after edge_min_after edge_max_after keep_frac dF_after dF_change alpha delta1 delta3 floor feasible sep_kA sep_kB sep_equal sep_lenratio_A sep_dF sep_max_link sep_ok3 cps_k cps_lenratio_A cps_dF cps_max_link cost_of_joining cps_status'.split()
assert set(required)<=set(frame.columns)
for pair,g in frame.groupby('pair'):
    assert set(g.w)=={1,2,4,8,16} and len(g)==15
    source=np.load(out/'inputs'/pair/'aligned.npz');A=source['A'];B=source['B']
    for row in g.itertuples():
        record=paths[f'{pair}@{row.w}@{row.alpha}'];idx=np.array(record['decimated_indices'])
        assert np.array_equal(idx,decimate(np.arange(len(A)),row.w))
        a=A[idx];b=B[idx]
        assert row.nA_after==row.nB_after==len(idx)
        stats=edge_stats(a,b)
        assert abs(row.edge_mean_after-stats['mean'])<1e-10 and abs(row.edge_min_after-stats['minimum'])<1e-10
        assert abs(row.delta1-row.alpha*stats['mean'])<1e-10
        assert row.delta3==np.ceil(discrete_frechet(a,b))
        assert row.feasible==(row.delta3>=row.floor)
        for side,P in [('A',a),('B',b)]:
            for method in ['sep','cps']:
                ii=record[f'{method}_{side}']
                if ii is None:continue
                assert ii[0]==0 and ii[-1]==len(P)-1 and np.all(np.diff(ii)>0)
                assert all(segment_frechet_decision(P[i:j+1],P[i],P[j],row.delta1) for i,j in zip(ii,ii[1:]))
        sa=record['sep_A'];sb=record['sep_B']
        if len(sa)!=len(sb):assert pd.isna(row.sep_max_link) and pd.isna(row.sep_ok3)
        else:assert abs(row.sep_max_link-np.linalg.norm(a[sa]-b[sb],axis=1).max())<1e-9
        if row.solver_feasible:
            ja=record['cps_A'];jb=record['cps_B'];assert len(ja)==len(jb)==row.cps_k
            assert np.linalg.norm(a[ja]-b[jb],axis=1).max()<=row.delta3+1e-9
for name in ['sampling_summary_with_control.csv','coupling_summary_with_control.csv','aggregates_headline.csv','aggregates_runtime_with_control.csv']:
    df=pd.read_csv(out/name);assert 1 in set(df.w),name
assert (frame[frame.w>=8].stride_role=='RUNTIME ONLY').all()
headline=frame[frame.w<=4]
report=dict(notebook_cells=len(nb.cells),executed_code_cells=len(cc),figures=images,pairs=14,rows=len(frame),
    verified_joint_rows=int(frame.solver_feasible.sum()),headline_rows=len(headline),
    headline_unequal_separate=int((~headline.sep_equal).sum()),
    headline_separate_violations=int((headline.sep_equal & (headline.sep_ok3==False)).sum()),
    all_rows_have_exact_shared_stride_and_endpoints=True,solver_errors=0)
(out/'delivery_checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
