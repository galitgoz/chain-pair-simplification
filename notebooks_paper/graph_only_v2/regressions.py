"""Small-instance optimality and checkpoint lifecycle checks; no pilot diagnostics here."""
from common import ROOT,OUT,solve,validate,sha,atomic,BASE_GRAPH_LIMITS,BUDGETS
from supervisor import attempt
from unittest.mock import patch
import numpy as np,pandas as pd,json,time
import cps_paper_algorithms as alg
from verify_cps_papers import candidates
from curve_algorithms import discrete_frechet
rows=[]
def record(name,fn):
    start=time.perf_counter()
    try:detail=fn();status='PASS';reason=''
    except Exception as error:detail='';status='FAIL';reason=repr(error)
    rows.append(dict(check=name,outcome=status,reason=reason,details=detail,seconds=time.perf_counter()-start))
    pd.DataFrame(rows).to_csv(OUT/'regressions.csv',index=False)
P=np.array([[0.,0.],[1,1],[2,0.]])
fixture=OUT/'regressions/triangle.npz';np.savez(fixture,A=P,B=P)
base=dict(diagnostic='regression',pair='synthetic_triangle',case='synthetic:triangle',alpha=None,input_file=str(fixture.relative_to(ROOT)),
          input_sha256=sha(fixture),method='CPS-2F',delta1=1.,delta2=0.,delta3=0.,graph_limits=BASE_GRAPH_LIMITS,
          phase_budgets=BUDGETS,memory_mib=512,algorithms={x:sha(ROOT/x) for x in ['cps_paper_algorithms.py','curve_algorithms.py']},
          certificate=False,parameter_file_sha256=None)
for dim in [2,3]:
    A=np.pad(P,((0,0),(0,dim-2)));B=A.copy()
    for method in ['Independent continuous','Independent discrete','CPS-2F','CPS-3F']:
        for d1,d2,d3 in [(1.,0.,0.),(0.,0.,0.)]:
            def check(A=A,B=B,m=method,d1=d1,d2=d2,d3=d3):
                c={**base,'method':m,'delta1':d1,'delta2':d2,'delta3':d3}
                kind='CPS-3F' if m in ['Independent discrete','CPS-3F'] else 'CPS-2F'
                ca=list(candidates(A,d1,kind));cb=list(candidates(B,d2,kind))
                expected=max(min(map(len,ca)),min(map(len,cb))) if m.startswith('Independent') else min(max(len(a),len(b)) for a in ca for b in cb if discrete_frechet(A[a],B[b])<=d3+1e-9)
                with patch.object(alg,'alternate_certificate',side_effect=AssertionError('Forbidden certificate search')):
                    raw,trace=solve(A,B,c,certificate=False)
                assert raw['k']==expected
                raw=json.loads(json.dumps(raw,default=lambda x:x.tolist()))
                result=validate(A,B,raw,c);assert result['feasibility_passed']
                if m.startswith('CPS'):assert trace['configuration_graph_entered'] and trace['route']=='configuration_graph'
                return 'Exhaustive optimum agrees; validation and graph route checked.'
            record(f'exhaustive_{dim}d_{method}_{d1}_{d2}',check)

def check_worker(mode):
    c={**base,'diagnostic':'regression_'+mode,'test_mode':mode,
       'phase_budgets':{**BUDGETS,'validation':.4 if mode=='validation_timeout' else 30.}}
    folder=attempt(c,category='regressions')
    assert (folder/'raw_solver.json').exists(),folder
    raw=json.loads((folder/'raw_solver.json').read_text());assert raw['raw_result']['k']==3
    assert raw['external_validation_performed'] is False
    events=[json.loads(s) for s in (folder/'phases.jsonl').read_text().splitlines()]
    saved=next(e for e in events if e.get('raw_checkpoint_saved'))
    validation_start=next(e for e in events if e['phase']=='validation')
    assert saved['monotonic']<validation_start['monotonic']
    assert saved['raw_checkpoint_sha256']==sha(folder/'raw_solver.json')
    validation=json.loads((folder/'validation.json').read_text())
    expected={'normal':'validated','validation_failure':'validation_failed','validation_timeout':'validation_incomplete','optional_failure':'validated'}[mode]
    assert validation['status']==expected,validation
    if mode=='optional_failure':assert json.loads((folder/'measurements.json').read_text())['status']=='optional_measurement_failed'
    return str(folder.relative_to(ROOT))+'; raw checkpoint survived; '+expected
for mode in ['normal','validation_failure','validation_timeout','optional_failure']:
    record('checkpoint_'+mode,lambda mode=mode:check_worker(mode))

def independent_coupling():
    c={**base,'method':'Independent continuous'}
    raw,_=solve(P,P,c,certificate=False);raw=json.loads(json.dumps(raw,default=lambda x:x.tolist()))
    checked=validate(P,P,raw,c)
    assert checked['feasibility_passed'] and not checked['coupling_feasible'] and not checked['coupling_constraint_required']
    return 'Independent baseline coupling exceeds delta3 and remains feasible under its own definition.'
record('independent_delta3_not_required',independent_coupling)
assert all(r['outcome']=='PASS' for r in rows),rows
atomic(OUT/'regression_gate.json',dict(checks=len(rows),passed=len(rows),failed=0))
print('Relevant regression checks passed:',len(rows))
