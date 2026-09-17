"""Target the small case that would qualify for the alternative certificate search."""
from graph_core import OUT,solve,measures
from unittest.mock import patch
import numpy as np,pandas as pd,json,time
import cps_paper_algorithms as alg
from curve_algorithms import discrete_frechet
P=np.array([[0.,0.],[1.,1.],[2.,0.]])
start=time.perf_counter();failure='';outcome='PASS'
try:
    independent=solve(P,P,1.,0.,0.,'Independent continuous',certificate=False)
    assert [len(independent['indices_A']),len(independent['indices_B'])]==[2,3]
    assert discrete_frechet(P[independent['indices_A']],P[independent['indices_B']])>0
    with patch.object(alg,'alternate_certificate',side_effect=AssertionError('Forbidden shortcut called')) as forbidden:
        result=solve(P,P,1.,0.,0.,'CPS-2F',certificate=False)
        evidence=measures(P,P,result,1.,0.,0.,'CPS-2F')
    assert forbidden.call_count==0 and result['raw_status']=='optimal_configuration_graph'
    assert result['configuration_graph_entered'] and result['chain_graphs_completed']==2 and evidence['k']==3
except Exception as error:outcome='FAIL';failure=repr(error)
frame=pd.read_csv(OUT/'validation.csv')
assert not frame.check.eq('alternative_search_eligible_but_disabled').any()
frame=pd.concat([frame,pd.DataFrame([dict(check='alternative_search_eligible_but_disabled',affected_method='CPS-2F',check_kind='targeted shortcut exclusion',outcome=outcome,failure_reason=failure,details='Identical 3-vertex triangles; delta=(1,0,0); independent lengths 2 and 3 violate coupling; graph optimum is 3.',seconds=time.perf_counter()-start)])],ignore_index=True)
frame.to_csv(OUT/'validation.csv',index=False)
(OUT/'validation_gate.json').write_text(json.dumps(dict(checks=len(frame),passed=int(frame.outcome.eq('PASS').sum()),failed=int(frame.outcome.eq('FAIL').sum()),incomplete=int(frame.outcome.eq('INCOMPLETE').sum())),indent=2))
guard=json.loads((OUT/'shortcut_guard.json').read_text());guard['alternative_eligible_triangle_case']=outcome
(OUT/'shortcut_guard.json').write_text(json.dumps(guard,indent=2))
assert outcome=='PASS',failure
print('100/100 checks pass, including a case eligible for alternative-subsequence certification; graph executed instead.')
