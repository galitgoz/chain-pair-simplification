"""Independent small-instance checks and a protein-scale timing probe."""
import itertools, json, time
from pathlib import Path
import numpy as np
from curve_algorithms import (shortest_independent,continuous_decision,continuous_distance,
    discrete_frechet,segment_frechet_decision,local_shortcut_path,cps_discrete_reference)

def brute(P,d):
    for k in range(2,len(P)+1):
        for middle in itertools.combinations(range(1,len(P)-1),k-2):
            ids=(0,)+middle+(len(P)-1,)
            if continuous_decision(P,P[list(ids)],d): return np.array(ids)

def run_checks():
    start=time.perf_counter(); checks=0
    line=np.array([[0.,0,0],[1,0,0],[2,0,0]])
    assert list(shortest_independent(line,0))==[0,2]; checks+=1
    assert abs(discrete_frechet(line,line[[0,2]])-1)<1e-9; checks+=1
    assert continuous_decision(line,line[[0,2]],0); checks+=1
    repeated=np.array([[0.,0,0],[0,0,0],[1,0,0],[1,0,0]])
    assert len(shortest_independent(repeated,0))==2; checks+=1
    tangent=np.array([[0.,0,0],[1,1,0],[2,0,0]])
    assert continuous_decision(tangent,tangent[[0,2]],1); checks+=1
    assert not continuous_decision(tangent,tangent[[0,2]],0.99); checks+=1
    backtrack=np.array([[0.,0,0],[2,0,0],[1,0,0],[3,0,0]])
    assert not segment_frechet_decision(backtrack,backtrack[0],backtrack[-1],0.49); checks+=1
    assert segment_frechet_decision(backtrack,backtrack[0],backtrack[-1],0.5); checks+=1
    rng=np.random.default_rng(749)
    counterexample=None
    for trial in range(160):
        P=np.cumsum(rng.normal(size=(int(rng.integers(3,9)),3)),axis=0)
        d=float(rng.uniform(.1,2.5))
        actual=shortest_independent(P,d); oracle=brute(P,d)
        assert len(actual)==len(oracle),(P,d,actual,oracle)
        assert actual[0]==0 and actual[-1]==len(P)-1 and np.all(np.diff(actual)>0)
        assert continuous_decision(P,P[actual],d)
        assert continuous_decision(P,P[actual],discrete_frechet(P,P[actual])+1e-8)
        local=local_shortcut_path(P,d)
        if len(actual)<len(local) and counterexample is None:
            counterexample=dict(P=P.tolist(),delta=d,global_ids=actual.tolist(),local_ids=local.tolist())
        # Segment-specific oracle uses a different 1D monotone-interval argument.
        assert continuous_decision(P,P[[0,-1]],d)==segment_frechet_decision(P,P[0],P[-1],d)
        checks+=5
    result=cps_discrete_reference(line,line,1,1,0)
    assert result['k']==1; checks+=1
    result=cps_discrete_reference(line,line,1,1,0,fix_endpoints=True)
    assert result['k']==2; checks+=1
    assert cps_discrete_reference(line,line+100,1,1,0) is None; checks+=1
    report=dict(checks=checks,random_seed=749,random_curves=160,elapsed_seconds=time.perf_counter()-start,
                local_global_counterexample=counterexample)
    Path('output').mkdir(exist_ok=True)
    Path('output/algorithm_checks.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report

if __name__=='__main__': run_checks()
