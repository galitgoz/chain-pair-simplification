"""Exhaustive independent oracle for both paper configuration solvers."""
from pathlib import Path
import itertools,json,time
import numpy as np
from cps_paper_algorithms import solve_cps,independent_discrete,extended_curve
from curve_algorithms import shortest_independent,continuous_decision,discrete_frechet

def candidates(P,d,kind):
    for k in range(2,len(P)+1):
        for middle in itertools.combinations(range(1,len(P)-1),k-2):
            ids=np.array((0,)+middle+(len(P)-1,))
            ok=continuous_decision(P,P[ids],d) if kind=='CPS-2F' else discrete_frechet(P,P[ids])<=d+1e-9
            if ok:yield ids

def run_checks(trials=24):
    rng=np.random.default_rng(20260914);checks=0;start=time.perf_counter()
    for trial in range(trials):
        A=np.cumsum(rng.normal(size=(rng.integers(3,6),3)),axis=0)
        B=np.cumsum(rng.normal(size=(rng.integers(3,6),3)),axis=0)
        d1,d2=rng.uniform(.1,2.4,2);d3=max(np.linalg.norm(A[0]-B[0]),np.linalg.norm(A[-1]-B[-1]))+rng.uniform(-.1,.8)
        for kind in ('CPS-3F','CPS-2F'):
            ca=list(candidates(A,d1,kind));cb=list(candidates(B,d2,kind))
            expected=min((max(len(a),len(b)) for a in ca for b in cb if discrete_frechet(A[a],B[b])<=d3+1e-9),default=None)
            actual=solve_cps(A,B,d1,d2,d3,kind,certificate=False,seconds=60)
            assert actual.get('k')==expected,(trial,kind,expected,actual)
            independent=shortest_independent(A,d1) if kind=='CPS-2F' else independent_discrete(A,d1)
            assert len(independent)==min(map(len,ca))
            checks+=2
    P=np.array([[0.,0,0],[1,0,0],[2,0,0]])
    assert len(shortest_independent(P,0))==2
    assert len(independent_discrete(P,0))==3;checks+=2
    for kind in ('CPS-2F','CPS-3F'):
        r=solve_cps(P,P,0,0,0,kind,certificate=False)
        assert r['k']==(2 if kind=='CPS-2F' else 3);checks+=1
    report=dict(checks=checks,trials=trials,seconds=time.perf_counter()-start,seed=20260914)
    Path('output/cps_papers').mkdir(parents=True,exist_ok=True)
    Path('output/cps_papers/algorithm_checks.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True);return report

if __name__=='__main__':run_checks()
