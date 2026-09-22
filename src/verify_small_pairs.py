from pathlib import Path
import itertools,json,time
import numpy as np
from small_pair_solvers import synchronous_cps,_synchronous_dp
from curve_algorithms import discrete_frechet,local_shortcut_path,segment_frechet_decision

def exhaustive(A,B,d):
    def candidates(P):
        return [np.array((0,)+middle+(len(P)-1,)) for k in range(len(P)-1)
            for middle in itertools.combinations(range(1,len(P)-1),k)
            if discrete_frechet(P,P[list((0,)+middle+(len(P)-1,))])<=d+1e-9]
    feasible=[]
    for ia in candidates(A):
        for ib in candidates(B):
            if len(ia)==len(ib): feasible.append((float(np.linalg.norm(A[ia]-B[ib],axis=1).max()),len(ia)))
    return feasible

def run_tests():
    rng=np.random.default_rng(711); tests=0; start=time.perf_counter()
    for case in range(60):
        A=np.cumsum(rng.normal(size=(int(rng.integers(3,8)),3)),axis=0)
        B=np.cumsum(rng.normal(size=(int(rng.integers(3,8)),3)),axis=0)
        delta=float(rng.uniform(.1,3.))
        oracle=exhaustive(A,B,delta)
        floor=min((x for x,k in oracle),default=np.inf)
        actual,_,_=_synchronous_dp(A,B,delta,delta,np.inf,True)
        assert (np.isinf(floor) and np.isinf(actual)) or abs(actual-floor)<1e-8
        tests+=1
        for threshold in [float(discrete_frechet(A,B)),floor if np.isfinite(floor) else 100.]:
            solution=synchronous_cps(A,B,delta,delta,threshold)
            optimum=min((k for d,k in oracle if d<=threshold+1e-9),default=None)
            if optimum is None: assert solution['indices_A'] is None
            else: assert solution['k']==optimum
            tests+=1
        ids=local_shortcut_path(A,delta)
        possible=[]
        for k in range(len(A)-1):
            for middle in itertools.combinations(range(1,len(A)-1),k):
                path=(0,)+middle+(len(A)-1,)
                if all(segment_frechet_decision(A[i:j+1],A[i],A[j],delta) for i,j in zip(path,path[1:])):
                    possible.append(len(path))
        assert len(ids)==min(possible); tests+=1
    # A threshold can meet endpoint bounds yet fail the synchronous feasibility floor.
    A=np.array([[0.,0,0],[1,0,0],[2,0,0],[3,0,0]])
    B=np.array([[0.,0,0],[.1,0,0],[1,0,0],[3,0,0]])
    sol=synchronous_cps(A,B,0,0,0.9)
    assert sol['indices_A'] is None and sol['floor']>0.9; tests+=1
    report=dict(checks=tests,random_seed=711,random_cases=60,seconds=time.perf_counter()-start,
        tested='independent local optimality; synchronous CPS optimality; exact minimax floor; infeasibility')
    out=Path('output/small_pairs');out.mkdir(parents=True,exist_ok=True)
    (out/'checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
    return report

if __name__=='__main__': run_tests()
