from pathlib import Path
import itertools,json,time
import numpy as np
from decimation_algorithms import decimate,lockstep_decimate,shortcut_dag,shortest_path,strict_pairing_floor,product_dag,joint_local

def paths(E):
    n=len(E)
    return [np.array((0,)+middle+(n-1,)) for k in range(n-1)
        for middle in itertools.combinations(range(1,n-1),k)
        if all(E[i,j] for i,j in zip((0,)+middle,middle+(n-1,)))]

def run_checks():
    start=time.perf_counter();checks=0;rng=np.random.default_rng(20260914)
    for n in [2,3,17,30,33,325]:
        A=rng.normal(size=(n,3));B=rng.normal(size=(n,3))
        assert decimate(A,1) is A;checks+=1
        for w in (1,2,4,8,16):
            expected=list(range(0,n,w))
            if expected[-1]!=n-1:expected.append(n-1)
            aa,bb,ids=lockstep_decimate(A,B,w)
            assert list(ids)==expected and np.array_equal(aa,A[expected]) and np.array_equal(bb,B[expected]);checks+=1
    for trial in range(80):
        A=np.cumsum(rng.normal(size=(int(rng.integers(3,8)),3)),axis=0)
        B=np.cumsum(rng.normal(size=(int(rng.integers(3,8)),3)),axis=0)
        d=float(rng.uniform(.1,2.5));EA=shortcut_dag(A,d);EB=shortcut_dag(B,d)
        pa=paths(EA);pb=paths(EB)
        assert len(shortest_path(EA))==min(map(len,pa));checks+=1
        floor=float(strict_pairing_floor(A,B))
        assert abs(floor-max(np.linalg.norm(A[0]-B[0]),np.linalg.norm(A[-1]-B[-1])))<1e-8;checks+=1
        candidates=[(float(np.linalg.norm(A[a]-B[b],axis=1).max()),len(a)) for a in pa for b in pb if len(a)==len(b)]
        exact_floor=min((d for d,k in candidates),default=np.inf)
        computed,_,_=product_dag(A,B,EA,EB,np.inf,True)
        assert (np.isinf(exact_floor) and np.isinf(computed)) or abs(exact_floor-computed)<1e-8;checks+=1
        for threshold in [floor,exact_floor if np.isfinite(exact_floor) else floor+10]:
            expected=min((k for d,k in candidates if d<=threshold+1e-9),default=None)
            result=joint_local(A,B,EA,EB,threshold)
            if expected is None:assert result['ia'] is None and result['status'].startswith('infeasible')
            else:assert result['k']==expected
            checks+=1
    A=np.array([[0.,0,0],[0,100,0],[0,0,0]])
    B=np.array([[0.,0,0],[100,0,0],[0,0,0]])
    result=joint_local(A,B,shortcut_dag(A,0),shortcut_dag(B,0),0)
    assert result['floor']==0 and result['status']=='infeasible_local_fidelity';checks+=1
    result=joint_local(A,B+1,shortcut_dag(A,0),shortcut_dag(B,0),0)
    assert result['status']=='infeasible_endpoint_floor';checks+=1
    report=dict(checks=checks,random_instances=80,seed=20260914,seconds=time.perf_counter()-start,
        coverage='exact stride/lockstep/endpoints; local shortest path; strict pairing floor; joint optimum and both infeasibility statuses')
    out=Path('output/decimation');out.mkdir(parents=True,exist_ok=True)
    (out/'checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
    return report

if __name__=='__main__':run_checks()
