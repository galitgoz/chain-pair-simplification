"""Local independent simplification and exact synchronous discrete CPS.

The joint method requires equal-cardinality subsequences with indexwise pairing.
This is a stricter endpoint-preserving variant of Fan et al.'s CPS-3F.
"""
import numpy as np
from numba import njit
from curve_algorithms import segment_frechet_decision,discrete_frechet,local_shortcut_path

EPS=1e-9

def arc_length(P): return float(np.linalg.norm(np.diff(P,axis=0),axis=1).sum())

@njit(cache=True)
def _synchronous_dp(A,B,d1,d2,d3,minimax):
    """Four-walk DAG; prefix minima compress simultaneous retained-vertex jumps.

    i,j: consumed original vertices; p,q: current retained vertices.
    Original indices advance by 0/1 independently. Retained indices either both
    stay or both strictly advance, costing one new paired vertex. Start/end
    retained vertices are fixed. Prefix tables store value and argmin over p,q.
    In minimax mode minimize the largest paired distance (exact feasibility floor).
    """
    n=len(A); m=len(B)
    costs=np.full((n,m,n,m),np.inf)
    prefix=np.full_like(costs,np.inf)
    arg=np.full((n,m,n,m),-1,np.int32)
    predecessor=np.full((n,m,n,m),-1,np.int32)
    for i in range(n):
        for j in range(m):
            for p in range(n):
                for q in range(m):
                    pair_distance=np.linalg.norm(A[p]-B[q])
                    valid=(np.linalg.norm(A[i]-A[p])<=d1+EPS and np.linalg.norm(B[j]-B[q])<=d2+EPS
                           and (minimax or pair_distance<=d3+EPS))
                    value=np.inf; parent=-1
                    if valid:
                        if i==0 and j==0 and p==0 and q==0:
                            value=pair_distance if minimax else 1.
                        for di in range(2):
                            for dj in range(2):
                                ii=i-di; jj=j-dj
                                if ii<0 or jj<0: continue
                                if di+dj>0:
                                    v=costs[ii,jj,p,q]
                                    if minimax: v=max(v,pair_distance)
                                    if v<value:
                                        value=v; parent=((ii*m+jj)*n+p)*m+q
                                if p>0 and q>0:
                                    v=prefix[ii,jj,p-1,q-1]
                                    if minimax: v=max(v,pair_distance)
                                    else: v+=1.
                                    if v<value: value=v; parent=arg[ii,jj,p-1,q-1]
                    costs[i,j,p,q]=value; predecessor[i,j,p,q]=parent
                    best=value; bestarg=((i*m+j)*n+p)*m+q if np.isfinite(value) else -1
                    if p>0 and prefix[i,j,p-1,q]<best:
                        best=prefix[i,j,p-1,q]; bestarg=arg[i,j,p-1,q]
                    if q>0 and prefix[i,j,p,q-1]<best:
                        best=prefix[i,j,p,q-1]; bestarg=arg[i,j,p,q-1]
                    prefix[i,j,p,q]=best; arg[i,j,p,q]=bestarg
    value=costs[n-1,m-1,n-1,m-1]
    if not np.isfinite(value): return value,np.empty(0,np.int64),np.empty(0,np.int64)
    i=n-1;j=m-1;p=n-1;q=m-1
    aa=[p];bb=[q]
    while predecessor[i,j,p,q]>=0:
        state=predecessor[i,j,p,q]
        qq=state%m; state//=m; pp=state%n; state//=n; jj=state%m; ii=state//m
        if pp!=p or qq!=q: aa.append(pp);bb.append(qq)
        i=ii;j=jj;p=pp;q=qq
    return value,np.array(aa[::-1]),np.array(bb[::-1])

def synchronous_cps(A,B,delta1,delta2,delta3):
    A=np.ascontiguousarray(A,float); B=np.ascontiguousarray(B,float)
    if len(A)<2 or len(B)<2 or min(delta1,delta2,delta3)<0: raise ValueError('Invalid curves or thresholds')
    floor,fa,fb=_synchronous_dp(A,B,delta1,delta2,np.inf,True)
    if delta3+EPS<floor:
        return dict(status='infeasible: delta3 below exact feasibility floor',floor=float(floor),indices_A=None,indices_B=None)
    value,ia,ib=_synchronous_dp(A,B,delta1,delta2,delta3,False)
    if not np.isfinite(value): raise AssertionError('Solver failed despite feasible minimax certificate')
    assert len(ia)==len(ib)==int(value)
    assert ia[0]==ib[0]==0 and ia[-1]==len(A)-1 and ib[-1]==len(B)-1
    assert discrete_frechet(A,A[ia])<=delta1+EPS and discrete_frechet(B,B[ib])<=delta2+EPS
    assert np.linalg.norm(A[ia]-B[ib],axis=1).max()<=delta3+EPS
    assert discrete_frechet(A[ia],B[ib])<=delta3+EPS
    return dict(status='feasible',floor=float(floor),indices_A=ia,indices_B=ib,k=int(value))

def cps_2f(*args,**kwargs):
    raise NotImplementedError('CPS-2F definition and further data are reserved for a later notebook.')
