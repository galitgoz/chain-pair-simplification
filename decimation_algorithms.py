"""Exact user-specified stride decimation and local continuous shortcut solvers."""
import numpy as np
from numba import njit
from curve_algorithms import segment_frechet_decision,discrete_frechet

def decimate(P, w):
    if w == 1:
        return P
    idx = list(range(0, len(P), w))
    if idx[-1] != len(P) - 1:
        idx.append(len(P) - 1)
    return P[np.array(idx)]

def lockstep_decimate(A,B,w):
    """Choose indices ONCE, by the supplied decimate function, then use both sides."""
    if len(A)!=len(B): raise ValueError('Lockstep input chains must have the same length')
    if len(A)<2 or w not in (1,2,4,8,16): raise ValueError('Need >=2 points and a specified stride')
    indices=decimate(np.arange(len(A)),w)
    aa=A[indices];bb=B[indices]
    assert indices[0]==0 and indices[-1]==len(A)-1
    assert np.array_equal(aa,decimate(A,w)) and np.array_equal(bb,decimate(B,w))
    return aa,bb,indices

@njit(cache=True)
def shortcut_dag(P,delta):
    n=len(P);edges=np.zeros((n,n),np.bool_)
    for i in range(n):
        for j in range(i+1,n):
            edges[i,j]=segment_frechet_decision(P[i:j+1],P[i],P[j],delta)
    return edges

@njit(cache=True)
def shortest_path(edges):
    n=len(edges);cost=np.full(n,n+1,np.int32);parent=np.full(n,-1,np.int32);cost[0]=1
    for j in range(1,n):
        for i in range(j):
            if edges[i,j] and cost[i]+1<cost[j]:cost[j]=cost[i]+1;parent[j]=i
    if parent[-1]<0: raise ValueError('No independent path')
    path=[n-1]
    while path[-1]>0:path.append(parent[path[-1]])
    return np.array(path[::-1])

@njit(cache=True)
def strict_pairing_floor(A,B):
    """Bottleneck over all strictly increasing pairings, without fidelity constraints."""
    n=len(A);m=len(B);value=np.full((n,m),np.inf);prefix=np.full((n,m),np.inf)
    for i in range(n):
        for j in range(m):
            d=np.linalg.norm(A[i]-B[j])
            if i==0 and j==0:value[i,j]=d
            elif i>0 and j>0:value[i,j]=max(d,prefix[i-1,j-1])
            best=value[i,j]
            if i>0:best=min(best,prefix[i-1,j])
            if j>0:best=min(best,prefix[i,j-1])
            prefix[i,j]=best
    return value[-1,-1]

@njit(cache=True)
def product_dag(A,B,EA,EB,delta3,minimax=False):
    """Shortest paired path, or fidelity-aware bottleneck precheck, on two DAGs.

    Vertices are retained pairs (i,j); edges require i and j both strictly
    increasing and each local continuous shortcut admissible. Endpoints fixed.
    Memory O(n*m+n^2+m^2); time O(|EA|*|EB|) in the worst case.
    """
    n=len(A);m=len(B);cost=np.full((n,m),np.inf)
    pi=np.full((n,m),-1,np.int32);pj=np.full((n,m),-1,np.int32)
    d0=np.linalg.norm(A[0]-B[0])
    if minimax or d0<=delta3+1e-9:cost[0,0]=d0 if minimax else 1.
    for i in range(1,n):
        for j in range(1,m):
            d=np.linalg.norm(A[i]-B[j])
            if not minimax and d>delta3+1e-9:continue
            best=np.inf;bu=-1;bv=-1
            for u in range(i):
                if not EA[u,i]:continue
                for v in range(j):
                    if not EB[v,j]:continue
                    c=max(cost[u,v],d) if minimax else cost[u,v]+1
                    if c<best:best=c;bu=u;bv=v
            cost[i,j]=best;pi[i,j]=bu;pj[i,j]=bv
    final=cost[-1,-1]
    if not np.isfinite(final):return final,np.empty(0,np.int64),np.empty(0,np.int64)
    aa=[n-1];bb=[m-1]
    while aa[-1]>0 or bb[-1]>0:
        i=aa[-1];j=bb[-1];aa.append(pi[i,j]);bb.append(pj[i,j])
    return final,np.array(aa[::-1]),np.array(bb[::-1])

def joint_local(A,B,EA,EB,delta3):
    """Precheck before the cardinality solver; distinguish infeasibility from bugs."""
    floor=float(strict_pairing_floor(A,B))
    if delta3+1e-9<floor:
        return dict(status='infeasible_endpoint_floor',floor=floor,local_floor=None,ia=None,ib=None)
    local_floor,_,_=product_dag(A,B,EA,EB,np.inf,True)
    if delta3+1e-9<local_floor:
        return dict(status='infeasible_local_fidelity',floor=floor,local_floor=float(local_floor),ia=None,ib=None)
    k,ia,ib=product_dag(A,B,EA,EB,delta3,False)
    if not np.isfinite(k):raise AssertionError('Cardinality solver failed after a feasible bottleneck certificate')
    assert len(ia)==len(ib)==int(k)
    assert np.all(np.diff(ia)>0) and np.all(np.diff(ib)>0)
    assert ia[0]==ib[0]==0 and ia[-1]==len(A)-1 and ib[-1]==len(B)-1
    assert all(EA[i,j] for i,j in zip(ia,ia[1:]))
    assert all(EB[i,j] for i,j in zip(ib,ib[1:]))
    assert np.linalg.norm(A[ia]-B[ib],axis=1).max()<=delta3+1e-9
    return dict(status='optimal_verified',floor=floor,local_floor=float(local_floor),ia=ia,ib=ib,k=int(k))

def arc_length(P):return float(np.linalg.norm(np.diff(P,axis=0),axis=1).sum())

def edge_stats(A,B):
    a=np.linalg.norm(np.diff(A,axis=0),axis=1);b=np.linalg.norm(np.diff(B,axis=0),axis=1)
    return dict(mean=float(.5*(a.mean()+b.mean())),minimum=float(min(a.min(),b.min())),
                maximum=float(max(a.max(),b.max())),mean_A=float(a.mean()),mean_B=float(b.mean()),
                min_A=float(a.min()),min_B=float(b.min()),last_A=float(a[-1]),last_B=float(b[-1]))
