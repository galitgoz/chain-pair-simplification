"""Continuous global vertex-restricted simplification and discrete CPS references.

The GCS solver represents Algorithm 1's cost lower envelope by its integer
cost layers: reach[v,k,i] is the earliest reachable point on original edge i
at spine v using exactly k shortcut links. No shortcut/subcurve anchoring.
"""
import itertools
import numpy as np
from numba import njit

TOL = 1e-10

@njit(cache=True)
def ball_segment_interval(a,b,c,delta):
    """Closed parameter interval where segment a+t(b-a) is within delta of c."""
    d=b-a; q=a-c
    aa=np.dot(d,d)
    if aa<1e-24:
        if np.dot(q,q)<=delta*delta+TOL: return 0.,1.
        return np.inf,-np.inf
    center=-np.dot(q,d)/aa
    closest=q+center*d
    slack=delta*delta-np.dot(closest,closest)
    if slack < -TOL: return np.inf,-np.inf
    radius=np.sqrt(max(0.,slack)/aa)
    lo=max(0.,center-radius); hi=min(1.,center+radius)
    if lo>hi+TOL: return np.inf,-np.inf
    return min(lo,hi),hi

@njit(cache=True)
def segment_frechet_decision(P,a,b,delta):
    """Exact continuous decision for a polygonal curve versus one segment."""
    if np.linalg.norm(P[0]-a)>delta+TOL or np.linalg.norm(P[-1]-b)>delta+TOL: return False
    previous=0.
    for p in P:
        lo,hi=ball_segment_interval(a,b,p,delta)
        previous=max(previous,lo)
        if previous>hi+TOL: return False
    return True

@njit(cache=True)
def continuous_decision(P,Q,delta):
    """Full continuous free-space reachability, including degenerate edges."""
    n=len(P); m=len(Q)
    if n==1: return np.max(np.sqrt(np.sum((Q-P[0])**2,axis=1)))<=delta+TOL
    if m==1: return np.max(np.sqrt(np.sum((P-Q[0])**2,axis=1)))<=delta+TOL
    if np.linalg.norm(P[0]-Q[0])>delta+TOL or np.linalg.norm(P[-1]-Q[-1])>delta+TOL: return False
    right=np.full((n-1,m-1),np.inf)
    top=np.full((n-1,m-1),np.inf)
    for i in range(n-1):
        for j in range(m-1):
            if i==0:
                lo,hi=ball_segment_interval(Q[j],Q[j+1],P[0],delta)
                left=0. if lo<=TOL and (j==0 or top[0,j-1]<=TOL) else np.inf
            else: left=right[i-1,j]
            if j==0:
                lo,hi=ball_segment_interval(P[i],P[i+1],Q[0],delta)
                bottom=0. if lo<=TOL and (i==0 or right[i-1,0]<=TOL) else np.inf
            else: bottom=top[i,j-1]
            rlo,rhi=ball_segment_interval(Q[j],Q[j+1],P[i+1],delta)
            tlo,thi=ball_segment_interval(P[i],P[i+1],Q[j+1],delta)
            r=rlo if np.isfinite(bottom) else max(rlo,left)
            t=tlo if np.isfinite(left) else max(tlo,bottom)
            if r<=rhi+TOL: right[i,j]=r
            if t<=thi+TOL: top[i,j]=t
    return right[-1,-1]<=1.+TOL or top[-1,-1]<=1.+TOL

def continuous_distance(P,Q,absolute_tolerance=1e-5):
    """Certified decision bracket [lower, upper], not an exact critical value."""
    P=np.asarray(P,float); Q=np.asarray(Q,float)
    lower=max(np.linalg.norm(P[0]-Q[0]),np.linalg.norm(P[-1]-Q[-1]))
    upper=discrete_frechet(P,Q)
    if continuous_decision(P,Q,lower): return float(lower),float(lower)
    while upper-lower>absolute_tolerance:
        mid=(lower+upper)/2
        if continuous_decision(P,Q,mid): upper=mid
        else: lower=mid
    return float(lower),float(upper)

@njit(cache=True)
def local_shortcut_path(P,delta):
    """Restricted local baseline; used ONLY as a feasible link-count upper bound."""
    n=len(P); cost=np.full(n,n+1,np.int32); pred=np.full(n,-1,np.int32); cost[0]=0
    for v in range(1,n):
        for u in range(v):
            if cost[u]+1<cost[v] and segment_frechet_decision(P[u:v+1],P[u],P[v],delta):
                cost[v]=cost[u]+1; pred[v]=u
    path=[n-1]
    while path[-1]>0: path.append(pred[path[-1]])
    return np.array(path[::-1])

@njit(cache=True)
def _gcs(P,delta,cap):
    n=len(P); rows=n-1
    spine_lo=np.empty((n,rows)); spine_hi=np.empty((n,rows))
    for v in range(n):
        for i in range(rows):
            spine_lo[v,i],spine_hi[v,i]=ball_segment_interval(P[i],P[i+1],P[v],delta)
    reach=np.full((n,cap+1,rows),np.inf)
    parent_u=np.full((n,cap+1,rows),-1,np.int32)
    parent_row=np.full((n,cap+1,rows),-1,np.int32)
    # Starting spine: move along P while the simplification stays at P[0].
    for i in range(rows):
        if spine_lo[0,i]>TOL: break
        reach[0,0,i]=0.
        if spine_hi[0,i]<1.-TOL: break
    for v in range(1,n):
        for u in range(v):
            top_lo=np.empty(rows); top_hi=np.empty(rows)
            for i in range(rows):
                top_lo[i],top_hi[i]=ball_segment_interval(P[u],P[v],P[i+1],delta)
            for k in range(1,min(cap,v)+1):
                # Each strip cell is convex. Bottom arrival frees all of its right
                # boundary; left arrival constrains its right boundary from below.
                bottom=np.inf; origin=-1
                for i in range(rows):
                    left=reach[u,k-1,i]
                    if np.isfinite(bottom):
                        r=spine_lo[v,i]; source=origin
                    else:
                        r=max(spine_lo[v,i],left); source=i
                    if r<=spine_hi[v,i]+TOL and r<reach[v,k,i]:
                        reach[v,k,i]=r
                        parent_u[v,k,i]=u; parent_row[v,k,i]=source
                    if np.isfinite(left):
                        t=top_lo[i]; origin=i
                    else: t=max(top_lo[i],bottom)
                    bottom=t if t<=top_hi[i]+TOL else np.inf
        # A terminal point may lie at the top of the last original edge.
    best=cap
    for k in range(1,cap+1):
        if reach[n-1,k,rows-1]<=1.+TOL:
            best=k; break
    if not np.isfinite(reach[n-1,best,rows-1]): raise ValueError('No reachable terminal state')
    path=np.empty(best+1,np.int64); path[best]=n-1
    v=n-1; i=rows-1
    for k in range(best,0,-1):
        u=parent_u[v,k,i]; r=parent_row[v,k,i]
        path[k-1]=u; v=u; i=r
    return path

def shortest_independent(P,delta):
    """Shortest endpoint-preserving vertex subsequence under global continuous F.

    Algorithm 1 / Lemma 4 of the supplied full GCS paper, using integer cost
    layers instead of explicitly subdividing all elementary intervals. The layers
    here represent EXACT link counts (the minimum over them is the lower envelope).
    Time O(n^3 K), memory O(n^2 K), K <= n-1 is a feasible local upper bound.
    Floating-point geometric predicates have TOL=1e-10; no resampling is used.
    """
    P=np.ascontiguousarray(P,dtype=float)
    if delta<0 or not np.isfinite(delta) or len(P)==0 or not np.isfinite(P).all(): raise ValueError('Invalid curve/tolerance')
    if len(P)==1: return np.array([0],dtype=int)
    local=local_shortcut_path(P,delta)
    result=_gcs(P,float(delta),len(local)-1)
    if not continuous_decision(P,P[result],delta): raise AssertionError('Independent continuous validation failed')
    return result

@njit(cache=True)
def discrete_frechet(P,Q):
    prev=np.full(len(Q),np.inf)
    for i in range(len(P)):
        current=np.full(len(Q),np.inf)
        for j in range(len(Q)):
            if i==0 and j==0: preceding=0.
            else:
                preceding=prev[j]
                if j>0: preceding=min(preceding,current[j-1],prev[j-1])
            current[j]=max(preceding,np.linalg.norm(P[i]-Q[j]))
        prev=current
    return prev[-1]

def cps_discrete_reference(A,B,delta1,delta2,delta3,fix_endpoints=False,max_vertices=10):
    """Exact exhaustive small-instance CPS-3F oracle; not scalable Algorithm 2.

    Fan et al. allow the first/last selected vertices to differ from the input
    endpoints. They still require all original vertices to be traversed by each
    discrete fidelity coupling. Minimize max(len(A'), len(B')).
    """
    if max(len(A),len(B))>max_vertices:
        raise ValueError('Reference CPS is exponential; use only small contiguous windows')
    def candidates(P,d):
        result=[]
        for size in range(1,len(P)+1):
            for ids in itertools.combinations(range(len(P)),size):
                if fix_endpoints and (ids[0]!=0 or ids[-1]!=len(P)-1): continue
                if discrete_frechet(P,P[list(ids)])<=d+TOL: result.append(ids)
        return result
    ca,cb=candidates(A,delta1),candidates(B,delta2)
    for k in range(1,max(len(A),len(B))+1):
        for ia in ca:
            if len(ia)>k: continue
            for ib in cb:
                if max(len(ia),len(ib))!=k: continue
                if discrete_frechet(A[list(ia)],B[list(ib)])<=delta3+TOL:
                    return dict(indices_A=np.array(ia),indices_B=np.array(ib),k=k,
                                method='exact exhaustive discrete CPS-3F reference')
    return None

def cps_2f(*args,**kwargs):
    raise NotImplementedError('CPS-2F is reserved; its constraints and objective must be specified before implementation.')
