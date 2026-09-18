"""Paper configuration DAGs; continuous GCS is in curve_algorithms.py.

CPS-3F: Fan et al. section 4, Algorithm 1 (explicit graph version).
CPS-2F: Anchored endpoints; asynchronous discrete coupling; objective max(kA,kB).
"""
import heapq
import itertools
import time
import numpy as np
from numba import njit
from curve_algorithms import (ball_segment_interval, segment_frechet_decision,
    continuous_decision, discrete_frechet, shortest_independent)

EPS = 1e-9

class ResourceLimit(RuntimeError):
    pass

def extended_curve(P, delta):
    """Original vertices plus exact sphere/edge intersections, ordered on P.

    Auxiliary points are owner positions ONLY, never output candidates.
    Store curve parameters to disambiguate coincident spatial points.
    """
    parameters = list(map(float, range(len(P))))
    for i in range(len(P)-1):
        a,b=P[i:i+2]; v=b-a; aa=float(v@v)
        if aa < 1e-24:
            continue  # zero edge has identical geometry; preserve its original ends
        for c in P:
            q=a-c; center=-float(q@v)/aa
            slack=delta*delta-float((q+center*v)@(q+center*v))
            if slack < -1e-10: continue
            radius=np.sqrt(max(0.,slack)/aa)
            for t in (center-radius,center+radius):
                if EPS<t<1-EPS: parameters.append(i+t)
    ordered=[]
    for t in sorted(parameters):
        if not ordered or t-ordered[-1]>1e-10: ordered.append(t)
    t=np.asarray(ordered); i=np.minimum(t.astype(int),len(P)-2)
    X=P[i]+(t-i)[:,None]*(P[i+1]-P[i])
    return X,t

@njit(cache=True)
def _discrete_edges(P,d,max_edges):
    n=len(P); edges=[]
    for i in range(n):
        for p in range(n):
            if np.linalg.norm(P[i]-P[p])>d+EPS: continue
            u=i*n+p
            for ii in range(i,min(n,i+2)):
                for pp in range(p,n):
                    if np.linalg.norm(P[ii]-P[pp])<=d+EPS:
                        edges.append((u,ii*n+pp,int(pp>p)))
                        if len(edges)>max_edges:return edges,False
    return edges,True

@njit(cache=True)
def _continuous_edges(P,X,d,max_edges):
    """Incremental exact curve/segment test; includes zero-cost owner motion."""
    n=len(P); h=len(X); edges=[]
    for p in range(n):
        for pp in range(p,n):
            lo=np.empty(h); hi=np.empty(h)
            for ii in range(h):
                lo[ii],hi[ii]=ball_segment_interval(P[p],P[pp],X[ii],d)
            for i in range(h):
                if np.linalg.norm(X[i]-P[p])>d+EPS: continue
                previous=0.
                for ii in range(i,h):
                    previous=max(previous,lo[ii])
                    if previous>hi[ii]+1e-10: break
                    if np.linalg.norm(X[ii]-P[pp])<=d+EPS:
                        edges.append((i*n+p,ii*n+pp,int(pp>p)))
                        if len(edges)>max_edges: return edges,False
    return edges,True

def chain_graph(P,d,kind,max_edges=2000000):
    P=np.ascontiguousarray(P,float)
    if kind=='CPS-2F':
        X,t=extended_curve(P,d)
        edges,finished=_continuous_edges(P,X,float(d),max_edges)
        if not finished: raise ResourceLimit('single-chain graph edge budget exceeded')
    elif kind=='CPS-3F':
        X=P; t=np.arange(len(P),dtype=float); edges,finished=_discrete_edges(P,float(d),max_edges)
        if not finished:raise ResourceLimit('single-chain graph edge budget exceeded')
    else: raise ValueError(kind)
    size=len(X)*len(P); out=[[] for _ in range(size)]; incoming=[[] for _ in range(size)]
    for u,v,c in edges:
        u,v,c=int(u),int(v),int(c)
        if u!=v: out[u].append((v,c));incoming[v].append(u)
    # Only states on an anchored source-to-terminal path can matter.
    forward=np.zeros(size,bool);forward[0]=True
    for u in range(size):
        if forward[u]:
            for v,c in out[u]: forward[v]=True
    backward=np.zeros(size,bool);backward[-1]=True
    for v in range(size-1,-1,-1):
        if backward[v]:
            for u in incoming[v]:backward[u]=True
    active=forward&backward
    for u in range(size):
        out[u]=([(u,0)]+[(v,c) for v,c in out[u] if active[v]]) if active[u] else []
    return dict(P=P,X=X,parameters=t,out=out,size=size,active=active,
                edges=sum(map(len,out)),auxiliary=len(X)-len(P))

def independent_discrete(P,d):
    """Exact shortest anchored discrete simplification, single configuration DAG."""
    g=chain_graph(P,d,'CPS-3F'); n=len(P)
    cost=np.full(g['size'],n+1,int); pred=np.full(g['size'],-1,int);cost[0]=0
    for u in range(g['size']):
        for v,c in g['out'][u]:
            if u!=v and cost[u]+c<cost[v]:cost[v]=cost[u]+c;pred[v]=u
    u=g['size']-1; path=[]
    while u>=0:path.append(u%n);u=pred[u]
    return np.asarray(list(dict.fromkeys(path[::-1])),int)

def discrete_coupling(A,B):
    n,m=len(A),len(B);D=np.full((n,m),np.inf);pred={}
    for i in range(n):
        for j in range(m):
            if i==j==0:D[i,j]=np.linalg.norm(A[i]-B[j]);continue
            candidates=[(D[u,v],u,v) for u,v in ((i-1,j),(i,j-1),(i-1,j-1)) if u>=0 and v>=0]
            value,u,v=min(candidates);D[i,j]=max(value,np.linalg.norm(A[i]-B[j]));pred[i,j]=(u,v)
    path=[(n-1,m-1)]
    while path[-1]!=(0,0):path.append(pred[path[-1]])
    return np.asarray(path[::-1],int),float(D[-1,-1])

def alternate_certificate(A,B,d1,d2,d3,ia,ib,limit=150000):
    """Optional exact lower-bound witness search; does not claim graph execution.

    Any feasible pair with max count equal to the independent lower bound is
    globally optimal. Failure to find one proves nothing and falls back to DP.
    """
    bound=max(len(ia),len(ib));tested=0;sets=[]
    for P,Q,d,base,other in [(A,B,d1,ia,ib),(B,A,d2,ib,ia)]:
        choices=[base]
        for k in range(len(base),bound+1):
            for middle in itertools.combinations(range(1,len(P)-1),k-2):
                ids=np.array((0,)+middle+(len(P)-1,));tested+=1
                if tested>limit:return None,tested
                if not continuous_decision(P,P[ids],d):continue
                if discrete_frechet(P[ids],Q[other])<=d3+EPS:
                    return ((ids,other) if len(sets)==0 else (other,ids)),tested
                choices.append(ids)
        sets.append(choices)
    for a in sets[0]:
        for b in sets[1]:
            tested+=1
            if tested>limit:return None,tested
            if discrete_frechet(A[a],B[b])<=d3+EPS:return (a,b),tested
    return None,tested

def solve_cps(A,B,d1,d2,d3,kind,*,certificate=True,max_states=300000,
              max_transitions=15000000,seconds=90.,max_chain_edges=2000000):
    """Same min/max-count recurrence as paper DP, sparse reachable states.

    X[state][r] stores minimum B hops at EXACT A hop count r. This is the
    equivalent exact-count formulation; a single scalar cost is insufficient.
    All non-self product transitions advance topological key u+v.
    Resource exits never claim infeasibility or optimality.
    """
    start=time.perf_counter(); A=np.ascontiguousarray(A,float);B=np.ascontiguousarray(B,float)
    floor=max(np.linalg.norm(A[0]-B[0]),np.linalg.norm(A[-1]-B[-1]))
    if floor>d3+EPS:return dict(status='infeasible_endpoints',floor=float(floor))
    single=shortest_independent if kind=='CPS-2F' else independent_discrete
    ia=single(A,d1);ib=single(B,d2); lower=max(len(ia),len(ib))
    if certificate and discrete_frechet(A[ia],B[ib])<=d3+EPS:
        return finish(A,B,ia,ib,d1,d2,d3,kind,'optimal_independent_certificate',start,
                      lower_bound=lower,states=0,transitions=0)
    if certificate and kind=='CPS-2F':
        alternate,tests=alternate_certificate(A,B,d1,d2,d3,ia,ib)
        if alternate is not None:
            return finish(A,B,*alternate,d1,d2,d3,kind,'optimal_alternate_certificate',start,
                          lower_bound=lower,states=0,transitions=0,certificate_tests=tests)
    ga=chain_graph(A,d1,kind,max_chain_edges);gb=chain_graph(B,d2,kind,max_chain_edges)
    na,nb=len(A),len(B);end=(ga['size']-1,gb['size']-1)
    allowed=np.linalg.norm(A[:,None]-B[None,:],axis=2)<=d3+EPS
    values={(0,0):{0:0}}; parents={};heap=[(0,0,0)];queued={(0,0)};transitions=0
    while heap:
        _,u,v=heapq.heappop(heap); state=(u,v)
        if state==end:break
        current=values[state]
        for uu,ca in ga['out'][u]:
            for vv,cb in gb['out'][v]:
                if (uu==u and vv==v) or not allowed[uu%na,vv%nb]:continue
                transitions+=1
                if transitions>max_transitions:raise ResourceLimit('product transition budget exceeded')
                target=(uu,vv)
                if target not in values:
                    if len(values)>=max_states:raise ResourceLimit('reachable configuration budget exceeded')
                    values[target]={}
                dest=values[target]
                for r,z in current.items():
                    rr,zz=r+ca,z+cb
                    if zz<dest.get(rr,nb+1):dest[rr]=zz;parents[(uu,vv,rr)]=(u,v,r)
                if target not in queued:
                    heapq.heappush(heap,(uu+vv,uu,vv));queued.add(target)
        if time.perf_counter()-start>seconds:raise ResourceLimit('wall-time budget exceeded')
    if end not in values:return dict(status='infeasible_graph',floor=float(floor),states=len(values),transitions=transitions)
    r=min(values[end],key=lambda r:(max(r,values[end][r]),r+values[end][r],r))
    trace=[(*end,r)]
    while trace[-1]!=(0,0,0):trace.append(parents[trace[-1]])
    trace.reverse();ia=np.array(list(dict.fromkeys(u%na for u,v,r in trace)),int)
    ib=np.array(list(dict.fromkeys(v%nb for u,v,r in trace)),int)
    return finish(A,B,ia,ib,d1,d2,d3,kind,'optimal_configuration_graph',start,
                  lower_bound=lower,states=len(values),transitions=transitions,
                  auxiliary_A=ga['auxiliary'],auxiliary_B=gb['auxiliary'])

def finish(A,B,ia,ib,d1,d2,d3,kind,status,start,**stats):
    assert ia[0]==ib[0]==0 and ia[-1]==len(A)-1 and ib[-1]==len(B)-1
    assert np.all(np.diff(ia)>0) and np.all(np.diff(ib)>0)
    if kind=='CPS-2F':
        assert continuous_decision(A,A[ia],d1) and continuous_decision(B,B[ib],d2)
    else:assert discrete_frechet(A,A[ia])<=d1+EPS and discrete_frechet(B,B[ib])<=d2+EPS
    coupling,d=discrete_coupling(A[ia],B[ib]);assert d<=d3+EPS
    return dict(status=status,indices_A=ia,indices_B=ib,coupling=coupling,k=max(len(ia),len(ib)),
                dDF=d,verified=True,seconds=time.perf_counter()-start,**stats)
