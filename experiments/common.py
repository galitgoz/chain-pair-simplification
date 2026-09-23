"""Orchestration and validation only; existing optimization algorithms remain unchanged."""
from pathlib import Path
import os,sys,json,hashlib,time,datetime
ROOT=Path(__file__).resolve().parents[1];OUT=Path(os.environ.get('CPS_EXPERIMENT_OUTPUT',str(ROOT/'reproduced/experiments'))).resolve()
sys.path[:0] = [str(ROOT / 'src'), str(ROOT)]
os.environ['NUMBA_CACHE_DIR']=str(OUT/'cache');os.environ['MPLCONFIGDIR']=str(OUT/'cache/matplotlib')
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def stamp():return dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),unix=time.time(),monotonic=time.monotonic())
def serial(value):
    if hasattr(value,'tolist'):return value.tolist()
    if isinstance(value,Path):return str(value)
    raise TypeError(type(value).__name__)
def atomic(path,data):
    path=Path(path);tmp=path.with_name(path.name+'.'+str(os.getpid())+'.tmp')
    with tmp.open('w',encoding='utf-8') as f:
        json.dump(data,f,indent=2,default=serial,allow_nan=False);f.flush();os.fsync(f.fileno())
    # Windows readers/antivirus may briefly deny replacement of an open file.
    # Retry this I/O operation for at most 0.7 s; all phase budgets still apply.
    for attempt_number in range(15):
        try:
            os.replace(tmp,path);break
        except PermissionError:
            if attempt_number==14:raise
            time.sleep(.05)

BUDGETS=dict(startup=60.,input_loading=15.,warmup=120.,solver=60.,checkpoint=5.,validation=30.,optional_measurements=30.,reporting=10.)
BASE_GRAPH_LIMITS=dict(seconds=55.,max_states=30000,max_transitions=500000,max_chain_edges=200000)

class Phases:
    def __init__(self,directory):self.directory=Path(directory);self.phase='startup';self.start=time.monotonic()
    def begin(self,phase,**details):
        now=stamp();self.phase=phase;self.start=now['monotonic']
        event=dict(phase=phase,phase_started_monotonic=self.start,**now,**details)
        with (self.directory/'phases.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(event,default=serial)+'\n');f.flush()
        atomic(self.directory/'progress.json',event)
    def note(self,**details):
        event=dict(phase=self.phase,phase_started_monotonic=self.start,**stamp(),**details)
        with (self.directory/'phases.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(event,default=serial)+'\n');f.flush()
        atomic(self.directory/'progress.json',event)

def classify(error):
    text=str(error)
    if 'product transition budget' in text:return 'transition_budget_limit'
    if 'single-chain graph edge budget' in text:return 'single_chain_edge_budget_limit'
    if 'reachable configuration budget' in text:return 'state_budget_limit'
    if 'wall-time budget' in text:return 'solver_timeout'
    if isinstance(error,MemoryError):return 'memory_limit'
    return 'software_error'

def exception_counters(error):
    result={};tb=error.__traceback__
    while tb:
        f=tb.tb_frame;local=f.f_locals
        if f.f_code.co_name=='solve_cps':
            if 'transitions' in local:result['product_transitions']=int(local['transitions'])
            if 'values' in local:result['states_discovered']=len(local['values'])
            if 'heap' in local:result['queue_size']=len(local['heap'])
        if f.f_code.co_name=='solve' and 'trace' in local:
            for key in ['states_dequeued','completed_states_lower_bound']:
                if key in local['trace']:result[key]=local['trace'][key]
        if f.f_code.co_name=='chain_graph' and 'edges' in local:
            result['single_chain_edges_generated']=len(local['edges'])
        tb=tb.tb_next
    return result

def solve(A,B,config,*,certificate,observe=None):
    """Explicit False at the actual CPS call; observes phases without changing recurrence."""
    assert certificate is False
    import cps_paper_algorithms as alg
    from curve_algorithms import shortest_independent
    from unittest.mock import patch
    method=config['method'];lim=config['graph_limits']
    observed_start=time.perf_counter();last_report=observed_start
    trace=dict(solver_subphase='independent_optimization' if method.startswith('Independent') else 'preprocessing',
               configuration_graph_entered=False,chain_graph_calls=0,chain_graphs_completed=0,
               states_dequeued=0,completed_states_lower_bound=0)
    def emit():
        if observe:observe(dict(trace))
    emit();original=alg.chain_graph;pop=alg.heapq.heappop
    def graph(*args,**kwargs):
        direct=sys._getframe(1).f_code is alg.solve_cps.__code__
        if direct:trace.update(solver_subphase='graph_construction',chain_graph_calls=trace['chain_graph_calls']+1);emit()
        g=original(*args,**kwargs)
        if direct:
            trace['chain_graphs_completed']+=1
            trace['last_completed_chain_edges']=g['edges'];trace['last_completed_chain_auxiliary']=g['auxiliary']
            # Boundary-only metadata: no changes to graph construction or search.
            side='A' if trace['chain_graphs_completed']==1 else 'B'
            trace['single_graph_'+side]=dict(dense_states=int(g['size']),
                active_states=int(g['active'].sum()),edges_including_waits=int(g['edges']),
                auxiliary=int(g['auxiliary']),extended_vertices=len(g['X']),
                original_vertices=len(args[0]),construction_completed=True)
            if trace['chain_graphs_completed']==2:trace['solver_subphase']='graph_initialization'
            emit()
        return g
    def observed_pop(heap):
        nonlocal last_report
        frame=sys._getframe(1)
        if frame.f_code is not alg.solve_cps.__code__:return pop(heap)
        item=pop(heap)  # unchanged heap operation and return value
        first=not trace['configuration_graph_entered']
        trace['completed_states_lower_bound']=trace['states_dequeued']
        trace['states_dequeued']+=1
        trace.update(solver_subphase='configuration_graph',configuration_graph_entered=True,
                     queue_size=len(heap),product_transitions=int(frame.f_locals.get('transitions',0)),
                     states_discovered=len(frame.f_locals.get('values',{})))
        now=time.perf_counter()
        if first or now-last_report>=1.:
            trace['observed_solver_seconds']=now-observed_start;emit();last_report=now
        return item
    with patch.object(alg,'chain_graph',graph),patch.object(alg.heapq,'heappop',observed_pop):
        if method.startswith('Independent'):
            fn=shortest_independent if method=='Independent continuous' else alg.independent_discrete
            ia,ib=fn(A,config['delta1']),fn(B,config['delta2'])
            raw=dict(status='optimal_independent_baseline',indices_A=ia,indices_B=ib,k=max(len(ia),len(ib)))
        else:
            raw=alg.solve_cps(A,B,config['delta1'],config['delta2'],config['delta3'],method,certificate=False,**lim)
    # No external feasibility or achieved-distance checks before returning this result.
    status=raw['status']
    if status in ['optimal_configuration_graph','infeasible_graph']:
        trace.update(product_transitions=raw.get('transitions'),states_discovered=raw.get('states'),
                     observed_solver_seconds=time.perf_counter()-observed_start)
    route='configuration_graph' if status in ['optimal_configuration_graph','infeasible_graph'] else 'preprocessing_exit' if status=='infeasible_endpoints' else 'independent_optimization'
    return raw,dict(route=route,**trace)

def warmup(dim,method):
    import numpy as np
    import cps_paper_algorithms as alg
    from curve_algorithms import continuous_decision,discrete_frechet,continuous_distance
    P=np.zeros((3,dim));P[:,0]=[0,1,2];P[1,1]=1
    c=dict(method=method,delta1=.6,delta2=.6,delta3=2.,graph_limits={**BASE_GRAPH_LIMITS,'seconds':float('inf')})
    solve(P,P,c,certificate=False)
    continuous_decision(P,P,0.);discrete_frechet(P,P);alg.discrete_coupling(P,P)
    continuous_distance(P,P)  # optional-distance code compilation also occurs outside solver timing

def validate(A,B,raw,config):
    """Feasibility only. Does not independently certify the solver's optimality claim."""
    import numpy as np
    import cps_paper_algorithms as alg
    from curve_algorithms import continuous_decision,discrete_frechet
    # Preserve every coordinate and predicate, but use the C-layout signature
    # already compiled during warmup rather than compiling F/C during validation.
    A=np.ascontiguousarray(A);B=np.ascontiguousarray(B)
    checks={};chosen=[]
    for side,P in [('A',A),('B',B)]:
        values=raw['indices_'+side]
        assert isinstance(values,(list,np.ndarray)) and len(values)>0,'empty or missing indices'
        ids=np.asarray(values)
        assert np.issubdtype(ids.dtype,np.integer),'indices must be integers'
        assert ((ids>=0)&(ids<len(P))).all(),'index bounds'
        assert (np.diff(ids)>0).all(),'strict index order'
        assert ids[0]==0 and ids[-1]==len(P)-1,'endpoints'
        checks[side+'_original_vertex_membership']=True;chosen.append(P[ids])
    a,b=chosen;counts=dict(kA=len(a),kB=len(b),k=max(len(a),len(b)))
    assert raw['k']==counts['k'],'reported objective differs from selected vertex counts'
    for field in ['kA','kB']:
        if field in raw:assert raw[field]==counts[field]
    discrete=config['method'] in ['CPS-3F','Independent discrete']
    for side,P,Q,d in [('A',A,a,config['delta1']),('B',B,b,config['delta2'])]:
        passed=discrete_frechet(P,Q)<=d+1e-9 if discrete else continuous_decision(P,Q,d)
        assert passed,side+' fidelity constraint';checks[side+'_fidelity']=True
    coupling,distance=alg.discrete_coupling(a,b)
    def coupling_check(path):
        path=np.asarray(path)
        assert path.ndim==2 and path.shape[1]==2 and np.issubdtype(path.dtype,np.integer),'coupling structure'
        assert ((path>=0)&(path<np.array([len(a),len(b)]))).all(),'coupling bounds'
        assert np.array_equal(path[0],[0,0]) and np.array_equal(path[-1],[len(a)-1,len(b)-1]),'coupling endpoints'
        steps=np.diff(path,axis=0)
        assert ((steps>=0)&(steps<=1)).all() and (steps.sum(1)>=1).all(),'coupling monotonicity'
        return float(np.linalg.norm(a[path[:,0]]-b[path[:,1]],axis=1).max())
    computed_cost=coupling_check(coupling)
    assert abs(computed_cost-distance)<1e-7
    joint=config['method'].startswith('CPS')
    if 'coupling' in raw:
        reported_cost=coupling_check(raw['coupling'])
        if joint:assert reported_cost<=config['delta3']+1e-9,'reported coupling budget'
    if joint:assert distance<=config['delta3']+1e-9,'coupling budget'
    return dict(status='validated',feasibility_passed=True,optimality_independently_proven=False,
        solver_optimality_claim=raw['status'].startswith('optimal'),checks=checks,**counts,
        coupling_distance=distance,coupling_feasible=bool(distance<=config['delta3']+1e-9),
        coupling_constraint_required=joint,coupling=coupling.tolist())
