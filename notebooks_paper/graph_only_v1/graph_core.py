"""Graph-only thesis adapter; the existing optimization implementations are unchanged."""
from pathlib import Path
import os,sys,time,json
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'notebooks_paper'));sys.path.insert(0,str(ROOT))
import paper_core as previous
OUT=ROOT/'output/paper_v1/graph_only_v1'
os.environ['NUMBA_CACHE_DIR']=str(OUT/'cache')
os.environ['MPLCONFIGDIR']=str(OUT/'cache/matplotlib')
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
METHODS=previous.METHODS;COLORS=previous.COLORS;measures=previous.measures;sha=previous.sha
LIMITS=dict(seconds=8.,max_states=30000,max_transitions=500000,max_chain_edges=200000)
PROTOCOL=dict(certificate=False,hard_seconds=20.,startup_seconds=90.,memory_mib=512,
              solver_limits=LIMITS,repetitions=1,configurations=32)
TRACE={}

def solve(A,B,d1,d2,d3,method,*,certificate,progress=None):
    """Every caller must explicitly opt out; accepting True is forbidden here."""
    if certificate is not False:raise ValueError('This evaluation requires certificate=False')
    import cps_paper_algorithms as alg
    from curve_algorithms import shortest_independent
    from unittest.mock import patch
    TRACE.clear();TRACE.update(phase='preprocessing',chain_graph_calls=0,chain_graphs_completed=0,preprocessing_chain_graph_calls=0,
                              configuration_graph_entered=False,certificate_enabled=False)
    def emit():
        if progress:progress(dict(TRACE))
    emit()
    if method.startswith('Independent'):
        TRACE['phase']='independent_optimization';emit()
        fn=shortest_independent if method=='Independent continuous' else alg.independent_discrete
        ia,ib=fn(A,d1),fn(B,d2)
        return dict(status='optimal',raw_status='optimal_independent_baseline',route='independent_optimization',indices_A=ia,indices_B=ib,**TRACE)
    original=alg.chain_graph
    def observed_graph(*args,**kwargs):
        if sys._getframe(1).f_code is not alg.solve_cps.__code__:
            TRACE['preprocessing_chain_graph_calls']+=1
            return original(*args,**kwargs)
        TRACE['phase']='graph_construction';TRACE['chain_graph_calls']+=1;emit()
        graph=original(*args,**kwargs)
        TRACE['chain_graphs_completed']+=1
        if TRACE['chain_graphs_completed']==2:
            TRACE['phase']='graph_initialization'
        emit();return graph
    original_pop=alg.heapq.heappop
    def observed_pop(heap):
        if not TRACE['configuration_graph_entered'] and sys._getframe(1).f_code is alg.solve_cps.__code__:
            TRACE['phase']='configuration_graph';TRACE['configuration_graph_entered']=True;emit()
        return original_pop(heap)
    # Instrumentation only: no replacement of any graph transition or recurrence.
    with patch.object(alg,'chain_graph',observed_graph),patch.object(alg.heapq,'heappop',observed_pop):
        raw=alg.solve_cps(A,B,d1,d2,d3,method,certificate=False,**LIMITS)
    status=raw['status']
    assert 'certificate' not in status,status
    if status=='optimal_configuration_graph':
        assert TRACE['configuration_graph_entered'] and TRACE['chain_graphs_completed']==2
        route='configuration_graph';normalized='optimal'
    elif status=='infeasible_endpoints':route='preprocessing_exit';normalized='infeasible'
    elif status=='infeasible_graph':route='configuration_graph';normalized='infeasible'
    else:raise AssertionError('Unexpected solver status: '+status)
    raw.update(raw_status=status,status=normalized,route=route,**TRACE)
    return raw

def warmup(method):
    import numpy as np
    P=np.array([[0.,0.],[1.,0.],[2.,0.]])
    # Compilation is outside the measured solve. All CPS warm-up calls are graph-only too.
    return solve(P,P,.5,.5,1.,method,certificate=False)

def plot_domain(domain):
    import numpy as np,pandas as pd,matplotlib.pyplot as plt
    from unittest.mock import patch
    from IPython.display import display,Image
    previous.OUT=OUT
    # Reuse the paired comparison. Choose a useful geometric illustration separately.
    with patch.object(plt,'show',return_value=None):previous.plot_domain(domain)
    folder=OUT/'figures'/('proteins' if domain=='protein' else 'hurricanes')
    display(Image(filename=str(folder/'paired_comparison.png')))
    frame=pd.read_csv(OUT/f'{domain}_results.csv')
    coverage=frame.groupby('case',sort=False).status.apply(lambda x:x.eq('optimal').sum())
    case=coverage.idxmax();r=frame[frame.case.eq(case)].iloc[0]
    paths=json.loads((OUT/f'{domain}_paths.json').read_text())
    inp=np.load(ROOT/r.input_file);A,B=inp['A'],inp['B'];points=np.vstack([A,B])
    mid=(points.min(0)+points.max(0))/2;half=.55*np.ptp(points,axis=0).max()
    with plt.rc_context({'figure.facecolor':'white','axes.facecolor':'white','text.color':'black'}):
        fig=plt.figure(figsize=(12,4))
        for i,m in enumerate(METHODS):
            ax=fig.add_subplot(1,4,i+1,projection='3d' if A.shape[1]==3 else None)
            path=paths.get(case+'|'+m,{})
            for P,side,gray in [(A,'A','#333333'),(B,'B','#999999')]:
                ax.plot(*P.T,color=gray,alpha=.3,lw=1)
                ids=path.get('indices_'+side)
                if ids is not None:ax.plot(*P[ids].T,'o',ls='-' if side=='A' else '--',color=COLORS[m],ms=3)
            ax.set_xlim(mid[0]-half,mid[0]+half);ax.set_ylim(mid[1]-half,mid[1]+half)
            if A.shape[1]==3:ax.set_zlim(mid[2]-half,mid[2]+half);ax.set_box_aspect([1,1,1])
            else:ax.set_aspect('equal')
            ax.set_title(m+(' — missing' if not path else ''),fontsize=9);ax.tick_params(labelsize=7)
            ax.set_xlabel('x ('+r.units+')',fontsize=7);ax.set_ylabel('y ('+r.units+')',fontsize=7)
        fig.suptitle(f'{r.pair.replace("__"," / ")}; alpha={r.alpha:g}; delta=({r.delta1:.4g}, {r.delta2:.4g}, {r.delta3:.4g}) {r.units}\nIllustration chosen by completed-method coverage; cohort unchanged. Faint inputs; solid A / dashed B.',fontsize=10)
        fig.tight_layout(rect=[0,0,1,.83]);fig.savefig(folder/'geometry_example.png',dpi=180,bbox_inches='tight');plt.show();plt.close(fig)

compact_results=previous.compact_results
