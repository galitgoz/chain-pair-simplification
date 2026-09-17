"""Reuse pure existing oracles; never call legacy test runners that write old results."""
import itertools,json,time,traceback
import numpy as np
import pandas as pd
from paper_core import ROOT,OUT,METHODS,solve,measures
from verify_cps_papers import candidates
from verify_algorithms import brute
from curve_algorithms import (shortest_independent,continuous_decision,discrete_frechet,
                              segment_frechet_decision,local_shortcut_path)
from cps_paper_algorithms import extended_curve,ResourceLimit

def exhaustive_outputs(P,d,kind):
    return [np.array([0])] if len(P)==1 else list(candidates(P,d,kind))

def check_all():
    rows=[]
    def record(name,method,fn,kind='targeted'):
        start=time.perf_counter()
        try:detail=fn();outcome='PASS';reason=''
        except ResourceLimit as error:detail='validation did not finish';outcome='INCOMPLETE';reason=str(error)
        except Exception:detail='';outcome='FAIL';reason=traceback.format_exc()
        rows.append(dict(check=name,affected_method=method,check_kind=kind,outcome=outcome,
                         failure_reason=reason,details=str(detail),seconds=time.perf_counter()-start))
        pd.DataFrame(rows).to_csv(OUT/'validation.csv',index=False)
    def compare(A,B,d1,d2,d3,method):
        kind='CPS-2F' if method in ['Independent continuous','CPS-2F'] else 'CPS-3F'
        ca=exhaustive_outputs(A,d1,kind);cb=exhaustive_outputs(B,d2,kind)
        if method.startswith('Independent'):expected=max(min(map(len,ca)),min(map(len,cb)))
        else:expected=min((max(len(a),len(b)) for a in ca for b in cb if discrete_frechet(A[a],B[b])<=d3+1e-9),default=None)
        result=solve(A,B,d1,d2,d3,method,certificate=False)
        if expected is None:
            assert result['status']=='infeasible',result
            return 'Exhaustive infeasibility agrees; no resource exit treated as infeasible.'
        assert result['status']=='optimal',result
        actual=measures(A,B,result,d1,d2,d3,method)
        assert actual['k']==expected,(actual['k'],expected)
        if method.startswith('CPS'):assert result['route']=='configuration_graph',result
        return f"optimum k={expected}; endpoints, counts, fidelity, output membership, monotone coupling validated; route={result['route']}"
    # Analytic geometry controls independent of global-DP optimization.
    line=np.array([[0.,0],[1,0],[2,0]])
    def analytic():
        assert continuous_decision(line,line[[0,2]],0)
        assert discrete_frechet(line,line[[0,2]])==1
        tangent=np.array([[0.,0],[1,1],[2,0]])
        assert continuous_decision(tangent,tangent[[0,2]],1)
        assert not continuous_decision(tangent,tangent[[0,2]],.99)
        back=np.array([[0.,0],[2,0],[1,0],[3,0]])
        assert not continuous_decision(back,back[[0,-1]],.49)
        assert continuous_decision(back,back[[0,-1]],.5)
        return 'Known thresholds: straight 0, discrete straight 1, tangency 1, backtracking .5.'
    record('analytic_frechet_controls','all',analytic)
    def grid_oracle():
        A=np.array([[0.,0],[1,1],[2,0]]);B=np.array([[0.,0],[.5,0],[2,0]])
        values=[]
        def visit(i,j,peak):
            peak=max(peak,float(np.linalg.norm(A[i]-B[j])))
            if (i,j)==(2,2):values.append(peak);return
            for di,dj in [(1,0),(0,1),(1,1)]:
                if i+di<3 and j+dj<3:visit(i+di,j+dj,peak)
        visit(0,0,0);assert np.isclose(min(values),discrete_frechet(A,B))
        return f'Enumerated {len(values)} monotone grid paths.'
    record('discrete_distance_all_grid_paths','all',grid_oracle)
    rng=np.random.default_rng(20260916)
    for dim in [2,3]:
        for trial in range(8):
            A=np.cumsum(rng.normal(size=(int(rng.integers(3,6)),dim)),axis=0)
            B=np.cumsum(rng.normal(size=(int(rng.integers(3,6)),dim)),axis=0)
            d1,d2=rng.uniform(.15,2.,2)
            d3=float(discrete_frechet(A,B)+.05) if trial%3 else max(np.linalg.norm(A[0]-B[0]),np.linalg.norm(A[-1]-B[-1]))-.05
            for method in METHODS:
                record(f'exhaustive_{dim}d_{trial:02d}',method,lambda A=A,B=B,d1=d1,d2=d2,d3=d3,m=method:compare(A,B,d1,d2,d3,m),'exhaustive; forced graph for joint methods')
    for name,A,B,d1,d2,d3 in [
        ('zero_threshold_count_vertices',line,line,0.,0.,0.),
        ('singleton',line[:1],line[:1],0.,0.,0.),
        ('duplicate_vertices',np.array([[0.,0],[0,0],[1,0],[1,0]]),line,0.,0.,1.),
        ('unequal_lengths',line,np.array([[0.,0],[1,1],[2,0]]),0.,0.,2.),
        ('infeasible_endpoints',line,line+10,0.,0.,1.),
        ('infeasible_internal',line,np.array([[0.,0],[1,5],[2,0]]),0.,0.,.1)]:
        for m in METHODS:record(name,m,lambda A=A,B=B,d1=d1,d2=d2,d3=d3,m=m:compare(A,B,d1,d2,d3,m))
    for m in METHODS:
        def repeated(method=m):
            B=np.array([[0.,0],[1,1],[2,0]])
            result=solve(line,B,0,0,2,method,certificate=False)
            evidence=measures(line,B,result,0,0,2,method)
            if method in ['Independent continuous','CPS-2F']:
                assert evidence['unequal_output_lengths'] and evidence['repeated_coupling_index']
            else:
                B=np.array([[0.,0],[.5,0],[1,0],[1.5,0],[2,0]])
                result=solve(line,B,0,0,.5,method,certificate=False)
                evidence=measures(line,B,result,0,0,.5,method)
                assert evidence['unequal_output_lengths'] and evidence['repeated_coupling_index']
            return f"kA={evidence['kA']}, kB={evidence['kB']}, coupling repeats an index."
        record('unequal_outputs_repeated_coupling',m,repeated)
    prior=json.loads((ROOT/'output/algorithm_checks.json').read_text())['local_global_counterexample']
    P=np.array(prior['P']);d=prior['delta']
    def global_regression():
        ids=shortest_independent(P,d);oracle=brute(P,d)
        assert len(ids)==len(oracle)==3 and len(local_shortcut_path(P,d))==4
        return 'Saved counterexample re-executed: global optimum 3 versus local shortcut count 4. Local tests are only a feasible upper bound internally.'
    record('global_not_local_regression','Independent continuous',global_regression)
    record('global_not_local_forced_graph','CPS-2F',lambda:compare(P,P,d,d,0,'CPS-2F'))
    def auxiliary():
        P=np.array([[0.,0],[1,1],[2,0]])
        X,t=extended_curve(P,.6)
        assert len(X)>len(P) and np.all(np.diff(t)>0)
        i=np.minimum(t.astype(int),len(P)-2)
        assert np.allclose(X,P[i]+(t-i)[:,None]*(P[i+1]-P[i]))
        result=solve(P,P,.6,.6,0,'CPS-2F',certificate=False)
        assert np.issubdtype(np.asarray(result['indices_A']).dtype,np.integer)
        assert max(result['indices_A'])<len(P)
        return 'Auxiliary positions subdivide inputs; output indices remain original-vertex indices.'
    record('auxiliary_not_output_candidates','CPS-2F',auxiliary)
    for m in ['CPS-2F','CPS-3F']:
        def route(method=m):
            result=solve(line,line,1,1,1,method,certificate=True)
            assert result['route']=='independent_certificate'
            assert result.get('states')==result.get('transitions')==0
            return 'Independent certificate is labeled separately; graph is not reported as executed.'
        record('certificate_route_label',m,route)
    report=pd.DataFrame(rows)
    (OUT/'validation_gate.json').write_text(json.dumps(dict(checks=len(rows),passed=int(report.outcome.eq('PASS').sum()),
        failed=int(report.outcome.eq('FAIL').sum()),incomplete=int(report.outcome.eq('INCOMPLETE').sum()),note='Finite tests are implementation evidence, not a proof of the draft discretization theorem.'),indent=2))
    return report
