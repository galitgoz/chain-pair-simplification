"""Measured comparison with explicit solver certificates and resource statuses."""
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from curve_algorithms import shortest_independent,continuous_distance,discrete_frechet
from cps_paper_algorithms import solve_cps,ResourceLimit,discrete_coupling

def arc(P):return float(np.linalg.norm(np.diff(P,axis=0),axis=1).sum())

def run_comparison(data,out,alphas=(1,2,4),budgets=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[];paths={}
    budgets=budgets or dict(max_states=300000,max_transitions=15000000,seconds=90.)
    for key,item in data['instances'].items():
        A=item['A'];B=item['B_xyz'];ea=np.linalg.norm(np.diff(A,axis=0),axis=1);eb=np.linalg.norm(np.diff(B,axis=0),axis=1)
        e=.5*(ea.mean()+eb.mean());original=discrete_frechet(A,B);d3=float(np.ceil(original))
        for alpha in alphas:
            d=float(alpha*e);case=f'{key}@{alpha}';start=time.perf_counter()
            ia=shortest_independent(A,d);ib=shortest_independent(B,d);elapsed=time.perf_counter()-start
            gcs=dict(indices_A=ia,indices_B=ib,status='optimal_global_GCS',k=max(len(ia),len(ib)),seconds=elapsed)
            results={'GCS':gcs}
            for method in ('CPS-3F','CPS-2F'):
                begin=time.perf_counter()
                try:results[method]=solve_cps(A,B,d,d,d3,method,**budgets)
                except ResourceLimit as error:results[method]=dict(status='resource_limit',reason=str(error),seconds=time.perf_counter()-begin)
            for method,result in results.items():
                row=dict(case=case,pair=item['pair'],B=item['B'],tables=item['tables'],w=item['w'],role=item['role'],
                    alpha=alpha,method=method,nA=len(A),nB=len(B),arc_A=arc(A),arc_B=arc(B),
                    edge_mean=float(e),edge_min=float(min(ea.min(),eb.min())),edge_max=float(max(ea.max(),eb.max())),
                    edge_mean_before=item['edge_mean_before'],edge_growth=float(e/item['edge_mean_before']),
                    delta1=d,delta2=d,delta3=d3,dDF_original=float(original),endpoint_floor=float(max(np.linalg.norm(A[0]-B[0]),np.linalg.norm(A[-1]-B[-1]))),
                    status=result['status'],seconds=result.get('seconds'),reason=result.get('reason',''),
                    states=result.get('states'),transitions=result.get('transitions'),
                    auxiliary_A=result.get('auxiliary_A'),auxiliary_B=result.get('auxiliary_B'))
                if 'indices_A' in result:
                    a=result['indices_A'];b=result['indices_B'];sa=A[a];sb=B[b];coupling,dist=discrete_coupling(sa,sb)
                    ca=continuous_distance(A,sa);cb=continuous_distance(B,sb)
                    maxlink=float(np.linalg.norm(sa-sb,axis=1).max()) if len(sa)==len(sb) else None
                    row.update(kA=len(a),kB=len(b),k=max(len(a),len(b)),arc_simplified_A=arc(sa),arc_simplified_B=arc(sb),
                        lenratio_A=arc(sa)/arc(A),lenratio_B=arc(sb)/arc(B),dDF_simplified=dist,
                        dDF_change=dist-original,discrete_fidelity_A=float(discrete_frechet(A,sa)),discrete_fidelity_B=float(discrete_frechet(B,sb)),
                        continuous_fidelity_A_lower=ca[0],continuous_fidelity_A_upper=ca[1],
                        continuous_fidelity_B_lower=cb[0],continuous_fidelity_B_upper=cb[1],
                        coupling_steps=len(coupling),coupling_max_link=dist,pair_bound_verified=bool(dist<=d3+1e-9),
                        equal_size=len(a)==len(b),indexwise_max_link=maxlink,
                        vertices_over_GCS=max(len(a),len(b))-gcs['k'])
                    paths[f'{case}@{method}']=dict(indices_A=a.tolist(),indices_B=b.tolist(),coupling=coupling.tolist())
                rows.append(row)
            print(f"{item['B']} w={item['w']} alpha={alpha}: "+', '.join(f"{m}={r.get('k',r['status'])}" for m,r in results.items()),flush=True)
            pd.DataFrame(rows).to_csv(out/'comparison.csv',index=False)
            (out/'paths.json').write_text(json.dumps(paths,indent=2))
    frame=pd.DataFrame(rows)
    for name in ('inputs','audit','paper','qualification','residues'):data[name].to_csv(out/f'{name}.csv',index=False)
    return frame,paths
