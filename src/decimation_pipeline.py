from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from decimation_algorithms import (lockstep_decimate,edge_stats,arc_length,shortcut_dag,shortest_path,joint_local,discrete_frechet)

def run_sweep(pairs,out,ws=(1,2,4,8,16),alphas=(1,2,4),verbose=True):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[];all_paths={};maps=[]
    for pairid,pair in pairs.items():
        A=pair['A'];B=pair['B'];before=edge_stats(A,B);df_before=float(discrete_frechet(A,B))
        for w in ws:
            a,b,indices=lockstep_decimate(A,B,w)
            after=edge_stats(a,b);df_after=float(discrete_frechet(a,b));delta3=float(np.ceil(df_after))
            stub=bool((len(A)-1)%w!=0)
            for alpha in alphas:
                delta1=alpha*after['mean'];t=time.perf_counter()
                ea=shortcut_dag(a,delta1);eb=shortcut_dag(b,delta1);graph_time=time.perf_counter()-t
                t=time.perf_counter();sa=shortest_path(ea);sb=shortest_path(eb);sep_time=time.perf_counter()-t
                equal=len(sa)==len(sb);sep_max=float(np.linalg.norm(a[sa]-b[sb],axis=1).max()) if equal else None
                t=time.perf_counter();cps=joint_local(a,b,ea,eb,delta3);cps_time=time.perf_counter()-t
                row=dict(pair=pairid,B_printed=pair['B_printed'],B_verified=pair['B_verified'],regime=pair['regime'],
                    w=w,stride_role='undecimated control / headline' if w==1 else ('headline' if w<=4 else 'RUNTIME ONLY'),
                    nA_full=pair['nA_full'],nB_full=pair['nB_full'],source_index_start=pair['source_index_start'],source_index_end=pair['source_index_end'],
                    nA_before=len(A),nB_before=len(B),edge_mean_before=before['mean'],dF_before=df_before,
                    nA_after=len(a),nB_after=len(b),edge_mean_after=after['mean'],edge_min_after=after['minimum'],edge_max_after=after['maximum'],
                    edge_mean_A_after=after['mean_A'],edge_mean_B_after=after['mean_B'],edge_last_A_after=after['last_A'],edge_last_B_after=after['last_B'],
                    appended_stub=stub,last_stride=int(indices[-1]-indices[-2]),keep_frac=len(a)/len(A),dF_after=df_after,dF_change=df_after-df_before,
                    edge_growth_factor=after['mean']/before['mean'],growth_per_stride=after['mean']/before['mean']/w,coarse_trace=after['mean']>10,
                    dec_arc_ratio_A=arc_length(a)/arc_length(A),dec_arc_ratio_B=arc_length(b)/arc_length(B),
                    alpha=alpha,delta1=delta1,delta2=delta1,delta3=delta3,floor=cps['floor'],feasible=bool(delta3>=cps['floor']),
                    local_fidelity_floor=cps['local_floor'],solver_feasible=cps['ia'] is not None,
                    sep_kA=len(sa),sep_kB=len(sb),sep_equal=equal,sep_lenratio_A=arc_length(a[sa])/arc_length(a),
                    sep_lenratio_B=arc_length(b[sb])/arc_length(b),sep_dF=float(discrete_frechet(a[sa],b[sb])),
                    sep_max_link=sep_max,sep_ok3=(sep_max<=delta3+1e-9) if equal else None,
                    sep_coupling='coupled' if equal else "no coupling, |A'| != |B'|",
                    cps_k=None,cps_lenratio_A=None,cps_lenratio_B=None,cps_dF=None,cps_max_link=None,
                    cost_of_joining=None,cps_status=cps['status'],cps_verified=None,
                    graph_seconds=graph_time,sep_seconds=sep_time,cps_seconds=cps_time)
                path=dict(decimated_indices=indices.tolist(),sep_A=sa.tolist(),sep_B=sb.tolist(),cps_A=None,cps_B=None)
                if cps['ia'] is not None:
                    ja=cps['ia'];jb=cps['ib'];ca=a[ja];cb=b[jb]
                    link=float(np.linalg.norm(ca-cb,axis=1).max());df=float(discrete_frechet(ca,cb))
                    assert link<=delta3+1e-9 and df<=delta3+1e-9
                    row.update(cps_k=cps['k'],cps_lenratio_A=arc_length(ca)/arc_length(a),cps_lenratio_B=arc_length(cb)/arc_length(b),
                        cps_dF=df,cps_max_link=link,cost_of_joining=cps['k']-max(len(sa),len(sb)),cps_verified=True)
                    assert row['cost_of_joining']>=0
                    path.update(cps_A=ja.tolist(),cps_B=jb.tolist())
                rows.append(row);all_paths[f'{pairid}@{w}@{alpha}']=path
                for side in ('A','B'):
                    fullmap=pair[f'map_{side}']
                    for stage,ii in [('raw_paired_input',np.arange(len(A))),('decimated',indices),
                        ('separate',indices[path[f'sep_{side}']]),('joint_CPS',None if path[f'cps_{side}'] is None else indices[path[f'cps_{side}']])]:
                        if ii is None:continue
                        mapping=fullmap.iloc[ii].copy();mapping.insert(0,'stage',stage);mapping.insert(0,'side',side)
                        mapping.insert(0,'alpha',alpha);mapping.insert(0,'w',w);mapping.insert(0,'pair',pairid)
                        mapping['stage_index']=np.arange(len(ii));maps.append(mapping)
                if verbose:print(f'{pair["B_verified"]} w={w} alpha={alpha}: {len(A)}->{len(a)}->{len(sa)}/{cps.get("k")}',flush=True)
        pd.DataFrame(rows).to_csv(out/'pipeline.csv',index=False)
        (out/'paths.json').write_text(json.dumps(all_paths,indent=2))
    frame=pd.DataFrame(rows)
    pd.concat(maps,ignore_index=True).to_csv(out/'stage_residue_mapping.csv',index=False)
    frame.to_csv(out/'pipeline.csv',index=False)
    (out/'paths.json').write_text(json.dumps(all_paths,indent=2))
    for scope,mask in [('headline',frame.w<=4),('runtime_with_control',(frame.w==1)|(frame.w>=8))]:
        subset=frame[mask & frame.solver_feasible]
        subset.groupby(['w','alpha']).agg(pairs=('pair','size'),sep_coupled=('sep_equal','sum'),
            mean_joining_cost=('cost_of_joining','mean'),mean_edge_growth=('edge_growth_factor','mean'),
            mean_cps_seconds=('cps_seconds','mean')).to_csv(out/f'aggregates_{scope}.csv')
    return frame,all_paths

if __name__=='__main__':
    from decimation_data import prepare_decimation_pairs
    root=Path(__file__).resolve().parents[1];pairs,_,_=prepare_decimation_pairs(root)
    frame,_=run_sweep(pairs,root/'output/decimation')
    print(frame.groupby(['w','alpha']).size().to_string())
