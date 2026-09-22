"""One complete, auditable pipeline row per aligned pair and fidelity threshold."""
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from curve_algorithms import local_shortcut_path,discrete_frechet,segment_frechet_decision
from small_pair_solvers import synchronous_cps,arc_length

def run_pipeline(pairs,output,thresholds=(2.,4.,8.),delta3_mode='exact'):
    if delta3_mode not in ['exact','ceil']: raise ValueError('delta3_mode must be exact or ceil')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    rows=[];paths={}; selected=[]
    for pairid,pair in pairs.items():
        A=pair['A'];B=pair['B']
        original=float(discrete_frechet(A,B)); paper_threshold=float(np.ceil(original))
        d3=original if delta3_mode=='exact' else paper_threshold
        la=arc_length(A);lb=arc_length(B)
        for delta in thresholds:
            t=time.perf_counter(); ia=local_shortcut_path(A,delta);ib=local_shortcut_path(B,delta);separate_time=time.perf_counter()-t
            aa=A[ia];bb=B[ib]; sa=arc_length(aa);sb=arc_length(bb)
            coupled=len(ia)==len(ib)
            index_max=float(np.linalg.norm(aa-bb,axis=1).max()) if coupled else None
            separate_df=float(discrete_frechet(aa,bb))
            t=time.perf_counter();joint=synchronous_cps(A,B,delta,delta,d3);cps_time=time.perf_counter()-t
            row=dict(pair=pairid,A_printed=pair['A_printed'],B_printed=pair['B_printed'],delta1=delta,delta2=delta,
                delta3_mode=delta3_mode,delta3=d3,paper_ceil_delta3=paper_threshold,
                original_n_A=len(A),original_n_B=len(B),original_arc_A=la,original_arc_B=lb,original_dDF=original,
                separate_k_A=len(ia),separate_k_B=len(ib),separate_arc_A=sa,separate_arc_B=sb,
                separate_arc_ratio_A=sa/la if la else np.nan,separate_arc_ratio_B=sb/lb if lb else np.nan,
                separate_dDF=separate_df,separate_dDF_minus_original=separate_df-original,
                separate_further_apart=separate_df>original+1e-9,
                separate_coupled=coupled,separate_index_max=index_max,
                separate_index_le_delta3=(index_max<=d3+1e-9) if coupled else None,
                separate_coupling_status='index-wise coupling' if coupled else "no coupling, |A'| != |B'|",
                feasibility_floor=joint['floor'],status=joint['status'],included_in_aggregates=joint['status']=='feasible',
                CPS_k=None,CPS_arc_A=None,CPS_arc_B=None,CPS_arc_ratio_A=None,CPS_arc_ratio_B=None,
                CPS_index_max=None,CPS_dDF=None,CPS_solver_verified=None,CPS_fidelity_A=None,CPS_fidelity_B=None,
                cost_of_joining=None,separate_seconds=separate_time,CPS_seconds=cps_time)
            path=dict(separate_A=ia.tolist(),separate_B=ib.tolist(),CPS_A=None,CPS_B=None)
            if joint['status']=='feasible':
                ja=joint['indices_A'];jb=joint['indices_B'];ca=A[ja];cb=B[jb]
                da=float(discrete_frechet(A,ca));db=float(discrete_frechet(B,cb))
                paired=float(np.linalg.norm(ca-cb,axis=1).max());df=float(discrete_frechet(ca,cb))
                verified=paired<=d3+1e-9 and df<=d3+1e-9 and da<=delta+1e-9 and db<=delta+1e-9
                assert verified
                row.update(CPS_k=joint['k'],CPS_arc_A=arc_length(ca),CPS_arc_B=arc_length(cb),
                    CPS_arc_ratio_A=arc_length(ca)/la if la else np.nan,CPS_arc_ratio_B=arc_length(cb)/lb if lb else np.nan,
                    CPS_index_max=paired,CPS_dDF=df,CPS_solver_verified=verified,CPS_fidelity_A=da,CPS_fidelity_B=db,
                    cost_of_joining=joint['k']-max(len(ia),len(ib)))
                path.update(CPS_A=ja.tolist(),CPS_B=jb.tolist())
            for side,P,ids in [('A',A,ia),('B',B,ib)]:
                assert all(segment_frechet_decision(P[i:j+1],P[i],P[j],delta) for i,j in zip(ids,ids[1:]))
            for method in ['separate','CPS']:
                for side in ['A','B']:
                    ids=path[f'{method}_{side}']
                    if ids is not None:
                        mapping=pair[f'map_{side}'].iloc[ids].copy()
                        mapping.insert(0,'method',method);mapping.insert(0,'side',side)
                        mapping.insert(0,'delta1',delta);mapping.insert(0,'pair',pairid)
                        mapping['simplified_index']=np.arange(len(mapping));selected.append(mapping)
            rows.append(row); paths[f'{pairid}@{delta:g}']=path
            print(f'{pair["B_key"]} delta={delta:g}: independent={len(ia)}/{len(ib)}, CPS={row["CPS_k"]}, floor={joint["floor"]:.3f}, delta3={d3:.3f}',flush=True)
    frame=pd.DataFrame(rows)
    frame.to_csv(output/'pipeline.csv',index=False)
    (output/'paths.json').write_text(json.dumps(paths,indent=2))
    pd.concat(selected,ignore_index=True).to_csv(output/'selected_residues.csv',index=False)
    feasible=frame[frame.included_in_aggregates]
    aggregate=feasible.groupby('delta1').agg(feasible_pairs=('pair','size'),
        independent_further_apart=('separate_further_apart','sum'),mean_independent_dDF_change=('separate_dDF_minus_original','mean'),
        mean_cost_of_joining=('cost_of_joining','mean'),verified_CPS=('CPS_solver_verified','sum'))
    aggregate.to_csv(output/'aggregates_feasible_only.csv')
    return frame,paths,aggregate

if __name__=='__main__':
    from small_pair_data import prepare_small_pairs
    root=Path(__file__).resolve().parents[1]
    _,pairs,_=prepare_small_pairs(root)
    frame,_,aggregate=run_pipeline(pairs,root/'output/small_pairs')
    print(aggregate.to_string())
