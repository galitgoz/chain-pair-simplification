"""Postprocess saved simplifications: compression, auxiliary counts and contacts."""
import numpy as np
import pandas as pd
from cps_paper_algorithms import extended_curve


def contact_metrics(P, ids, sequence, cutoff=8., separation=5):
    """Compare contacts on the SAME input-vertex grid, excluding local residues.

    Lift a shortcut to discarded input vertices by original subcurve arc fraction.
    This is a fixed evaluation correspondence, not a Frechet matching or atom model.
    """
    P=np.asarray(P,float);ids=np.asarray(ids,int);sequence=np.asarray(sequence,int)
    lifted=np.empty_like(P)
    for a,b in zip(ids[:-1],ids[1:]):
        lengths=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(P[a:b+1],axis=0),axis=1))]
        t=lengths/lengths[-1] if lengths[-1]>0 else np.linspace(0,1,b-a+1)
        lifted[a:b+1]=P[a]+t[:,None]*(P[b]-P[a])
    if len(ids)==1:lifted[:]=P[ids[0]]
    mask=np.triu(np.abs(sequence[:,None]-sequence[None,:])>=separation,1)
    before=np.linalg.norm(P[:,None]-P[None,:],axis=2)[mask]<=cutoff
    after=np.linalg.norm(lifted[:,None]-lifted[None,:],axis=2)[mask]<=cutoff
    tp=int(np.sum(before&after));fp=int(np.sum(~before&after));fn=int(np.sum(before&~after))
    return dict(contact_TP=tp,contact_FP=fp,contact_FN=fn,contacts_input=int(before.sum()),
        contacts_output=int(after.sum()),contact_candidates=int(mask.sum()),
        contact_precision=tp/(tp+fp) if tp+fp else np.nan,
        contact_recall=tp/(tp+fn) if tp+fn else np.nan,
        contact_F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan)


def enrich_results(results,data,paths,out,cutoff=8.,separation=5):
    rows=[];aux_cache={}
    for row in results.to_dict('records'):
        key=f"{row['pair']}@{row['w']}";item=data['instances'][key]
        row['A']=row['pair'].split('__')[0]
        row.update(contact_cutoff=cutoff,contact_sequence_separation=separation)
        row['auxiliary_A']=row['auxiliary_B']=0
        row['extended_A']=row['nA'];row['extended_B']=row['nB']
        row['auxiliary_graph_used']=row['method']=='CPS-2F' and row['status']=='optimal_configuration_graph'
        if row['method']=='CPS-2F':
            for side,curve,delta in [('A',item['A'],row['delta1']),('B',item['B_xyz'],row['delta2'])]:
                ck=(key,side,delta)
                if ck not in aux_cache:aux_cache[ck]=len(extended_curve(curve,delta)[0])-len(curve)
                row[f'auxiliary_{side}']=aux_cache[ck]
                row[f'extended_{side}']=len(curve)+aux_cache[ck]
        row['auxiliary_total']=row['auxiliary_A']+row['auxiliary_B']
        p=paths.get(f"{row['case']}@{row['method']}")
        if p is not None:
            counts=[]
            for side,P in [('A',item['A']),('B',item['B_xyz'])]:
                ids=p[f'indices_{side}'];n=len(P);k=len(ids)
                row[f'compression_factor_{side}']=n/k
                row[f'removed_pct_{side}']=100*(1-k/n)
                metrics=contact_metrics(P,ids,item[f'map_{side}'].label_seq_id.to_numpy(),cutoff,separation)
                row.update({f'{name}_{side}':value for name,value in metrics.items()})
                counts.append(metrics)
            tp=sum(m['contact_TP'] for m in counts);fp=sum(m['contact_FP'] for m in counts);fn=sum(m['contact_FN'] for m in counts)
            row['contact_F1_pooled']=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan
            row['removed_pct_total']=100*(1-(row['kA']+row['kB'])/(row['nA']+row['nB']))
            row['max_continuous_error']=max(row['continuous_fidelity_A_upper'],row['continuous_fidelity_B_upper'])
        rows.append(row)
    result=pd.DataFrame(rows)
    result.to_csv(out/'comparison.csv',index=False)
    return result


def paired_comparison(results,out):
    keys=['case','pair','w','alpha','delta1','delta2','delta3','nA','nB']
    cols=['status','kA','kB','k','removed_pct_total','contact_F1_pooled','max_continuous_error','dDF_simplified']
    a=results[results.method=='CPS-2F'][keys+cols]
    b=results[results.method=='CPS-3F'][keys+cols]
    paired=a.merge(b,on=keys,suffixes=('_2F','_3F'),validate='one_to_one')
    paired['both_solved']=paired.status_2F.str.startswith('optimal')&paired.status_3F.str.startswith('optimal')
    paired['k_saved_by_2F']=paired.k_3F-paired.k_2F
    paired['extra_removal_pp_2F']=paired.removed_pct_total_2F-paired.removed_pct_total_3F
    paired['contact_F1_gain_2F']=paired.contact_F1_pooled_2F-paired.contact_F1_pooled_3F
    paired['continuous_error_change_2F']=paired.max_continuous_error_2F-paired.max_continuous_error_3F
    paired.to_csv(out/'cps2f_vs_cps3f.csv',index=False)
    return paired
