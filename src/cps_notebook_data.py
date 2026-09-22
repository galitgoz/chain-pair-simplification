"""Load supplied decimation unchanged; explicitly qualify gaps and supplement Table 3."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from protein_data import paper_experiments

def load_inputs(root,ws=(16,),include_table3=True,scope='supplied'):
    root=Path(root);source=root/'data/decimated protein backbones';prior=root/'output/decimation'
    curves=pd.read_csv(source/'decimated_curves.csv')
    chain_summary=pd.read_csv(source/'decimated_chains.csv')
    supplied_pairs=pd.read_csv(source/'decimated_pairs.csv')
    audit=pd.read_csv(prior/'identifier_audit.csv');qualification=pd.read_csv(prior/'input_qualification.csv')
    for row in audit.itertuples():
        path=root/'data/pdb_cache'/f'{row.verified.split(".")[0]}.cif'
        assert hashlib.sha256(path.read_bytes()).hexdigest()==row.sha256,'PDB cache changed: re-audit required'
    for (_,frame,w),g in curves.groupby(['chain','frame','w']):
        g=g.sort_values('t');idx=np.arange(0,325,int(w))
        if idx[-1]!=324:idx=np.r_[idx,324]
        assert np.array_equal(g.source_row,idx)
        assert np.isfinite(g[['x','y','z']]).all().all()
    if scope not in ('supplied','qualified'): raise ValueError('scope must be supplied or qualified')
    if not ws or any(w not in set(curves.w) for w in ws): raise ValueError('Choose supplied stride levels')
    lookup=audit.set_index('verified');instances={};rows=[];residues=[]
    available=set(curves.chain)
    for q in qualification.itertuples():
        in_csv=q.B_verified in available
        if not in_csv and not include_table3:continue
        folder=prior/'inputs'/q.pair
        cache=np.load(folder/'aligned.npz')
        maps={side:pd.read_csv(folder/f'{side}.residues.csv') for side in ['A','B']}
        if scope=='supplied':
            maps={side:pd.read_csv(root/'data/prepared'/f'{chain}.residues.csv')
                  for side,chain in [('A',q.A_verified),('B',q.B_verified)]}
        fit=json.loads((folder/'prepared.json').read_text())['alignment']
        R=np.array(fit['rotation_B_to_A']);t=np.array(fit['translation_B_to_A'])
        edge_before=.5*sum(np.linalg.norm(np.diff(maps[s][['x','y','z']].to_numpy(float),axis=0),axis=1).mean() for s in ('A','B'))
        for w in ws:
            if scope=='supplied':
                indices={};coordinates={}
                for side,chain in [('A',q.A_verified),('B',q.B_verified)]:
                    if in_csv:
                        g=curves[(curves.chain==chain)&(curves.frame=='aligned')&(curves.w==w)].sort_values('t')
                        indices[side]=g.source_row.to_numpy(int)
                        coordinates[side]=g[['x','y','z']].to_numpy(float)
                    else:
                        ix=np.arange(0,len(maps[side]),w)
                        if ix[-1]!=len(maps[side])-1: ix=np.r_[ix,len(maps[side])-1]
                        indices[side]=maps[side].curve_index.to_numpy(int)[ix]
                        coordinates[side]=maps[side][['x','y','z']].to_numpy(float)[ix]
                        if side=='B': coordinates[side]=coordinates[side]@R.T+t
                A=coordinates['A'];B=coordinates['B'];idx=indices['A'];dropped=0
                origin=('supplied aligned CSV; all vertices retained' if in_csv else
                        'Table 3 supplement: full audited chains; independent uniform strides; cached rigid US-align transform')
            elif in_csv:
                ga=curves[(curves.chain==q.A_verified)&(curves.frame=='deposited')&(curves.w==w)].sort_values('t')
                gb=curves[(curves.chain==q.B_verified)&(curves.frame=='deposited')&(curves.w==w)].sort_values('t')
                assert np.array_equal(ga.source_row,gb.source_row)
                keep=(ga.source_row>=q.source_index_start)&(ga.source_row<=q.source_index_end)
                dropped=int((~keep).sum());idx=ga.source_row[keep].to_numpy(int)
                A=ga.loc[keep,['x','y','z']].to_numpy(float)
                B=gb.loc[keep.to_numpy(),['x','y','z']].to_numpy(float)@R.T+t
                origin='supplied decimated CSV; no second decimation'
            else:
                relative=np.arange(0,len(cache['A']),w)
                if relative[-1]!=len(cache['A'])-1:relative=np.r_[relative,len(cache['A'])-1]
                idx=cache['source_indices'][relative]
                A=cache['A'][relative];B=cache['B'][relative];dropped=0
                origin='Table 3 supplement: cached audited continuous interval; exact uniform stride'
            if len(A)<2:raise ValueError('Insufficient vertices after explicit gap qualification')
            am=maps['A'].set_index('curve_index').loc[idx].reset_index()
            bm=maps['B'].set_index('curve_index').loc[indices['B'] if scope=='supplied' else idx].reset_index()
            if scope=='qualified':
                assert am.fragment_id.nunique()==bm.fragment_id.nunique()==1
                assert np.max(np.abs(A-am[['x','y','z']].to_numpy(float)))<.002
                assert np.max(np.abs(B-(bm[['x','y','z']].to_numpy(float)@R.T+t)))<.002
            key=f'{q.pair}@{w}'
            row=dict(instance=key,pair=q.pair,B=q.B_verified,A_printed=q.A_printed,B_printed=q.B_printed,
                tables='1,2' if in_csv else '3',w=w,nA=len(A),nB=len(B),origin=origin,
                dropped_decimated_vertices_per_chain=dropped,source_first=int(idx[0]),source_last=int(idx[-1]),
                nA_full=q.nA_full,nB_full=q.nB_full,qualified_raw_n=q.nA_before,
                edge_mean_before=float(edge_before),
                full_chain_endpoints_preserved=True if scope=='supplied' else bool(q.full_chain_endpoints_preserved),
                input_scope=scope,fragments_A=int(am.fragment_id.nunique()),fragments_B=int(bm.fragment_id.nunique()),
                regime=q.regime,TM_score=q.TM_score,
                role='runtime/coarse trace only' if w>=8 else ('undecimated control' if w==1 else 'headline'))
            instances[key]=dict(**row,A=A,B_xyz=B,map_A=am,map_B=bm)
            rows.append(row)
            for side,mp in [('A',am),('B',bm)]:
                mp=mp.copy();mp.insert(0,'input_index',np.arange(len(mp)));mp.insert(0,'side',side);mp.insert(0,'instance',key);residues.append(mp)
    return dict(curves=curves,chain_summary=chain_summary,supplied_pairs=supplied_pairs,audit=audit,
                paper=paper_experiments(),qualification=qualification,instances=instances,
                inputs=pd.DataFrame(rows),residues=pd.concat(residues,ignore_index=True))
