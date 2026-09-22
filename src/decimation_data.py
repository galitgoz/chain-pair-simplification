"""Explicit data qualification for exact same-index, endpoint-preserving decimation."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from protein_data import prepare
from small_pair_data import structural_fit

def common_runs(Amap,Bmap):
    """Maximal original ROW-index intervals that are continuous in both chains.

    This is input qualification, never a decimation or geometry-based selection.
    """
    n=min(len(Amap),len(Bmap));starts=[0]
    for i in range(1,n):
        if Amap.iloc[i].fragment_id!=Amap.iloc[i-1].fragment_id or Bmap.iloc[i].fragment_id!=Bmap.iloc[i-1].fragment_id:
            starts.append(i)
    return list(zip(starts,starts[1:]+[n]))

def prepare_decimation_pairs(root,recompute_alignment=False):
    root=Path(root);out=root/'output'/'decimation';out.mkdir(parents=True,exist_ok=True)
    chains,audit,breaks,missing=prepare(root)
    audit.to_csv(out/'identifier_audit.csv',index=False);breaks.to_csv(out/'chain_breaks.csv',index=False)
    missing.to_csv(out/'missing_residues.csv',index=False)
    exe=next((root/'tools'/'usalign').glob('*.exe'))
    akey='1o7j.a';pairs={};qualification=[];excluded=[]
    for bkey in chains:
        if bkey==akey:continue
        ca=chains[akey];cb=chains[bkey];ma=ca['mapping'];mb=cb['mapping']
        runs=common_runs(ma,mb)
        start,end=max(runs,key=lambda interval:interval[1]-interval[0])
        if end-start<2:raise ValueError(f'{bkey}: no common continuous input interval')
        indices=np.arange(start,end)
        A=ca['xyz'][indices].copy();B=cb['xyz'][indices].copy()
        am=ma.iloc[indices].copy().reset_index(drop=True);bm=mb.iloc[indices].copy().reset_index(drop=True)
        for mapping in (am,bm):mapping.insert(0,'paired_input_index',np.arange(len(mapping)))
        pairid=f'{akey}__{bkey}';folder=out/'inputs'/pairid;folder.mkdir(parents=True,exist_ok=True)
        signature=hashlib.sha256(A.tobytes()+B.tobytes()+am.to_csv(index=False).encode()+bm.to_csv(index=False).encode()+
            exe.read_bytes()+Path(__file__).read_bytes()).hexdigest()
        cache=folder/'prepared.json';stored=json.loads(cache.read_text()) if cache.exists() else {}
        if recompute_alignment or stored.get('signature')!=signature:
            aligned,alignment=structural_fit(A,B,am,bm,folder/'structural_alignment',exe)
            np.savez_compressed(folder/'aligned.npz',A=A,B=aligned,B_deposited=B,source_indices=indices)
            stored=dict(signature=signature,alignment=alignment)
            cache.write_text(json.dumps(stored,indent=2))
        else:
            data=np.load(folder/'aligned.npz');aligned=data['B'];alignment=stored['alignment']
        am.to_csv(folder/'A.residues.csv',index=False);bm.to_csv(folder/'B.residues.csv',index=False)
        score=alignment['tm_scores'][0]
        regime='similar' if score>=.5 else ('intermediate' if score>=.3 else 'dissimilar')
        ap=audit.set_index('verified').loc[akey,'printed'];bp=audit.set_index('verified').loc[bkey,'printed']
        row=dict(pair=pairid,A_printed=ap,B_printed=bp,A_verified=akey,B_verified=bkey,
            nA_full=len(ma),nB_full=len(mb),source_index_start=start,source_index_end=end-1,
            nA_before=len(A),nB_before=len(B),excluded_A=len(ma)-len(A),excluded_B=len(mb)-len(B),
            A_author_first=am.iloc[0].author_seq_id,A_author_last=am.iloc[-1].author_seq_id,
            B_author_first=bm.iloc[0].author_seq_id,B_author_last=bm.iloc[-1].author_seq_id,
            full_chain_endpoints_preserved=(start==0 and end==len(ma) and end==len(mb)),
            qualification='longest shared original-row interval continuous in both chains; earliest tie',
            TM_score=score,regime=regime)
        qualification.append(row)
        for side,mapping in [('A',ma),('B',mb)]:
            drop=mapping[~mapping.curve_index.isin(indices)].copy()
            if len(drop):
                drop.insert(0,'pair',pairid);drop.insert(1,'side',side)
                drop['exclusion_stage']='input qualification, NOT stride decimation';excluded.append(drop)
        pairs[pairid]=dict(**row,A=A,B=aligned,map_A=am,map_B=bm,alignment=alignment)
        print(f'{bkey}: full {len(ma)}/{len(mb)} -> qualified {len(A)}/{len(B)}, source rows {start}..{end-1}, {regime}',flush=True)
    qualification=pd.DataFrame(qualification);qualification.to_csv(out/'input_qualification.csv',index=False)
    if excluded:pd.concat(excluded).to_csv(out/'qualification_excluded_residues.csv',index=False)
    return pairs,audit,qualification

if __name__=='__main__':prepare_decimation_pairs(Path(__file__).resolve().parents[1])
