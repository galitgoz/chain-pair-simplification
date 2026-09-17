"""Auditable PDB/mmCIF to fragment preparation. No interpolation across gaps."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

SELECTIONS = [
 ('107j.a','1o7j','A',325,'Probable letter-o/zero typo; local coordinates match 1O7J A; literal 107J returns 404.'),
 ('1hfj.c','1hfj','C',325,'Author chain C; label chain B.'),
 ('1qd1.b','1qd1','B',325,'Author chain B; author numbering starts at 2002.'),
 ('1toh','1toh','A',325,'Chain omitted in paper; A is the only protein chain and matches the local file.'),
 ('4eca.c','4eca','C',325,'Local file omits a deposited polymer residue; retain raw subset and split there.'),
 ('1d9q.d','1d9q','D',297,'Reported count differs from both local raw and current deposited C-alpha count.'),
 ('4cea.b','4eca','B',325,'Probable transposition; 4CEA is HIV integrase (B: 150 C-alpha), local file matches 4ECA B.'),
 ('4cea.d','4eca','D',325,'Probable transposition; 4CEA has no protein chain D, local file matches 4ECA D.'),
 ('3ntx.a','3ntx','A',322,'No supplied raw file; current deposition used without trimming.'),
 ('1wls.a','1wls','A',316,'No supplied raw file; current deposition used without trimming.'),
 ('2eq5.a','2eq5','A',215,'No supplied raw file; current deposition used without trimming.'),
 ('2zsk.a','2zsk','A',219,'No supplied raw file; current deposition used without trimming.'),
 ('1zq1.a','1zq1','A',418,'No supplied raw file; current deposition used without trimming.'),
 ('3jq0.a','3jq0','A',457,'No supplied raw file; current deposition used without trimming.'),
 ('2fep.a','2fep','A',273,'No supplied raw file; current deposition used without trimming.'),
]

def paper_experiments():
    similar = [s[0] for s in SELECTIONS[1:8]]
    t1 = [(4,4,z,a,b) for z,a,b in zip([1,21,21,6,20,5,5],[109,126,149,111,130,111,113],[83,82,84,83,82,82,84])]
    t2 = [(x,x,z,a,b) for x,z,a,b in zip([12,15,16,12,15,12,12],[1,12,13,3,13,3,3],[26,21,22,27,24,26,32],[15,11,11,16,12,15,16])]
    t3 = [(10,10,5,39,25),(15,13,6,22,14),(8,6,19,58,32),(12,8,17,38,19),(10,12,19,45,23),(12,12,26,70,36),(12,12,10,11,6)]
    rows=[]
    for table,names,params in [(1,similar,t1),(2,similar,t2),(3,[s[0] for s in SELECTIONS[8:]],t3)]:
        for name,(d1,d2,d3,old,opt) in zip(names,params):
            rows.append(dict(table=table,A_printed='107j.a',B_printed=name,delta1=d1,delta2=d2,delta3=d3,paper_CPS3F_plus=old,paper_CPS3F=opt))
    return pd.DataFrame(rows)

def read_raw(path):
    rows=[row.strip() for row in path.read_text().strip().split(';') if row.strip()]
    xyz=np.array([[float(x) for x in row.split(',')] for row in rows])
    if xyz.ndim!=2 or xyz.shape[1]!=3 or not np.isfinite(xyz).all():
        raise ValueError(f'Invalid raw coordinates: {path}')
    return xyz

def atom_table(path):
    cif=MMCIF2Dict(str(path))
    prefix='_atom_site.'
    frame=pd.DataFrame({k[len(prefix):]:v for k,v in cif.items() if k.startswith(prefix)})
    for col in ['Cartn_x','Cartn_y','Cartn_z','occupancy']:
        frame[col]=pd.to_numeric(frame[col])
    return cif,frame

def prepare(root, ca_break=4.5, peptide_break=2.0):
    root=Path(root); chains={}; audits=[]; gaps=[]; all_missing=[]
    for printed,pid,author,reported,note in SELECTIONS:
        key=f'{pid}.{author.lower()}'
        path=root/'data'/'pdb_cache'/f'{pid}.cif'
        cif,atoms=atom_table(path)
        polymer=atoms[(atoms.auth_asym_id==author)&(~atoms.label_seq_id.isin(['.','?']))].copy()
        if polymer.empty: raise ValueError(f'No polymer author chain {key}')
        models=sorted(polymer.pdbx_PDB_model_num.unique(),key=int)
        model=models[0]
        polymer=polymer[polymer.pdbx_PDB_model_num==model]
        labels=polymer.label_asym_id.unique()
        if len(labels)!=1: raise ValueError(f'Ambiguous author/label chain mapping: {key}')
        entity=polymer.label_entity_id.iloc[0]
        seq=[int(s) for e,s in zip(cif['_entity_poly_seq.entity_id'],cif['_entity_poly_seq.num']) if e==entity]
        ca=polymer[polymer.label_atom_id=='CA'].copy()
        ca['seq']=ca.label_seq_id.astype(int)
        ca['alt_rank']=ca.label_alt_id.map(lambda a: '0' if a in ['.','?'] else ('1' if a=='A' else '2'+a))
        ca=ca.sort_values(['seq','occupancy','alt_rank'],ascending=[True,False,True])
        canonical=ca.drop_duplicates('seq').sort_values('seq').copy()
        rawpath=root/'data'/'raw protein backbones'/f'{key}.raw'
        match_error=None
        if rawpath.exists():
            xyz=read_raw(rawpath)
            # Match all alternate C-alpha atoms, never infer residue numbers from row counts.
            deposited=ca[['Cartn_x','Cartn_y','Cartn_z']].to_numpy(float)
            distances=np.linalg.norm(xyz[:,None,:]-deposited[None,:,:],axis=2)
            best=distances.argmin(axis=1)
            match_error=float(distances[np.arange(len(xyz)),best].max())
            if match_error>0.002: raise ValueError(f'{key}: raw cannot be mapped exactly (max {match_error} A)')
            for drow in distances:
                candidates=ca.loc[drow<=0.002,'seq'].unique()
                if len(candidates)!=1: raise ValueError(f'{key}: ambiguous residue coordinate match')
            chosen=ca.iloc[best].copy().reset_index(drop=True)
            if not np.all(np.diff(chosen.seq)>0): raise ValueError(f'{key}: repeated/out-of-order raw residues')
            source='local raw, exact coordinate mapping'
        else:
            chosen=canonical.reset_index(drop=True)
            xyz=chosen[['Cartn_x','Cartn_y','Cartn_z']].to_numpy(float)
            source='current mmCIF, selected C-alpha conformer'
        canonical_by_seq=canonical.set_index('seq')
        residue_rows=[]
        for index,row in chosen.iterrows():
            residue_rows.append(dict(curve_index=int(index),raw_index=int(index) if rawpath.exists() else None,
                pdb_id=pid,model=int(model),author_chain=author,label_chain=row.label_asym_id,
                label_seq_id=int(row.seq),author_seq_id=row.auth_seq_id,insertion_code='' if row.pdbx_PDB_ins_code in ['.','?'] else row.pdbx_PDB_ins_code,
                residue_name=row.label_comp_id,atom_id=row.id,altloc=row.label_alt_id,occupancy=float(row.occupancy),
                canonical_altloc=canonical_by_seq.loc[row.seq,'label_alt_id'],
                x=float(xyz[index,0]),y=float(xyz[index,1]),z=float(xyz[index,2])))
        mapping=pd.DataFrame(residue_rows)
        def backbone_atom(seqid,name,alt):
            rows=polymer[(polymer.label_seq_id==str(seqid))&(polymer.label_atom_id==name)]
            rows=rows[rows.label_alt_id.isin(['.','?',alt])].sort_values('occupancy',ascending=False)
            if rows.empty: return None
            return rows.iloc[0][['Cartn_x','Cartn_y','Cartn_z']].to_numpy(float)
        starts=[0]
        for i in range(1,len(mapping)):
            prev,curr=mapping.iloc[i-1],mapping.iloc[i]
            reasons=[]
            if curr.label_seq_id!=prev.label_seq_id+1: reasons.append('absent selected polymer residues')
            ca_distance=float(np.linalg.norm(xyz[i]-xyz[i-1]))
            if ca_distance>ca_break: reasons.append('C-alpha distance exceeds threshold')
            c=backbone_atom(prev.label_seq_id,'C',prev.altloc)
            n=backbone_atom(curr.label_seq_id,'N',curr.altloc)
            peptide_distance=None if c is None or n is None else float(np.linalg.norm(c-n))
            if peptide_distance is None: reasons.append('peptide connectivity unverifiable: C or N absent for chosen conformer')
            elif peptide_distance>peptide_break: reasons.append('C-N distance exceeds threshold')
            if reasons:
                starts.append(i)
                gaps.append(dict(chain=key,before_index=i-1,after_index=i,
                    before_residue=f'{prev.author_seq_id}{prev.insertion_code}',after_residue=f'{curr.author_seq_id}{curr.insertion_code}',
                    missing_label_seq_ids=list(range(prev.label_seq_id+1,curr.label_seq_id)),ca_distance=ca_distance,
                    peptide_distance=peptide_distance,reason='; '.join(reasons)))
        fragments=[]
        for f,(start,end) in enumerate(zip(starts,starts[1:]+[len(mapping)])):
            ix=np.arange(start,end)
            mapping.loc[ix,'fragment_id']=f
            mapping.loc[ix,'fragment_index']=np.arange(len(ix))
            fragments.append(dict(id=f'{key}/f{f:02d}',chain=key,xyz=xyz[ix].copy(),chain_indices=ix,
                                  residues=mapping.iloc[ix].copy()))
        deposited_seq=set(canonical.seq)
        observed_seq=set(polymer.label_seq_id.astype(int))
        selected_seq=set(chosen.seq)
        for s in seq:
            reason=('unresolved residue (no deposited atoms)' if s not in observed_seq else
                    'residue present but C-alpha missing' if s not in deposited_seq else
                    'deposited C-alpha omitted by local raw file' if s not in selected_seq else None)
            if reason:
                scheme_prefix='_pdbx_poly_seq_scheme.'
                scheme_keys=[k for k in cif if k.startswith(scheme_prefix)]
                scheme_rows=[dict(zip([k[len(scheme_prefix):] for k in scheme_keys],r)) for r in zip(*(cif[k] for k in scheme_keys))]
                scheme=next((r for r in scheme_rows if r['asym_id']==labels[0] and int(r['seq_id'])==s),{})
                all_missing.append(dict(chain=key,label_seq_id=s,author_seq_id=scheme.get('auth_seq_num','?'),
                    pdb_seq_num=scheme.get('pdb_seq_num','?'),insertion_code=scheme.get('pdb_ins_code','?'),
                    residue_name=scheme.get('mon_id','?'),reason=reason))
        audits.append(dict(printed=printed,verified=key,reported_length=reported,polymer_sequence_length=len(seq),
            deposited_ca_length=len(canonical),prepared_ca_length=len(xyz),local_raw_length=len(xyz) if rawpath.exists() else None,
            length_matches_paper=len(xyz)==reported,author_chain=author,label_chain=labels[0],model=int(model),models_available=','.join(models),
            alternate_ca_residues=int((ca.groupby('seq').size()>1).sum()),
            selected_zero_occupancy_ca=int((chosen.occupancy<=0).sum()),
            selected_noncanonical_altloc=int((mapping.altloc!=mapping.canonical_altloc).sum()),fragments=len(fragments),
            raw_max_coordinate_error_A=match_error,source=source,interpretation=note,
            pdb_url=f'https://www.rcsb.org/structure/{pid.upper()}',sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        chains[key]=dict(xyz=xyz,mapping=mapping,fragments=fragments,canonical_xyz=canonical[['Cartn_x','Cartn_y','Cartn_z']].to_numpy(float))
    return chains,pd.DataFrame(audits),pd.DataFrame(gaps),pd.DataFrame(all_missing)
