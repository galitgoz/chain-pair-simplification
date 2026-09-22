"""Contiguous, residue-traceable windows with US-align structural superposition."""
from pathlib import Path
import hashlib,json,re,subprocess
import numpy as np
import pandas as pd
from protein_data import SELECTIONS

def load_prepared(root):
    root=Path(root); chains={}
    for printed,pid,author,length,note in SELECTIONS:
        key=f'{pid}.{author.lower()}'
        mapping=pd.read_csv(root/'data'/'prepared'/f'{key}.residues.csv',keep_default_na=False)
        fragments=[]
        for f,group in mapping.groupby('fragment_id',sort=True):
            xyz=group[['x','y','z']].to_numpy(float)
            assert np.isfinite(xyz).all() and np.all(np.diff(group.label_seq_id)==1)
            fragments.append(dict(id=f'{key}/f{int(f):02d}',xyz=xyz,mapping=group.reset_index(drop=True)))
        chains[key]=dict(printed=printed,reported_n=length,fragments=fragments)
    return chains

def write_ca_pdb(path,P,mapping):
    lines=[]
    for i,(point,name) in enumerate(zip(P,mapping.residue_name),1):
        x,y,z=point
        lines.append(f'ATOM  {i:5d}  CA  {str(name):>3s} A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{20.:6.2f}           C')
    path.write_text('\n'.join(lines)+'\nTER\nEND\n')

def structural_fit(A,B,map_A,map_B,folder,exe):
    """US-align first structure B onto second structure A, rigid sequential mode."""
    folder=Path(folder); folder.mkdir(parents=True,exist_ok=True)
    a=folder/'A.pdb'; b=folder/'B.pdb'; mat=folder/'rotation.txt'
    write_ca_pdb(a,A,map_A); write_ca_pdb(b,B,map_B)
    command=[str(exe),str(b),str(a),'-m',str(mat),'-mol','prot']
    process=subprocess.run(command,capture_output=True,text=True,timeout=180)
    (folder/'stdout.txt').write_text(process.stdout)
    (folder/'stderr.txt').write_text(process.stderr)
    if process.returncode or not mat.exists(): raise RuntimeError(f'US-align failed: {process.stderr}\n{process.stdout}')
    rows=[]
    for line in mat.read_text().splitlines():
        tokens=line.split()
        if len(tokens)==5 and tokens[0] in ['0','1','2']:
            rows.append([float(x) for x in tokens[1:]])
    if len(rows)!=3: raise ValueError('Unrecognized US-align rotation matrix')
    transform=np.array(rows); t=transform[:,0]; R=transform[:,1:]
    assert np.allclose(R@R.T,np.eye(3),atol=1e-6) and np.linalg.det(R)>0.99999
    lines=process.stdout.splitlines()
    marker=next(i for i,s in enumerate(lines) if 'denotes residue pairs' in s)
    seq_B=lines[marker+1].strip(); seq_A=lines[marker+3].strip()
    if len(seq_A)!=len(seq_B): raise ValueError('Invalid structural alignment strings')
    ia=ib=-1; pairs=[]
    for cb,ca in zip(seq_B,seq_A):
        if cb!='-': ib+=1
        if ca!='-': ia+=1
        if ca!='-' and cb!='-': pairs.append([ia,ib])
    if ia+1!=len(A) or ib+1!=len(B): raise ValueError('US-align dropped a coordinate; residue mapping would be invalid')
    scores=[float(x) for x in re.findall(r'TM-score=\s*([0-9.]+)',process.stdout)]
    info=dict(method='US-align, rigid sequential sequence-independent structural alignment',command=command,
        rotation_B_to_A=R.tolist(),translation_B_to_A=t.tolist(),alignment_pairs=pairs,tm_scores=scores,
        aligned_rmsd=float(re.search(r'RMSD=\s*([0-9.]+)',process.stdout).group(1)))
    (folder/'alignment.json').write_text(json.dumps(info,indent=2))
    return B@R.T+t,info

def prepare_small_pairs(root,n=30,recompute=False):
    root=Path(root); directory=root/'output'/'small_pairs'; directory.mkdir(parents=True,exist_ok=True)
    exe=next((root/'tools'/'usalign').glob('*.exe'))
    chains=load_prepared(root); keys=list(chains); akey='1o7j.a'
    cache=directory/'prepared_pairs.json'
    deps=sorted((root/'data'/'prepared').glob('*.residues.csv'))+[exe,Path(__file__).resolve()]
    signature=hashlib.sha256(json.dumps({'n':n,'files':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps}},sort_keys=True).encode()).hexdigest()
    metadata=json.loads(cache.read_text()) if cache.exists() else {}
    if recompute or metadata.get('signature')!=signature:
        records=[]
        for bkey in keys:
            if bkey==akey: continue
            pairid=f'{akey}__{bkey}'; folder=directory/'alignments'/pairid
            # A geometry-independent fragment rule: longest, tie broken by first.
            fa=max(chains[akey]['fragments'],key=lambda f:len(f['xyz']))
            fb=max(chains[bkey]['fragments'],key=lambda f:len(f['xyz']))
            _,parent=structural_fit(fa['xyz'],fb['xyz'],fa['mapping'],fb['mapping'],folder/'parent',exe)
            pa,pb=parent['alignment_pairs'][len(parent['alignment_pairs'])//2]
            na=min(n,len(fa['xyz'])); nb=min(n,len(fb['xyz']))
            sa=max(0,min(pa-na//2,len(fa['xyz'])-na)); sb=max(0,min(pb-nb//2,len(fb['xyz'])-nb))
            A=fa['xyz'][sa:sa+na].copy(); B=fb['xyz'][sb:sb+nb].copy()
            ma=fa['mapping'].iloc[sa:sa+na].copy(); mb=fb['mapping'].iloc[sb:sb+nb].copy()
            aligned,window=structural_fit(A,B,ma,mb,folder/'window',exe)
            np.savez_compressed(folder/'aligned.npz',A=A,B=aligned,B_deposited=B)
            for side,mapping in [('A',ma),('B',mb)]:
                mapping.insert(0,'window_index',np.arange(len(mapping)))
                mapping.to_csv(folder/f'{side}.residues.csv',index=False)
            record=dict(pair_id=pairid,A_key=akey,B_key=bkey,A_printed=chains[akey]['printed'],B_printed=chains[bkey]['printed'],
                A_fragment=fa['id'],B_fragment=fb['id'],A_start=sa,B_start=sb,n_A=na,n_B=nb,
                selection='middle aligned parent residue; contiguous window; no geometric-result filtering',
                folder=str(folder.relative_to(root)),parent_tm_scores=parent['tm_scores'],window_alignment=window)
            records.append(record)
            print(f'{pairid}: {na}/{nb} residues; window TM scores {window["tm_scores"]}',flush=True)
        metadata=dict(signature=signature,n_target=n,records=records)
        cache.write_text(json.dumps(metadata,indent=2))
    pairs={}
    for record in metadata['records']:
        folder=root/record['folder']; data=np.load(folder/'aligned.npz')
        pairs[record['pair_id']]=dict(**record,A=data['A'],B=data['B'],
            map_A=pd.read_csv(folder/'A.residues.csv',keep_default_na=False),map_B=pd.read_csv(folder/'B.residues.csv',keep_default_na=False))
    return chains,pairs,metadata

if __name__=='__main__': prepare_small_pairs(Path(__file__).resolve().parents[1])
