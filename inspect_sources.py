from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path
import numpy as np
for p in Path('data/pdb_cache').glob('*.cif'):
    d = MMCIF2Dict(str(p))
    keys = ['label_atom_id','auth_asym_id','label_asym_id','label_seq_id','auth_seq_id','label_alt_id','pdbx_PDB_model_num','Cartn_x','Cartn_y','Cartn_z']
    rows = list(zip(*(d['_atom_site.'+k] for k in keys)))
    ca = [r for r in rows if r[0]=='CA' and r[3]!='.']
    print(p.stem, d.get('_struct.title'))
    print('CA chains', {ch:len(set((r[3],r[6]) for r in ca if r[1]==ch)) for ch in sorted(set(r[1] for r in ca))})
    print('sequence', {e:sum(x==e for x in d.get('_entity_poly_seq.entity_id',[])) for e in set(d.get('_entity_poly_seq.entity_id',[]))})
    for raw in Path('data/raw protein backbones').glob(p.stem+'.*.raw'):
        xyz = np.array([[float(v) for v in row.split(',')] for row in raw.read_text().strip().split(';') if row.strip()])
        chain=raw.name.split('.')[1].upper()
        sub=[r for r in ca if r[1]==chain]
        coords=np.array([[float(v) for v in r[7:]] for r in sub])
        dist=np.linalg.norm(xyz[:,None]-coords[None,:],axis=2)
        best=dist.argmin(axis=1)
        print(raw.name, len(xyz), 'direct match max',dist.min(axis=1).max(), 'first/last',sub[best[0]][:7],sub[best[-1]][:7])
