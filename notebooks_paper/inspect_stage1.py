"""Read-only audit of prior artifacts; write exclusively to the Stage 1 directory."""
from pathlib import Path
import hashlib,json,datetime,os
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/paper_v1';OUT.mkdir(parents=True,exist_ok=True)
for folder in ['references','inputs','figures/proteins','figures/hurricanes','figures/correctness','logs','cache','workers']:
    (OUT/folder).mkdir(parents=True,exist_ok=True)
os.environ['NUMBA_CACHE_DIR']=str(OUT/'cache')
import pandas as pd
import pymupdf
rows=[]
for pattern,kind in [('*.ipynb','notebook'),('*.py','implementation'),('output/**/*.csv','result_table'),('output/**/manifest.json','manifest')]:
    for p in sorted(ROOT.glob(pattern)):
        if OUT in p.parents:continue
        rows.append(dict(path=str(p.relative_to(ROOT)),kind=kind,bytes=p.stat().st_size,
            modified_utc=datetime.datetime.fromtimestamp(p.stat().st_mtime,datetime.timezone.utc).isoformat(),
            sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(rows).to_csv(OUT/'existing_inventory.csv',index=False)
refdir=Path(r'C:\Users\a\Documents\גלית\מחקר\סקירה\Simplification\CPS research')
refs=[refdir/'M_Sc__Thesis___Galit_Gozoltzani.pdf',*refdir.glob('Chain_Pair_Simplification_under_the_Continuous_*.pdf')]
assert len(refs)==2
records=[]
for i,p in enumerate(refs):
    reader=pymupdf.open(p)
    text='\n\n'.join(f'PDF PAGE {j+1}\n'+page.get_text() for j,page in enumerate(reader))
    (OUT/'references'/f'reference_{i+1}.txt').write_text(text,encoding='utf-8')
    records.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),pages=len(reader)))
(OUT/'references/manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
(OUT/'status.md').write_text('''# Stage 1 status

In progress: inventory and reference inspection; pilot cohort and thresholds will be frozen before optimization.
Only Stage 1 is authorized. Existing notebooks/results are read-only.
Planned: correctness checks, then two protein and two hurricane pairs at alpha 0.5 and 1 with four methods.
Not executed: graph-cost experiment (04), paper-summary conclusions (05), all later stages.
The thesis experiments chapter is not established evidence. Document annotations are reference material only.
''',encoding='utf-8')
print(json.dumps(records,ensure_ascii=True,indent=2))
print('Inventory entries:',len(rows))
