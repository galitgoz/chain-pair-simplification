"""Render saved graph-only evidence; no optimization or later stages."""
from graph_core import ROOT,OUT
from pathlib import Path
import os,sys,json,time
os.environ['JUPYTER_PATH']=str(OUT/'workers/jupyter')
k=OUT/'workers/jupyter/kernels/graph-stage1';k.mkdir(parents=True,exist_ok=True)
(k/'kernel.json').write_text(json.dumps({'argv':[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],'display_name':'Graph-only Stage 1','language':'python'}))
import nbformat
from nbclient import NotebookClient
for p in sorted(Path(__file__).parent.glob('0[0-3]_*.ipynb')):
    nb=nbformat.read(p,as_version=4);start=time.perf_counter()
    try:NotebookClient(nb,timeout=180,kernel_name='graph-stage1',resources={'metadata':{'path':str(ROOT)}}).execute()
    finally:nbformat.write(nb,p)
    print(p.name,round(time.perf_counter()-start,2),'seconds',flush=True)
