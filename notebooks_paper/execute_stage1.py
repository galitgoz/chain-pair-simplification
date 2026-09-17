"""Execute only explicitly named new notebooks in a separate kernel."""
import os,sys,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output/paper_v1'
os.environ['NUMBA_CACHE_DIR']=str(OUT/'cache');os.environ['MPLCONFIGDIR']=str(OUT/'cache/matplotlib')
os.environ['JUPYTER_PATH']=str(OUT/'workers/jupyter')
kernel=OUT/'workers/jupyter/kernels/paper-stage1';kernel.mkdir(parents=True,exist_ok=True)
(kernel/'kernel.json').write_text(json.dumps({'argv':[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],'display_name':'Paper Stage 1','language':'python'}))
import nbformat
from nbclient import NotebookClient
for name in sys.argv[1:]:
    assert name in ['00_data_and_protocol.ipynb','01_correctness.ipynb','02_proteins.ipynb','03_hurricanes.ipynb']
    path=ROOT/'notebooks_paper'/name;nb=nbformat.read(path,as_version=4);start=time.perf_counter()
    try:NotebookClient(nb,timeout=1800,kernel_name='paper-stage1',resources={'metadata':{'path':str(ROOT)}}).execute()
    finally:nbformat.write(nb,path)
    elapsed=time.perf_counter()-start
    with (OUT/'logs/notebook_execution.jsonl').open('a') as f:f.write(json.dumps({'notebook':name,'wall_seconds':elapsed})+'\n')
    print(name,round(elapsed,2),'seconds',flush=True)
