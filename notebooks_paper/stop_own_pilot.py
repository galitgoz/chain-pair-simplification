"""Stop only this task's pilot runner tree; leave the pre-existing kernel alone."""
import psutil,json
from pathlib import Path
victims={}
for p in psutil.process_iter(['pid','cmdline']):
    c=' '.join(p.info['cmdline'] or [])
    if ('notebooks_paper/execute_stage1.py 02_proteins' in c or 'notebooks_paper\\pilot_worker.py' in c) and p.pid not in [11420,22664]:
        for q in p.children(recursive=True)+[p]:
            assert q.pid not in [11420,22664]
            victims[q.pid]=q
for p in reversed(list(victims.values())):
    try:p.kill()
    except psutil.NoSuchProcess:pass
Path('output/paper_v1/logs/interrupted_runner.json').write_text(json.dumps({'stopped_owned_pids':list(victims),'reason':'Correct worker descendant memory/time accounting before proceeding'},indent=2))
print(list(victims))
