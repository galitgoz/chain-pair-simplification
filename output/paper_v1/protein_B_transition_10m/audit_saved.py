"""Read-only experiment audit and provenance packaging; does not solve."""
import sys,json,pathlib,psutil
out=pathlib.Path(__file__).resolve().parent
root=out.parents[2]
sys.path.insert(0,str(root/'notebooks_paper/graph_only_v3'))
from common import atomic,sha
d=root/json.loads((out/'attempt_index.json').read_text())['attempt_directory']
config=json.loads((d/'configuration.json').read_text())
summary=json.loads((out/'result_summary.json').read_text())
atomic(out/'solver_exit_checkpoint.json',dict(configuration=config,configuration_sha256=sha(d/'configuration.json'),
    solver_exit=json.loads((d/'solver_exit.json').read_text()),summary=summary,source_attempt=str(d),
    certificate_enabled=False,validation_status='not_run_no_solver_output'))
workers=[p.info for p in psutil.process_iter(['pid','cmdline'])
    if any(pathlib.Path(arg).name=='worker.py' for arg in (p.info['cmdline'] or []))]
assert not workers,workers
assert len(list((out/'attempts').glob('*/configuration.json')))==1
atomic(out/'completion_check.json',dict(active_workers=workers,attempts=1,verified_coverage=24,
    validation_executed=False,objective_missing=True,stopped_without_retry=True))
atomic(out/'artifact_hashes.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*')
    if p.is_file() and 'cache' not in p.parts and p.name!='artifact_hashes.json'})
print('One attempt, exit checkpoint saved, no active experiment worker, no retry.')
