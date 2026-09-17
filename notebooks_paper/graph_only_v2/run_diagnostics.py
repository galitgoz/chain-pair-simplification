"""Exactly A, B, C, once each. Never escalates or reruns automatically."""
from common import OUT,atomic,sha
from supervisor import attempt
import json
gate=json.loads((OUT/'regression_gate.json').read_text());assert gate['failed']==0
references=json.loads((OUT/'references.json').read_text())
assert all(r['reference_verified'] for r in references if r['diagnostic'] in ['A','B'])
marker=OUT/'diagnostic_execution_started.json'
assert not marker.exists(),'Diagnostic attempts already started; inspect saved directories, do not repeat.'
configs=json.loads((OUT/'diagnostics.json').read_text())
assert [c['diagnostic'] for c in configs]==['A','B','C']
atomic(marker,dict(configurations=3,configuration_sha256=sha(OUT/'diagnostics.json'),automatic_retries=False))
index=[]
for config in configs:
    folder=attempt(config)
    index.append(dict(diagnostic=config['diagnostic'],attempt_directory=str(folder),configuration_sha256=sha(folder/'configuration.json')))
    atomic(OUT/'attempt_index.json',index)
    print(config['diagnostic'],folder.name,'raw saved:',(folder/'raw_solver.json').exists(),flush=True)
assert len(index)==3
