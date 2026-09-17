"""One authorized configuration, one new attempt, no automatic retries."""
from common import ROOT,OUT,atomic,sha
from supervisor import attempt
from pathlib import Path
import json,ast
launch=OUT/'launch.json'
if launch.exists():raise SystemExit('This experiment was already launched. Inspect existing output instead of repeating it.')
config=json.loads((OUT/'configuration.json').read_text())
assert not list((OUT/'attempts').glob('*/configuration.json'))
assert config['graph_limits']==dict(seconds=55.,max_states=30000,max_transitions=10000000,max_chain_edges=200000)
assert config['phase_budgets']['solver']==60 and config['memory_mib']==512 and config['certificate'] is False
assert json.loads((OUT.parent/'graph_only_v2/regression_gate.json').read_text())['passed']==21
# Compile the copied orchestration and statically verify every CPS/adapter call opts out.
calls=[]
for p in Path(__file__).parent.glob('*.py'):
    source=p.read_text(encoding='utf-8');compile(source,str(p),'exec')
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node,ast.Call):continue
        name=node.func.id if isinstance(node.func,ast.Name) else node.func.attr if isinstance(node.func,ast.Attribute) else ''
        if name not in ['solve','solve_cps']:continue
        opt=[k.value for k in node.keywords if k.arg=='certificate']
        assert len(opt)==1 and isinstance(opt[0],ast.Constant) and opt[0].value is False
        calls.append(dict(file=str(p.relative_to(ROOT)),line=node.lineno,certificate=False))
atomic(OUT/'callsite_audit.json',calls)
atomic(launch,dict(configuration_sha256=sha(OUT/'configuration.json'),attempts_authorized=1,automatic_retries=False))
folder=attempt(config)
atomic(OUT/'attempt_index.json',[dict(attempt_directory=str(folder),configuration_sha256=sha(folder/'configuration.json'))])
print('Single attempt finished:',folder,flush=True)
print('Raw solution saved:',(folder/'raw_solver.json').exists(),flush=True)
