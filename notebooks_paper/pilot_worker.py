from pathlib import Path
import sys,json,time,traceback
from paper_core import ROOT,OUT,solve,measures,warmup
from cps_paper_algorithms import ResourceLimit
import numpy as np
task=json.loads(Path(sys.argv[1]).read_text());p=task['param'];method=task['method']
data=np.load(ROOT/p['input_file']);A,B=data['A'],data['B']
t=time.perf_counter();warmup(method);warm=time.perf_counter()-t
Path(task['ready']).write_text('warmup complete; one measured execution starts')
start=time.perf_counter()
try:
    result=solve(A,B,p['delta1'],p['delta2'],p['delta3'],method)
    elapsed=time.perf_counter()-start
    row=dict(status=result['status'],route=result['route'],solve_seconds=elapsed,warmup_seconds=warm,
        lower_bound=result.get('lower_bound'),states_discovered=result.get('states'),transitions_examined=result.get('transitions'))
    if result['status']=='optimal':row.update(measures(A,B,result,p['delta1'],p['delta2'],p['delta3'],method))
except ResourceLimit as e:
    row=dict(status='time_limit' if 'wall-time' in str(e) else 'other_resource_limit',route='none',
             solve_seconds=time.perf_counter()-start,warmup_seconds=warm,failure_reason=str(e))
except MemoryError:
    row=dict(status='memory_limit',route='none',solve_seconds=time.perf_counter()-start,failure_reason='MemoryError')
except Exception:
    row=dict(status='software_error',route='none',solve_seconds=time.perf_counter()-start,failure_reason=traceback.format_exc())
Path(task['result_path']).write_text(json.dumps(row,indent=2))
