"""One isolated graph-only or independent-baseline attempt, including startup reporting."""
from pathlib import Path
import sys,json,time,traceback
from graph_core import ROOT,OUT,solve,measures,warmup,TRACE
from cps_paper_algorithms import ResourceLimit
import numpy as np
task=json.loads(Path(sys.argv[1]).read_text());p=task['param'];method=task['method']
def progress(event):Path(task['progress']).write_text(json.dumps(event))
start=time.perf_counter();warm=None;began=False
try:
    data=np.load(ROOT/p['input_file']);A,B=data['A'],data['B']
    warmup(method);warm=time.perf_counter()-start
    # Write timestamps from the actual worker, not the Windows launcher.
    started=time.time();Path(task['ready']).write_text(json.dumps(dict(started_unix=started,pid=__import__('os').getpid())))
    began=True;start=time.perf_counter()
    result=solve(A,B,p['delta1'],p['delta2'],p['delta3'],method,certificate=False,progress=progress)
    elapsed=time.perf_counter()-start
    row={k:v for k,v in result.items() if k not in ['indices_A','indices_B','coupling','k']}
    row.update(solve_seconds=elapsed,warmup_seconds=warm,measured_execution_count=1)
    if result['status']=='optimal':row.update(measures(A,B,result,p['delta1'],p['delta2'],p['delta3'],method))
    row['states_discovered']=result.get('states');row['transitions_examined']=result.get('transitions')
except Exception as e:
    if isinstance(e,ResourceLimit):status=('time_limit' if 'wall-time' in str(e) else 'other_resource_limit') if began else 'startup_resource_limit'
    elif isinstance(e,MemoryError):status='memory_limit'
    else:status='software_error'
    trace=dict(TRACE) if began else dict(phase='startup',configuration_graph_entered=False,chain_graph_calls=0,chain_graphs_completed=0)
    route=trace.get('phase','startup')
    row=dict(status=status,raw_status=status,route=route,solve_seconds=time.perf_counter()-start if began else None,
             warmup_seconds=warm,measured_execution_count=int(began),failure_reason=str(e),traceback=traceback.format_exc(),**trace)
row['certificate_enabled']=False
Path(task['result_path']).write_text(json.dumps(row,indent=2))
