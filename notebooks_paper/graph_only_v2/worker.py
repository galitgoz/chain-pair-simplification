"""Phase-separated worker. Raw checkpoints survive subsequent validation/measurement failures."""
from pathlib import Path
import sys,time,json,traceback
DIRECTORY=Path(sys.argv[1]);sys.path.insert(0,str(Path(__file__).parent))
from common import ROOT,OUT,Phases,atomic,sha,solve,warmup,validate,classify,exception_counters,stamp
phases=Phases(DIRECTORY);phases.begin('input_loading')
config=json.loads((DIRECTORY/'configuration.json').read_text());test=config.get('test_mode')
try:
    import numpy as np
    assert sha(ROOT/config['input_file'])==config['input_sha256']
    for name,digest in config['algorithms'].items():assert sha(ROOT/name)==digest
    data=np.load(ROOT/config['input_file']);A,B=data['A'],data['B']
    phases.begin('warmup',internal_time_limit='disabled; external warm-up budget applies')
    started=time.perf_counter();warmup(A.shape[1],config['method']);warm_seconds=time.perf_counter()-started
    atomic(DIRECTORY/'warmup.json',dict(completed=True,seconds=warm_seconds,**stamp()))
    phases.begin('solver');start=time.perf_counter()
    try:
        raw,route=solve(A,B,config,certificate=False,observe=lambda event:phases.note(**event))
    except Exception as error:
        atomic(DIRECTORY/'solver_exit.json',dict(status=classify(error),solver_seconds=time.perf_counter()-start,
            failure_reason=str(error),counters=exception_counters(error),traceback=traceback.format_exc(),**stamp()))
        raise
    elapsed=time.perf_counter()-start
    phases.begin('checkpoint')
    checkpoint=dict(schema_version=2,configuration=config,configuration_sha256=sha(DIRECTORY/'configuration.json'),
        input_sha256=config['input_sha256'],method=config['method'],code_version=config['code_version'],
        solver_status=raw['status'],solution_route=route['route'],raw_result=raw,observation=route,
        solver_seconds=elapsed,resource_limits=dict(graph=config['graph_limits'],phase_budgets=config['phase_budgets'],memory_mib=config['memory_mib']),
        solver_internal_validation='CPS finish() checks endpoints/order, fidelity, and computes/validates discrete coupling before returning; included in solver time. Independent discrete does not perform this finish().',
        external_validation_performed=False,**stamp())
    assert not (DIRECTORY/'raw_solver.json').exists()
    atomic(DIRECTORY/'raw_solver.json',checkpoint)  # BEFORE any external validation or quality measurement
    phases.note(raw_checkpoint_saved=True,raw_checkpoint_sha256=sha(DIRECTORY/'raw_solver.json'))
    if not raw['status'].startswith('optimal'):
        atomic(DIRECTORY/'validation.json',dict(status='not_applicable',reason='No optimal solution returned; infeasibility requires investigation for this pilot'))
    else:
        phases.begin('validation');start=time.perf_counter()
        saved=json.loads((DIRECTORY/'raw_solver.json').read_text())
        try:
            if test=='validation_timeout':time.sleep(5.)
            candidate=saved['raw_result']
            if test=='validation_failure':candidate={**candidate,'indices_A':[-1]}  # in-memory regression injection; checkpoint unchanged
            checked=validate(A,B,candidate,config)
            checked.update(validation_seconds=time.perf_counter()-start,raw_checkpoint_sha256=sha(DIRECTORY/'raw_solver.json'),**stamp())
            atomic(DIRECTORY/'validation.json',checked)
        except Exception as error:
            atomic(DIRECTORY/'validation.json',dict(status='validation_failed',validation_seconds=time.perf_counter()-start,
                failure_reason=str(error),failure_kind=classify(error),raw_checkpoint_sha256=sha(DIRECTORY/'raw_solver.json'),**stamp()))
        if (DIRECTORY/'validation.json').exists() and json.loads((DIRECTORY/'validation.json').read_text())['status']=='validated':
            phases.begin('optional_measurements');start=time.perf_counter()
            try:
                if test=='optional_failure':raise RuntimeError('Injected optional measurement failure')
                from curve_algorithms import continuous_distance
                ia=saved['raw_result']['indices_A'];ib=saved['raw_result']['indices_B']
                lowerA,upperA=continuous_distance(A,A[ia]);lowerB,upperB=continuous_distance(B,B[ib])
                optional=dict(status='completed',continuous_A_lower=lowerA,continuous_A_upper=upperA,
                    continuous_B_lower=lowerB,continuous_B_upper=upperB,seconds=time.perf_counter()-start)
            except Exception as error:optional=dict(status='optional_measurement_failed',reason=str(error),seconds=time.perf_counter()-start)
            atomic(DIRECTORY/'measurements.json',optional)
    phases.begin('reporting');atomic(DIRECTORY/'worker_finished.json',dict(status='finished',**stamp()))
except Exception as error:
    atomic(DIRECTORY/'worker_error.json',dict(phase=phases.phase,status='startup_failure' if phases.phase in ['startup','input_loading','warmup'] else classify(error),reason=str(error),traceback=traceback.format_exc(),**stamp()))
