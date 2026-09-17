"""Monitor one uniquely named attempt tree with phase-specific budgets."""
from common import ROOT,OUT,atomic,sha,stamp
from pathlib import Path
import sys,json,time,uuid,subprocess,psutil,os

def attempt(config,category='attempts',*,output_directory=None,monitor=None):
    output=Path(output_directory).resolve() if output_directory is not None else OUT
    folder=output/category/(config['diagnostic']+'_'+uuid.uuid4().hex)
    folder.mkdir(parents=True,exist_ok=False)
    config={**config,'attempt_id':folder.name,'code_version':{str(p.relative_to(ROOT)):sha(p) for p in Path(__file__).parent.glob('*.py')}}
    atomic(folder/'configuration.json',config)  # all limits committed before launch
    start=time.monotonic();peak=0;limit=None;progress=dict(phase='startup',phase_started_monotonic=start)
    known={};samples=0;last_monitor=-float('inf')
    with (folder/'stdout.log').open('w') as log,(folder/'resources.jsonl').open('w') as resource:
        process=subprocess.Popen([sys.executable,str(Path(__file__).with_name('worker.py')),str(folder)],cwd=ROOT,
                                 stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW,
                                 env={**os.environ,'CPS_EXPERIMENT_OUTPUT':str(output)})
        while process.poll() is None:
            try:
                parent=psutil.Process(process.pid)
                for child in parent.children(recursive=True)+[parent]:known[child.pid]=child
            except psutil.NoSuchProcess:pass
            rss=0
            for child in known.values():
                try:rss+=child.memory_info().rss
                except psutil.NoSuchProcess:pass
            peak=max(peak,rss);samples+=1
            try:
                if (folder/'progress.json').exists():progress=json.loads((folder/'progress.json').read_text())
            except (OSError,ValueError):pass
            now=time.monotonic();phase=progress['phase'];elapsed=now-progress['phase_started_monotonic']
            if monitor and now-last_monitor>=5:
                monitor(dict(attempt_directory=str(folder.relative_to(ROOT)),phase=phase,phase_elapsed_seconds=elapsed,
                    phase_budget_seconds=config['phase_budgets'][phase],worker_progress=progress,tree_rss_mib=rss/1024**2,**stamp()))
                last_monitor=now
            resource.write(json.dumps(dict(**stamp(),phase=phase,phase_elapsed_seconds=elapsed,tree_rss_bytes=rss,owned_pids=list(known)))+'\n');resource.flush()
            if rss>config['memory_mib']*1024**2:limit='memory_limit'
            elif elapsed>config['phase_budgets'][phase]:
                limit={'startup':'startup_failure','input_loading':'startup_failure','warmup':'startup_failure',
                       'solver':'solver_timeout','validation':'validation_timeout','optional_measurements':'optional_measurement_timeout',
                       'checkpoint':'checkpoint_timeout','reporting':'reporting_timeout'}[phase]
            if limit:
                for child in reversed(list(known.values())):
                    try:child.kill()
                    except psutil.NoSuchProcess:pass
                process.wait();break
            time.sleep(.05)
    terminal=dict(status=limit or 'worker_exited',limiting_phase=progress['phase'] if limit else None,
        phase_elapsed_seconds=time.monotonic()-progress['phase_started_monotonic'],process_seconds=time.monotonic()-start,
        approximate_sampled_tree_rss_peak_mib=peak/1024**2,memory_scope='sum of sampled process-tree RSS, not exact algorithm memory',
        process_exit_code=process.returncode,samples=samples,owned_pids=list(known),**stamp())
    atomic(folder/'supervisor.json',terminal)
    if config.get('validation_replay_raw') and not (folder/'validation.json').exists():
        atomic(folder/'validation.json',dict(status='validation_incomplete',reason=limit or 'recovery worker exited',
            validation_seconds=terminal['phase_elapsed_seconds'] if progress['phase']=='validation' else None,
            raw_checkpoint_sha256=config['validation_replay_sha256'],**stamp()))
    if not (folder/'validation.json').exists() and (folder/'raw_solver.json').exists() and (limit=='validation_timeout' or progress['phase']=='validation'):
        atomic(folder/'validation.json',dict(status='validation_incomplete',reason=limit or 'worker exited during validation',
            validation_seconds=terminal['phase_elapsed_seconds'],raw_checkpoint_sha256=sha(folder/'raw_solver.json'),**stamp()))
    if limit=='optional_measurement_timeout' and not (folder/'measurements.json').exists():
        atomic(folder/'measurements.json',dict(status='optional_measurement_timeout',seconds=terminal['phase_elapsed_seconds']))
    return folder
