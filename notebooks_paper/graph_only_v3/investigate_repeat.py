"""Read-only historical/load investigation, then preparation of an identical repeat.
Uses existing controlled_followup.run/report; no algorithm implementations.
"""
import sys,os,json,time,platform,statistics,importlib.metadata as meta
from pathlib import Path
import psutil
from common import ROOT,atomic,sha,stamp
from controlled_followup import read
OUT=ROOT/'output/paper_v1/protein_B_states_100k_repeat1'
PREV=ROOT/'output/paper_v1/protein_B_states_100k'
OLDER=ROOT/'output/paper_v1/protein_B_transition_10m'

def active_check():
    for p in psutil.process_iter(['pid','cmdline']):
        args=p.info['cmdline'] or []
        if p.pid!=os.getpid() and any(Path(a).name in ['worker.py','run_once.py'] for a in args):
            raise SystemExit('Active experiment worker: '+str(p.info))

def prepare():
    active_check()
    if OUT.exists():raise SystemExit('Destination exists; inspect saved artifacts instead of repeating')
    historical=[];configs=[];attempt_configs=[];snapshot=[]
    for folder in [OLDER,PREV]:
        c=read(folder/'configuration.json');configs.append(c)
        d=ROOT/read(folder/'attempt_index.json')['attempt_directory'];ac=read(d/'configuration.json');attempt_configs.append(ac)
        phases=[json.loads(x) for x in (d/'phases.jsonl').read_text().splitlines()]
        resources=[json.loads(x) for x in (d/'resources.jsonl').read_text().splitlines()]
        gaps=[b['monotonic']-a['monotonic'] for a,b in zip(resources,resources[1:])]
        historical.append(dict(output=str(folder.relative_to(ROOT)),configuration_sha256=sha(folder/'configuration.json'),
            summary=read(folder/'result_summary.json'),supervisor=read(d/'supervisor.json'),solver_exit=read(d/'solver_exit.json'),
            sampling_gap_median_seconds=statistics.median(gaps),sampling_gap_max_seconds=max(gaps),
            sampling_gap_mean_seconds=statistics.mean(gaps),phase_transitions=[e for e in phases if e['monotonic']==e['phase_started_monotonic']],
            chain_completion_events=[e for e in phases if e.get('solver_subphase') in ['graph_construction','graph_initialization']],
            cache_path=str(folder/'cache'),cache_files=[dict(path=str(p.relative_to(folder)),size=p.stat().st_size,mtime=p.stat().st_mtime) for p in (folder/'cache').rglob('*') if p.is_file()],
            historical_interpreter='Project .venv/Scripts/python.exe in saved launch/provenance; Python 3.14 cache suffix. Full per-attempt environment/executable hash was not recorded.'))
        for p in folder.rglob('*'):
            if p.is_file():snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    assert {**configs[1],'graph_limits':{**configs[1]['graph_limits'],'max_states':30000}}==configs[0]
    protected=['common.py','worker.py','supervisor.py']
    matches=[]
    for name in protected:
        key=str(Path('notebooks_paper/graph_only_v3')/name)
        hashes=[c['code_version'][key] for c in attempt_configs]
        assert hashes[0]==hashes[1]==sha(ROOT/key)
        matches.append(dict(file=key,hash=hashes[0],identical=True))
    config=configs[1]
    for name,digest in config['algorithms'].items():assert sha(ROOT/name)==digest
    assert sha(ROOT/config['input_file'])==config['input_sha256']
    # Sample outside experiment and stop sampling before worker launch.
    processes={}
    for p in psutil.process_iter(['pid','name']):
        try:p.cpu_percent();processes[p.pid]=p
        except (psutil.AccessDenied,psutil.NoSuchProcess):pass
    psutil.cpu_percent();samples=[]
    for _ in range(6):
        time.sleep(2);vm=psutil.virtual_memory();top=[]
        for p in processes.values():
            try:top.append(dict(pid=p.pid,name=p.name(),cpu_percent_one_core=p.cpu_percent()))
            except (psutil.AccessDenied,psutil.NoSuchProcess):pass
        samples.append(dict(**stamp(),cpu_percent=psutil.cpu_percent(),memory_available_mib=vm.available/1024**2,
            memory_percent=vm.percent,top_processes=sorted(top,key=lambda p:p['cpu_percent_one_core'],reverse=True)[:10]))
    active_check()
    busy_python=[s for s in samples if any(p['name'].lower().startswith('python') and p['pid']!=os.getpid() and p['cpu_percent_one_core']>20 for p in s['top_processes'])]
    blocked=statistics.mean(s['cpu_percent'] for s in samples)>80 or min(s['memory_available_mib'] for s in samples)<1024 or len(busy_python)>=3
    OUT.mkdir(exist_ok=False)
    atomic(OUT/'investigation.json',dict(historical=historical,core_code_matches=matches,host_load=samples,repeat_informative=not blocked,
        limitations=['Historical CPU load, CPU time, disk I/O, page faults, clock frequency and full environments were not recorded; current load cannot establish prior cause.',
                    'Only outer prepare/run/report code changed between old attempts; worker execution, instrumentation, warmup and timing boundaries are byte-identical.',
                    'Both historical runs used output-local caches. Repeat uses the same fresh-output-local cache policy, without clearing or modifying historical caches.'],
        monitoring_change='Pre-run host sampling only. Existing worker-tree RSS monitoring and solver instrumentation unchanged. Sampling ends before warmup/solving.'))
    env=dict(**stamp(),executable=sys.executable,executable_sha256=sha(sys.executable),python=sys.version,platform=platform.platform(),logical_cpus=psutil.cpu_count(),
        packages={n:meta.version(n) for n in ['numpy','numba','llvmlite','psutil','pandas']},
        cache_configuration={'NUMBA_CACHE_DIR':str(OUT/'cache'),'MPLCONFIGDIR':str(OUT/'cache/matplotlib'),'initial_cache_exists':False},
        thread_environment={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']})
    atomic(OUT/'environment.json',env);atomic(OUT/'historical_results_snapshot.json',snapshot)
    if blocked:
        (OUT/'status.md').write_text('Repeat not launched: sustained heavy host load, insufficient available memory, or another active Python computation in pre-run observations. See investigation.json. Coverage remains 24/32.')
        raise SystemExit('BLOCKED: no experiment launched')
    atomic(OUT/'configuration.json',config)
    atomic(OUT/'provenance.json',dict(**stamp(),previous_output=str(PREV.relative_to(ROOT)),identical_repeat=True,experimental_configuration_changes=[],
        previous_configuration_sha256=sha(PREV/'configuration.json'),input_sha256=config['input_sha256'],reference_supplied_to_solver=False,
        existing_subsampling_w=16,inherited_verified=24,code_and_environment='environment.json and attempt configuration code_version',pre_run_observations='investigation.json'))
    import pandas as pd
    cov=pd.read_csv(PREV/'coverage_with_provenance.csv');assert cov.status.eq('optimal').sum()==24
    cov.to_csv(OUT/'inherited_coverage.csv',index=False);cov[cov.status.eq('optimal')].to_csv(OUT/'inherited_verified_results.csv',index=False)
    atomic(OUT/'reference_evidence.json',read(PREV/'reference_evidence.json'))
    print(json.dumps(dict(ready=True,cpu_mean=statistics.mean(s['cpu_percent'] for s in samples),available_mib_min=min(s['memory_available_mib'] for s in samples),historical_sample_gaps=[{k:h[k] for k in ['sampling_gap_median_seconds','sampling_gap_max_seconds','sampling_gap_mean_seconds']} for h in historical]),indent=2))

def finish():
    active_check()
    investigation=read(OUT/'investigation.json');summary=read(OUT/'result_summary.json')
    diagnosis='''\n## Slowdown investigation\n\nWorker, supervisor, warm-up, validation, instrumentation, and algorithm hashes match across the original attempts and repeat. The original configurations differ only in max_states. Outer preparation/reporting changes do not execute within the solver. All attempts use the same output-local cache policy and warm-up procedure; each new output starts without a cache. No caches were cleared.\n\nThe historical records do not include host CPU load, CPU time, paging or disk activity. The full historical environment and executable hashes were also not recorded, so identical package/runtime state cannot be proven retrospectively. Current interpreter, package versions, executable hash and cache settings are saved in environment.json. No observed core-code or parameter mismatch blocks this diagnostic repeat.\n\nHost load was sampled before the repeat only; no new monitoring was inserted in solver loops. Current light load does not prove the cause of the earlier slowdown. Longer resource-sampling gaps and warm-up time indicate a broader wall-clock slowdown, consistent with scheduling, system load or I/O effects, but do not identify a cause. The raised state cap was never reached in the earlier timed-out run and does not explain the unchanged warm-up procedure taking longer.\n\nThis is one diagnostic repeat of the same configuration, not a new configuration or a runtime benchmark. Failed solver times mean time to termination. Completed-state counts are conservative lower bounds, not completion percentages. Reference comparison occurs only after saved output passes validation.\n'''
    recommendation='Inspect host CPU scheduling and disk/paging activity in a separately authorized diagnostic before making algorithmic or runtime-performance claims.'
    diagnosis+='\nRecommended next step: '+recommendation+' No further run or diagnostic was performed.\n'
    with (OUT/'status.md').open('a',encoding='utf-8') as f:f.write(diagnosis)
    atomic(OUT/'completion_check.json',dict(active_experiment_workers=False,new_attempts=1,
        exact_configuration_equal=read(OUT/'configuration.json')==read(PREV/'configuration.json'),
        all_historical_files_unchanged=all(sha(ROOT/r['path'])==r['sha256'] for r in read(OUT/'historical_results_snapshot.json')),
        repeat_not_new_configuration=True,reference_used_by_solver=False))
    print('Audit and diagnosis saved. No further attempt.')

def blocked_report():
    import pandas as pd
    assert not (OUT/'launch.json').exists() and not list((OUT/'attempts').glob('*/configuration.json'))
    inv=read(OUT/'investigation.json')
    inv['repeat_informative']=False
    inv['decision']='Do not launch: manual review found sustained multi-core PyCharm computation and high/variable total CPU load. Initial simple >80% mean gate was insufficient; reviewed actual process samples before any launch.'
    atomic(OUT/'investigation.json',inv)
    fields=['warmup_seconds','single_curve_graphs_completed','configuration_graph_entered','graph_completed','status','kA','kB','k','transitions','states','states_dequeued','completed_states_lower_bound','queue_size','solver_seconds','validation_seconds','process_seconds','peak_sampled_tree_rss_mib','limiting_resource']
    rows=[]
    for label,folder in [('30,000-state run',OLDER),('First 100,000-state run',PREV)]:
        s=read(folder/'result_summary.json');rows.append(dict(attempt=label,**{k:s.get(k) for k in fields},provenance=str(folder.relative_to(ROOT))))
    rows.append(dict(attempt='100,000-state repeat1',status='not_executed_host_load',limiting_resource='pre-run competing computation; not a solver resource exit'))
    pd.DataFrame(rows).to_csv(OUT/'comparison.csv',index=False)
    cov=pd.read_csv(PREV/'coverage_with_provenance.csv');cov['new_execution_this_investigation']=False
    cov.to_csv(OUT/'coverage_with_provenance.csv',index=False)
    unchanged=all(sha(ROOT/r['path'])==r['sha256'] for r in read(OUT/'historical_results_snapshot.json'));assert unchanged
    atomic(OUT/'decision.json',dict(**stamp(),launch_authorized_by_conditions=False,launched=False,new_attempts=0,
        blocker='Sustained competing PyCharm computation and high/variable CPU load in final pre-run sample',
        mean_cpu_percent=statistics.mean(s['cpu_percent'] for s in inv['host_load']),
        cpu_range=[min(s['cpu_percent'] for s in inv['host_load']),max(s['cpu_percent'] for s in inv['host_load'])],
        inherited_verified=24,new_verified=0,total_verified=24,previous_artifacts_unchanged=unchanged,
        next_action='Wait until competing work subsides, then request one repeat with the same frozen configuration. No further sampling or repeat now.'))
    text='''# Protein CPS-2F slowdown investigation — repeat not executed

## Decision
No experiment worker was active at the initial process check; only the pre-existing Jupyter kernel was present. Destination had no previous attempt. No experiment was launched. The prepared repeat configuration is exactly equal to protein_B_states_100k/configuration.json; it was saved for review only.

Initial Windows samples were mostly light (4–8% CPU, one 37% sample), but the final pre-run six-sample interval changed the assessment: aggregate CPU 47.7–100%, mean 73.98%. PyCharm PID 11492 consumed 131.3–704.8% of one core across all samples (1.3–7.0 cores). Antivirus, Zoom, browsers and system activity also appeared. Available memory remained 2418–2649 MiB. This is competing computation; the repeat was withheld. A simple automated mean-CPU >80% gate initially passed, but manual inspection of the full samples overrode it before any launch. No processes were terminated, notebooks closed, caches cleared or settings changed.

## Observed historical differences
The two historical configurations differ only in max_states (30,000 vs 100,000). Input SHA-256, thresholds, tolerances (unchanged algorithm code), certificate=False and all other limits match. Worker, supervisor, common orchestration/instrumentation, warm-up and algorithms have identical recorded/current hashes. Only controlled_followup.py, the outer preparation/reporting script, changed to support the state-budget follow-up; it does not execute the optimization in the worker.

Both runs use project .venv/Scripts/python.exe according to saved launch/provenance and Python 3.14 cache suffixes. Historical per-attempt executable hashes and complete package environments were not recorded, so runtime/package identity cannot be proven retrospectively. Current interpreter, executable hash and versions are saved in environment.json. No observed unexplained core-code or parameter mismatch was found.

Both historical runs use fresh output-local NUMBA_CACHE_DIR and MPLCONFIGDIR, with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1, the same synthetic warm-up and unlimited internal warm-up time under its separate 120 s external cap. Cache inventories and modification times are saved in investigation.json. No historical cache was touched. The prepared repeat retains that cache policy but has never warmed or run.

Warm-up rose from 20.898 s to 91.761 s. Resource-monitor gaps grew from mean 0.0738 s / median 0.0698 s / maximum 0.307 s to mean 0.3247 s / median 0.1108 s / maximum 5.709 s, despite the same nominal 0.05 s sampling sleep. Timing boundaries and monitoring code are unchanged. The slow run saved an internal wall-time exit after 56.100 s, then the supervisor reached the external solver allowance at recorded phase elapsed 62.711 s and terminated the exiting worker. Both timeout records are retained.

## Diagnosis and uncertainties
The slowdown extends beyond graph work to warm-up and monitoring responsiveness. Scheduling, competing load, paging or I/O are possible explanations, not established causes. The current competing PyCharm process does not prove it caused the earlier slowdown. Historical host CPU load, worker CPU time, disk/paging activity and clock frequency were not captured. Increasing max_states is not an established cause; the first 100,000-state run stopped with only 17,567 discovered states, and warm-up does not use that experimental state cap.

The only extra monitoring in this investigation was read-only host sampling before any potential execution; no instrumentation was added to solver loops or to the existing worker resource monitor. Current samples include some investigation overhead. No general runtime claim follows from these diagnostic records.

## Results and coverage
See comparison.csv for both prior attempts and the repeat marked not_executed_host_load. Neither prior attempt returned a solution: kA, kB, k and validation times remain missing. Queue and processed-state lower bounds are not completion percentages. Failed-run times are time to termination, not solution time. The repeat has no solver outcome or resource-limit claim because it never started. Reference k=5 was not supplied to any solver or substituted for a missing output.

Coverage remains 24/32, entirely inherited; zero new attempts and zero new verified results. Historical artifacts, including caches, were hash-checked unchanged. Configuration, input, code, cache and environment provenance are saved alongside host samples and decision.json. Only controlled_followup.py was minimally adapted for an identical repeat and the new investigate_repeat.py was added in the existing code directory. No modules were copied and no notebook was changed.

## One next step
Wait until competing work subsides, then request one repeat with the same frozen configuration. No additional sampling, repeat, resource escalation or recommendation was executed after this decision.
'''
    (OUT/'status.md').write_text(text,encoding='utf-8')
    print('Repeat withheld; zero new attempts; coverage 24/32; historical artifacts unchanged.')

def readiness():
    active_check()
    assert not (OUT/'launch.json').exists() and not list((OUT/'attempts').glob('*/configuration.json')), 'Attempt already exists; inspect instead'
    current=read(OUT/'configuration.json');prior=read(PREV/'configuration.json')
    transition_followup=current['graph_limits']['max_transitions']==20000000
    if transition_followup:
        assert {**current,'graph_limits':{**current['graph_limits'],'max_transitions':10000000}}==prior
    else:assert current==prior
    config=read(OUT/'configuration.json')
    assert sha(ROOT/config['input_file'])==config['input_sha256']
    for name,digest in config['algorithms'].items():assert sha(ROOT/name)==digest
    old_attempt=ROOT/read(PREV/'attempt_index.json')['attempt_directory']
    old_code=read(old_attempt/'configuration.json')['code_version']
    for name in ['common.py','worker.py','supervisor.py']:
        key=str(Path('notebooks_paper/graph_only_v3')/name);assert sha(ROOT/key)==old_code[key]
    tracked=[]
    for p in psutil.process_iter(['pid','name']):
        try:
            if p.pid==0:continue
            p.cpu_percent();tracked.append((p,p.info['name']))
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    psutil.cpu_percent();samples=[];begin=time.monotonic()
    for _ in range(10):
        total=psutil.cpu_percent(interval=2);vm=psutil.virtual_memory();top=[]
        for p,name in tracked:
            try:top.append(dict(pid=p.pid,name=name,cpu_percent_one_core=p.cpu_percent()))
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
        samples.append(dict(**stamp(),cpu_percent=total,available_mib=vm.available/1024**2,
            top=sorted(top,key=lambda x:x['cpu_percent_one_core'],reverse=True)[:12]))
    active_check()
    # Review sustained competing work per process as well as aggregate CPU and memory.
    competitors=[]
    for p,name in tracked:
        if p.pid==os.getpid():continue
        values=[next((v['cpu_percent_one_core'] for v in s['top'] if v['pid']==p.pid),0) for s in samples]
        if sum(v>=50 for v in values)>=5:competitors.append(dict(pid=p.pid,name=name,per_sample_cpu=values))
    heavy=sum(s['cpu_percent']>=80 for s in samples)>=5
    ready=not competitors and not heavy and min(s['available_mib'] for s in samples)>=1024
    tag=time.strftime('%Y%m%dT%H%M%S',time.gmtime())
    decision=dict(**stamp(),ready=ready,elapsed_sampling_seconds=time.monotonic()-begin,samples=samples,
        sustained_competitors=competitors,reason='No active experiment, no sustained CPU-intensive competing process, adequate available memory' if ready else 'Sustained competing computation, heavy CPU load, or insufficient memory; do not launch',
        exact_configuration_equal=not transition_followup,only_change={'max_transitions':[10000000,20000000]} if transition_followup else {},core_code_hashes_match=True,input_hash_verified=True,attempt_exists=False,
        environment=dict(executable=sys.executable,executable_sha256=sha(sys.executable),python=sys.version,
            packages={n:meta.version(n) for n in ['numpy','numba','llvmlite','psutil','pandas']},
            cache_path=str(OUT/'cache'),cache_exists=(OUT/'cache').exists()),
        monitoring='Read-only pre-run sampling only; worker monitoring unchanged')
    path=OUT/('readiness_'+tag+'.json');atomic(path,decision)
    print(json.dumps(dict(readiness_file=str(path),ready=ready,reason=decision['reason'],
        cpu_samples=[s['cpu_percent'] for s in samples],min_available_mib=min(s['available_mib'] for s in samples),competitors=competitors),indent=2))
    if ready and not transition_followup:
        # Preserve exact prior reports before the existing report function updates summaries.
        history=OUT/('pre_launch_history_'+tag);history.mkdir(exist_ok=False)
        import shutil
        for name in ['status.md','comparison.csv','coverage_with_provenance.csv']:
            shutil.copy2(OUT/name,history/name)
        atomic(history/'preserved_record_hashes.json',{name:sha(OUT/name) for name in ['decision.json','investigation.json','environment.json','provenance.json']})
    elif not ready:
        with (OUT/'status.md').open('a',encoding='utf-8') as f:f.write('\n\n## Readiness recheck '+tag+'\nRepeat not launched. '+decision['reason']+'. See '+path.name+'. Prior decision and investigation remain unchanged.\n')

def prepare_transition_followup():
    active_check()
    assert not OUT.exists(), 'Destination exists; inspect instead of repeating'
    prior=read(PREV/'configuration.json');result=read(PREV/'result_summary.json')
    assert result['status']=='transition_budget_limit' and result['transitions']==10000001 and result['states']==57332
    config={**prior,'graph_limits':{**prior['graph_limits'],'max_transitions':20000000}}
    assert {**config,'graph_limits':{**config['graph_limits'],'max_transitions':10000000}}==prior
    OUT.mkdir(exist_ok=False)
    atomic(OUT/'configuration.json',config)
    atomic(OUT/'provenance.json',dict(**stamp(),previous_output=str(PREV.relative_to(ROOT)),transition_followup=True,
        only_configuration_change={'graph_limits.max_transitions':[10000000,20000000]},
        previous_configuration_sha256=sha(PREV/'configuration.json'),reference_supplied_to_solver=False,
        input_sha256=config['input_sha256'],existing_subsampling_w=16,inherited_verified=24))
    atomic(OUT/'reference_evidence.json',read(PREV/'reference_evidence.json'))
    import pandas as pd
    cov=pd.read_csv(PREV/'coverage_with_provenance.csv');assert cov.status.eq('optimal').sum()==24
    cov.to_csv(OUT/'inherited_coverage.csv',index=False)
    cov[cov.status.eq('optimal')].to_csv(OUT/'inherited_verified_results.csv',index=False)
    snapshot=[]
    for parent in [OLDER,ROOT/'output/paper_v1/protein_B_states_100k',PREV]:
        for p in parent.rglob('*'):
            if p.is_file():snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    atomic(OUT/'historical_results_snapshot.json',snapshot)
    readiness()

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='finish':finish()
    elif sys.argv[1]=='blocked_report':blocked_report()
    elif sys.argv[1]=='readiness':readiness()
    elif sys.argv[1]=='transition20':
        PREV=OUT
        OUT=ROOT/'output/paper_v1/protein_B_transitions_20m'
        prepare_transition_followup()
