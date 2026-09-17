"""Configuration-driven, single-attempt orchestration using the existing worker.

prepare verifies saved evidence; run launches once; report reads saved outputs only.
No optimization implementation or reference-based solver pruning lives here.
"""
import argparse, ast, json, os, time, sys
from pathlib import Path
from common import ROOT, atomic, sha, stamp

def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8'))

def prepare_states(out,previous):
    import numpy as np
    import pandas as pd
    import psutil
    active=[]
    for p in psutil.process_iter(['pid','name','cmdline']):
        info=p.info;args=info.get('cmdline') or []
        if 'python' in (info.get('name') or '').lower():
            active.append(info)
            if p.pid!=os.getpid() and any(Path(a).name in ['worker.py','run_once.py'] for a in args):
                raise SystemExit('Experiment worker active; no changes made')
    assert not out.exists(), 'Destination exists: inspect, do not repeat'
    old=read(previous/'configuration.json');prior=read(previous/'result_summary.json')
    assert prior['status']=='state_budget_limit' and prior['states']==30000
    assert old['graph_limits']==dict(seconds=55.,max_states=30000,max_transitions=10000000,max_chain_edges=800000)
    assert old['certificate'] is False and old['pair']=='1o7j.a__1hfj.c' and old['alpha']==1
    assert sha(ROOT/old['input_file'])==old['input_sha256']
    for name,digest in old['algorithms'].items():assert sha(ROOT/name)==digest
    data=np.load(ROOT/old['input_file']);assert data['A'].shape==data['B'].shape==(22,3)
    config={**old,'graph_limits':{**old['graph_limits'],'max_states':100000}}
    restored={**config,'graph_limits':{**config['graph_limits'],'max_states':30000}}
    assert restored==old
    coverage=pd.read_csv(previous/'coverage_with_provenance.csv');assert coverage.status.eq('optimal').sum()==24
    ref=read(previous/'reference_evidence.json');base=ref['baseline']
    assert ref['reference_verified'] and base['validation']['feasibility_passed'] and base['validation']['coupling_feasible']
    assert sha(ROOT/base['source_table'])==base['source_sha256']
    assert sha(ROOT/'output/paper_v1/graph_only_v1/protein_paths.json')==base['source_paths_sha256']
    snapshot=[]
    for parent in [ROOT/'output/paper_v1/graph_only_v2',ROOT/'output/paper_v1/graph_only_v3',previous]:
        for p in parent.rglob('*'):
            if p.is_file() and not any(x in p.parts for x in ['cache','notebook_kernel']):snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    out.mkdir(parents=True,exist_ok=False)
    atomic(out/'configuration.json',config);atomic(out/'reference_evidence.json',ref)
    coverage.to_csv(out/'inherited_coverage.csv',index=False)
    coverage[coverage.status.eq('optimal')].to_csv(out/'inherited_verified_results.csv',index=False)
    atomic(out/'historical_results_snapshot.json',snapshot)
    atomic(out/'provenance.json',dict(**stamp(),previous_output=str(previous.relative_to(ROOT)),
        previous_configuration_sha256=sha(previous/'configuration.json'),active_python_processes=active,
        only_configuration_change={'graph_limits.max_states':[30000,100000]},exact_other_parameters_equal=True,
        input_sha256=sha(ROOT/config['input_file']),input_shapes={'A':[22,3],'B':[22,3]},existing_subsampling_w=16,
        algorithm_changes=False,reference_supplied_to_solver=False,inherited_verified=24))
    print('Prepared: max_states is the only configuration change (30000 -> 100000).')

def prepare(out):
    import numpy as np
    import pandas as pd
    import psutil
    assert not out.exists(), 'Output already exists: inspect it instead of repeating'
    v2=ROOT/'output/paper_v1/graph_only_v2';v3=v2.parent/'graph_only_v3'
    old=next(c for c in read(v2/'diagnostics.json') if c['diagnostic']=='B')
    matching=[]
    for p in (ROOT/'output/paper_v1').rglob('configuration.json'):
        if 'cache' in p.parts:continue
        c=read(p)
        if c.get('case')==old['case'] and c.get('method')==old['method']:
            matching.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),transition_limit=c.get('graph_limits',{}).get('max_transitions')))
            if c.get('graph_limits',{}).get('max_transitions')==10000000:
                raise SystemExit('Matching follow-up already exists: '+str(p))
    active=[]
    for p in psutil.process_iter(['pid','ppid','name','cmdline']):
        try:
            info=p.info;cmd=' '.join(info.get('cmdline') or [])
            if 'python' in (info.get('name') or '').lower():
                active.append(info)
                if p.pid!=os.getpid() and any(x in cmd for x in ['worker.py','run_once.py','controlled_followup.py run']):
                    raise RuntimeError('An experiment worker/runner is active; stop before launch')
        except (psutil.NoSuchProcess,psutil.AccessDenied):continue
    coverage=pd.read_csv(v3/'coverage_with_provenance.csv')
    assert len(coverage)==32 and coverage.status.eq('optimal').sum()==24
    hurricane=read(v3/'result_summary.json')
    assert hurricane['verified_quality_eligible'] and hurricane['objective']==5 and hurricane['transitions']==7076530
    previous=pd.read_csv(v2/'diagnostic_summary.csv').query("diagnostic == 'B'").iloc[0]
    prevdir=ROOT/previous.attempt_directory
    previous_exit=read(prevdir/'solver_exit.json')
    assert previous_exit['counters']['product_transitions']==500001 and previous.single_chain_graphs_completed==2
    prior_config=read(prevdir/'configuration.json')
    assert all(prior_config[k]==v for k,v in old.items())
    assert sha(ROOT/old['input_file'])==old['input_sha256']
    for name,digest in old['algorithms'].items():assert sha(ROOT/name)==digest
    assert sha(v2/'pilot_parameters.csv')==old['parameter_file_sha256']
    data=np.load(ROOT/old['input_file']);assert data['A'].shape==data['B'].shape==(22,3)
    reference=read(v2/'references/B.json');base=reference['baseline']
    assert reference['reference_verified'] and reference['reference_k']==5
    assert sha(ROOT/base['source_table'])==base['source_sha256']
    paths=ROOT/'output/paper_v1/graph_only_v1/protein_paths.json'
    assert sha(paths)==base['source_paths_sha256']
    saved=read(paths)[old['case']+'|Independent continuous']
    assert all(saved[k]==base['raw_result'][k] for k in ['indices_A','indices_B','coupling'])
    assert base['validation']['status']=='validated' and base['validation']['feasibility_passed'] and base['validation']['coupling_feasible']
    assert base['validation']['coupling_distance']<=old['delta3']+1e-9
    baseline=pd.read_csv(ROOT/base['source_table'])
    row=baseline[baseline.case.eq(old['case']) & baseline.method.eq('Independent continuous')].iloc[0]
    assert row.k==5
    config={**old,'graph_limits':{**old['graph_limits'],'max_transitions':10000000}}
    assert config['graph_limits']==dict(seconds=55.,max_states=30000,max_transitions=10000000,max_chain_edges=800000)
    assert config['certificate'] is False
    out.mkdir(parents=True,exist_ok=False);(out/'cache').mkdir()
    atomic(out/'configuration.json',config)
    coverage.to_csv(out/'inherited_coverage.csv',index=False)
    coverage[coverage.status.eq('optimal')].to_csv(out/'inherited_verified_results.csv',index=False)
    atomic(out/'reference_evidence.json',reference)
    atomic(out/'provenance.json',dict(**stamp(),input_sha256=sha(ROOT/config['input_file']),
        input_shapes={k:list(data[k].shape) for k in ['A','B']},existing_subsampling_w=16,
        thresholds={k:config[k] for k in ['delta1','delta2','delta3']},
        previous_attempt=str(prevdir.relative_to(ROOT)),previous_configuration_sha256=sha(prevdir/'configuration.json'),
        reference_source=str((v2/'references/B.json').relative_to(ROOT)),reference_source_sha256=sha(v2/'references/B.json'),
        reference_verification='Saved feasibility validation and source hashes/selected indices verified; inherited independent optimality claim. Not a graph run.',
        matching_configurations=matching,active_python_processes=active,inherited_verified=24,
        algorithm_changes=False,orchestration_changes=['supervisor accepts explicit output_directory and passes cache environment to existing worker','common accepts CPS_EXPERIMENT_OUTPUT; default unchanged'],
        resource_change={'max_transitions':[500000,10000000]},reference_supplied_to_solver=False))
    snapshot=[]
    for parent in [ROOT/'output/paper_v1/graph_only_v2',ROOT/'output/paper_v1/graph_only_v3']:
        for p in parent.rglob('*'):
            if p.is_file() and not any(x in p.parts for x in ['cache','notebook_kernel']):snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    atomic(out/'historical_results_snapshot.json',snapshot)
    print('Prepared exact frozen B; 24/32 inherited; saved reference feasibility/source hashes verified.',flush=True)

def run(config_path,out):
    from supervisor import attempt
    assert config_path.resolve()==(out/'configuration.json').resolve()
    assert not (out/'launch.json').exists() and not list((out/'attempts').glob('*/configuration.json')), 'Already launched; no retries'
    config=read(config_path)
    assert config['certificate'] is False and config['graph_limits']==dict(seconds=55.,max_states=config['graph_limits']['max_states'],max_transitions=config['graph_limits']['max_transitions'],max_chain_edges=800000)
    assert config['graph_limits']['max_transitions'] in [10000000,20000000]
    assert config['graph_limits']['max_states'] in [30000,100000]
    if config['graph_limits']['max_states']==100000:
        provenance=read(out/'provenance.json')
        parent=ROOT/provenance['previous_output']
        old=read(parent/'configuration.json')
        if provenance.get('transition_followup'):
            assert {**config,'graph_limits':{**config['graph_limits'],'max_transitions':10000000}}==old
            decisions=[read(p) for p in out.glob('readiness_*.json')];assert decisions and decisions[-1]['ready']
        elif provenance.get('identical_repeat'):assert config==old
        else:assert {**config,'graph_limits':{**config['graph_limits'],'max_states':30000}}==old
    assert sha(ROOT/config['input_file'])==config['input_sha256']
    calls=[]
    for name in ['common.py','worker.py','supervisor.py']:
        p=Path(__file__).with_name(name);source=p.read_text();compile(source,str(p),'exec')
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node,ast.Call):continue
            fn=node.func.id if isinstance(node.func,ast.Name) else node.func.attr if isinstance(node.func,ast.Attribute) else ''
            if fn in ['solve','solve_cps']:
                flag=[k.value for k in node.keywords if k.arg=='certificate']
                assert len(flag)==1 and isinstance(flag[0],ast.Constant) and flag[0].value is False
                calls.append(dict(file=name,line=node.lineno,certificate=False))
    atomic(out/'callsite_audit.json',calls)
    atomic(out/'launch.json',dict(**stamp(),configuration_sha256=sha(config_path),authorized_attempts=1,automatic_retries=False))
    folder=attempt(config,output_directory=out)
    atomic(out/'attempt_index.json',dict(attempt_directory=str(folder.relative_to(ROOT))))
    print('One attempt finished:',folder,flush=True)

def report(out):
    import pandas as pd
    import nbformat
    d=ROOT/read(out/'attempt_index.json')['attempt_directory'];config=read(d/'configuration.json')
    def get(name):return read(d/name) if (d/name).exists() else {}
    raw=get('raw_solver.json');ex=get('solver_exit.json');val=get('validation.json');sup=get('supervisor.json');warm=get('warmup.json')
    events=[json.loads(s) for s in (d/'phases.jsonl').read_text().splitlines()]
    observed=[e for e in events if 'solver_subphase' in e];last=observed[-1] if observed else {}
    counters={**last,**raw.get('observation',{}),**ex.get('counters',{})}
    if not val:
        val=dict(status='not_run_no_solver_output',validation_seconds=None,reason='Solver returned no output to validate')
        atomic(d/'validation.json',val)
    result=raw.get('raw_result',{});status=raw.get('solver_status') or ex.get('status') or sup['status']
    valid=val.get('status')=='validated'
    if valid:assert val['raw_checkpoint_sha256']==sha(d/'raw_solver.json') and raw['monotonic']<val['monotonic']
    graph=bool(raw and raw.get('certificate_enabled') is False and raw['solution_route']=='configuration_graph' and status=='optimal_configuration_graph' and counters.get('configuration_graph_entered') and counters.get('chain_graphs_completed')==2)
    ref=read(out/'reference_evidence.json')
    agreement='agrees' if valid and result.get('k')==ref['reference_k'] else 'DISAGREES_investigate' if valid else 'unavailable'
    eligible=bool(graph and valid and agreement=='agrees')
    if agreement=='DISAGREES_investigate' or 'infeasible' in status:
        atomic(out/'investigation_required.json',dict(status=status,reference=ref,raw=raw))
    summary=dict(pair=config['pair'],alpha=config['alpha'],method=config['method'],status=status,
        warmup_completed=warm.get('completed',False),warmup_seconds=warm.get('seconds'),graph_completed=graph,
        configuration_graph_entered=counters.get('configuration_graph_entered',False),single_curve_graphs_completed=counters.get('chain_graphs_completed',0),
        raw_solution_saved=bool(raw),validation_status=val['status'],verified_quality_eligible=eligible,
        kA=result.get('kA',len(result['indices_A']) if 'indices_A' in result else None),
        kB=result.get('kB',len(result['indices_B']) if 'indices_B' in result else None),k=result.get('k'),reference_agreement=agreement,
        transitions=result.get('transitions',counters.get('product_transitions')),states=result.get('states',counters.get('states_discovered')),
        states_dequeued=counters.get('states_dequeued'),completed_states_lower_bound=counters.get('completed_states_lower_bound'),queue_size=counters.get('queue_size'),
        solver_seconds=raw.get('solver_seconds',ex.get('solver_seconds',sup.get('phase_elapsed_seconds') if sup.get('limiting_phase')=='solver' else None)),
        validation_seconds=val.get('validation_seconds'),peak_sampled_tree_rss_mib=sup['approximate_sampled_tree_rss_peak_mib'],
        process_seconds=sup['process_seconds'],limiting_resource=None if graph else status,
        limiting_phase=None if graph else sup.get('limiting_phase') or 'solver / '+counters.get('solver_subphase','unknown'),
        certificate_enabled=False,solution_route=raw.get('solution_route','configuration_graph' if counters.get('configuration_graph_entered') else 'preprocessing_exit'),
        reference_used_by_solver=False,new_executions=1,new_timing_repetitions=0,attempt_directory=str(d.relative_to(ROOT)))
    atomic(out/'result_summary.json',summary)
    previous=pd.read_csv(ROOT/'output/paper_v1/graph_only_v2/diagnostic_summary.csv').query("diagnostic == 'B'").iloc[0]
    comparison=pd.DataFrame([dict(attempt='Previous B',transition_budget=500000,status=previous.solver_status,k=None,
        transitions=previous.limiting_counter,states=previous.states_at_limit,solver_seconds=previous.solver_seconds,
        validation_seconds=None,peak_sampled_tree_rss_mib=previous.peak_sampled_tree_rss_mib),
        dict(attempt='Controlled follow-up',transition_budget=10000000,**{k:summary[k] for k in ['status','k','transitions','states','solver_seconds','validation_seconds','peak_sampled_tree_rss_mib']})])
    comparison.to_csv(out/'comparison.csv',index=False)
    provenance=read(out/'provenance.json')
    if 'previous_output' in provenance:
        previous=read(ROOT/provenance['previous_output']/'result_summary.json')
        fields=['status','kA','kB','k','transitions','states','states_dequeued','completed_states_lower_bound','queue_size','solver_seconds','validation_seconds','process_seconds','peak_sampled_tree_rss_mib']
        comparison=pd.DataFrame([dict(attempt='Previous 30000 states',state_budget=30000,**{k:previous[k] for k in fields}),
                                 dict(attempt='Follow-up 100000 states',state_budget=100000,**{k:summary[k] for k in fields})])
        comparison.to_csv(out/'comparison.csv',index=False)
    if provenance.get('transition_followup'):
        fields=['warmup_seconds','single_curve_graphs_completed','configuration_graph_entered','graph_completed','status','kA','kB','k','transitions','states','states_dequeued','completed_states_lower_bound','queue_size','solver_seconds','validation_seconds','process_seconds','peak_sampled_tree_rss_mib','limiting_resource']
        prior=read(ROOT/provenance['previous_output']/'result_summary.json')
        pd.DataFrame([dict(attempt='Completed 10m repeat',**{k:prior.get(k) for k in fields}),dict(attempt='20m transition follow-up',**{k:summary.get(k) for k in fields})]).to_csv(out/'comparison.csv',index=False)
    if provenance.get('identical_repeat'):
        fields=['warmup_seconds','single_curve_graphs_completed','configuration_graph_entered','graph_completed','status','kA','kB','k','transitions','states','states_dequeued','completed_states_lower_bound','queue_size','solver_seconds','validation_seconds','process_seconds','peak_sampled_tree_rss_mib','limiting_resource']
        rows=[]
        for label,source in [('30,000-state run',ROOT/'output/paper_v1/protein_B_transition_10m'),('First 100,000-state run',ROOT/provenance['previous_output'])]:
            saved=read(source/'result_summary.json');rows.append(dict(attempt=label,**{k:saved.get(k) for k in fields}))
        rows.append(dict(attempt='100,000-state repeat1',**{k:summary.get(k) for k in fields}))
        pd.DataFrame(rows).to_csv(out/'comparison.csv',index=False)
    if not raw:
        atomic(out/'solver_exit_checkpoint.json',dict(configuration=config,configuration_sha256=sha(d/'configuration.json'),
            solver_exit=ex,supervisor=sup,summary=summary,source_attempt=str(d.relative_to(ROOT))))
    coverage=pd.read_csv(out/'inherited_coverage.csv');coverage['new_execution_this_followup']=False
    mask=coverage.case.eq(config['case']) & coverage.method.eq(config['method']);assert mask.sum()==1
    coverage.loc[mask,'latest_followup_status']=status
    if eligible:coverage.loc[mask,['status','k','provenance','new_execution_this_followup']]=['optimal',result['k'],str(d.relative_to(ROOT)),True]
    coverage.to_csv(out/'coverage_with_provenance.csv',index=False)
    verified=int(coverage.status.eq('optimal').sum());assert verified==24+int(eligible)
    unchanged=all(sha(ROOT/r['path'])==r['sha256'] for r in read(out/'historical_results_snapshot.json'));assert unchanged
    for name,digest in config['algorithms'].items():assert sha(ROOT/name)==digest
    for name,digest in config['code_version'].items():assert sha(ROOT/name)==digest
    atomic(out/'final_audit.json',dict(original_results_unchanged=unchanged,algorithms_unchanged=True,
        attempt_code_hashes_unchanged=True,new_attempts=1,inherited_verified=24,new_verified=int(eligible),total_verified=verified,
        successful_CPS2F_graph_configurations=int((coverage.status.eq('optimal') & coverage.method.eq('CPS-2F')).sum())))
    text=f'''# Controlled protein CPS-2F follow-up

Pair {config['pair']}, alpha=1, existing w=16 inputs (22 vertices each), angstrom units. Frozen hashes match diagnostic B; no new subsampling/alignment/cropping. Delta1={config['delta1']!r}, delta2={config['delta2']!r}, delta3={config['delta3']!r}.

Exactly one new attempt. Only product-transition cap changed: 500,000 to 10,000,000. Limits: 800,000 single-chain edges, 30,000 states, 55 s internal soft limit, 60 s external solver allowance, 512 MiB sampled worker-tree RSS. Separate budgets unchanged: startup 60, loading 15, warmup 120, checkpoint 5, validation 30, optional measurements 30, reporting 10 seconds.

## Outcome
Warmup completed: {summary['warmup_completed']}. Both single-curve graphs completed: {summary['single_curve_graphs_completed']}. Product graph entered: {summary['configuration_graph_entered']}. Graph execution completed with an optimal solution: {graph}.

Status: {status}. Raw solution saved: {bool(raw)}. External validation: {val['status']}. kA={summary['kA']}, kB={summary['kB']}, k={summary['k']}. Reference agreement: {agreement}.

Transitions: {summary['transitions']}; discovered states: {summary['states']}; dequeued: {summary['states_dequeued']}; completed-state lower bound: {summary['completed_states_lower_bound']}; queue: {summary['queue_size']}. These are not completion percentages.

Solver: {summary['solver_seconds']} s; external validation: {summary['validation_seconds']} s; approximate peak sampled process-tree RSS: {summary['peak_sampled_tree_rss_mib']} MiB. Limiting phase: {summary['limiting_phase']}; resource: {summary['limiting_resource']}. Resource limits do not establish infeasibility.

## Interpretation and provenance
certificate=False at every CPS call; unchanged algorithms. Existing heap-boundary instrumentation retained with unmeasured diagnostic overhead. No new optimization instrumentation. Solver finish() internally validates fidelity, endpoints/order and coupling before returning; this remains included in solver time. External validation is a separate phase after atomic raw checkpointing, and establishes feasibility rather than independently proving optimality. Independent optimality is an inherited solver claim supported by prior correctness tests.

Reference k=5 is the saved independent continuous optimum whose saved vertex indices, source hashes and feasibility/coupling validation were verified before launch. It was not supplied to the worker as a bound, pruning rule or stop condition. Objective comparison occurs only after graph validation. No missing graph objective is replaced by the reference.

Existing supervisor/worker reused. Only output-location plumbing was added to common.py and supervisor.py; no modules were copied. Code hashes are stored in the attempt configuration. All historical output hashes are unchanged. See provenance.json, phase/resource logs, comparison.csv, validation.json inside the attempt, and coverage_with_provenance.csv.

Coverage: 24 inherited verified configurations + {int(eligible)} new = {verified}/32. Inherited rows are not new executions or timing repetitions. No automatic retry, other missing case, dataset expansion or repeated measurement was launched. Stop after this report.
'''
    if 'previous_output' in provenance:
        text=text.replace('Only product-transition cap changed: 500,000 to 10,000,000. Limits: 800,000 single-chain edges, 30,000 states',
            'Only max_states changed: 30,000 to 100,000; configuration equality after restoring this single field was asserted. Limits: 10,000,000 transitions, 800,000 single-chain edges, 100,000 states')
        text=text.replace('Only output-location plumbing was added to common.py and supervisor.py; no modules were copied.',
            'Only controlled_followup.py was adapted to prepare and report this state-budget configuration; worker, supervisor, common and optimization modules are unchanged in this follow-up. No modules were copied.')
        text+='\nTotal process time: '+str(summary['process_seconds'])+' s. Failed-attempt solver time is time to termination, not time to solution.\n'
        text+='\nRecommendation: review the new limiting resource and memory headroom before authorizing one further bounded diagnostic; do not retry automatically.\n'
    if provenance.get('identical_repeat'):
        text=text.replace('Only max_states changed: 30,000 to 100,000; configuration equality after restoring this single field was asserted.',
            'Exact repeat: configuration equals the first 100,000-state configuration in every field; no resource or experimental parameter changed.')
        text=text.replace('No automatic retry, other missing case, dataset expansion or repeated measurement was launched.',
            'Exactly one authorized diagnostic repeat was launched. No automatic retry, other missing case or dataset expansion was launched. This is not a new experimental configuration or a runtime study repetition.')
        text=text.replace('Only controlled_followup.py was adapted to prepare and report this state-budget configuration;',
            'The outer runner was adapted to accept and report an identical repeat; investigate_repeat.py records pre-run observations without solver-loop instrumentation;')
    if provenance.get('transition_followup'):
        text=text.replace('Only max_states changed: 30,000 to 100,000; configuration equality after restoring this single field was asserted. Limits: 10,000,000 transitions',
            'Only max_transitions changed: 10,000,000 to 20,000,000; configuration equality after restoring this single field was asserted. Limits: 20,000,000 transitions')
        text=text.replace('Only controlled_followup.py was adapted to prepare and report this state-budget configuration;',
            'Only outer preparation/readiness and reporting code was adapted for this transition-budget configuration;')
    (out/'status.md').write_text(text,encoding='utf-8')
    nbpath=Path(__file__).with_name('00_controlled_transition_budget.ipynb');nb=nbformat.read(nbpath,as_version=4)
    marker='## Protein diagnostic B: '+('20-million-transition follow-up' if provenance.get('transition_followup') else '100,000-state repeat1' if provenance.get('identical_repeat') else '100,000-state follow-up' if 'previous_output' in provenance else 'controlled follow-up')
    assert not any(marker in c.source for c in nb.cells), 'Notebook section already exists'
    relative=out.relative_to(ROOT).as_posix()
    nb.cells.append(nbformat.v4.new_markdown_cell(marker+'\n\n'+
        f"Saved result for **{config['pair']}, alpha=1, w=16 (22+22 vertices)**. Limits: {config['graph_limits']['max_states']:,} states, {config['graph_limits']['max_transitions']:,} transitions.\n\n"+
        f"Status: **{status}**. kA={summary['kA']}, kB={summary['kB']}, k={summary['k']}; validation: {val['status']}; reference: {agreement}.\n\n"+
        f"Transitions: {summary['transitions']}; states: {summary['states']}; solver: {summary['solver_seconds']} s; validation: {summary['validation_seconds']} s; sampled RSS: {summary['peak_sampled_tree_rss_mib']} MiB. Coverage: **{verified}/32**, including 24 inherited results.\n\n"+
        f"[Comparison CSV](../../{relative}/comparison.csv) · [Full status](../../{relative}/status.md). This section reads saved results; no optimization is launched."))
    nbformat.write(nb,nbpath)
    print(json.dumps(summary,indent=2));print('Verified coverage',verified,'/32')

COHORT=ROOT/'output/paper_v1/unified_pilot_results'

def cohort_rows():
    import csv
    with (ROOT/'output/paper_v1/unified_pilot_plan/run_manifest.csv').open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))

def cohort_active_check():
    import psutil
    for p in psutil.process_iter(['pid','cmdline']):
        args=p.info['cmdline'] or []
        if p.pid!=os.getpid() and any(Path(a).name in ['worker.py','run_once.py'] for a in args):
            raise RuntimeError('Another experiment worker is active: '+str(p.info))

def cohort_readiness(slot):
    import psutil
    cohort_active_check();tracked=[]
    for p in psutil.process_iter(['pid','name']):
        try:
            if p.pid==0:continue
            p.cpu_percent();tracked.append((p,p.info['name']))
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    samples=[]
    for _ in range(10):
        cpu=psutil.cpu_percent(interval=2);vm=psutil.virtual_memory();top=[]
        for p,name in tracked:
            try:top.append(dict(pid=p.pid,name=name,cpu_one_core=p.cpu_percent()))
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
        samples.append(dict(**stamp(),cpu_percent=cpu,available_mib=vm.available/1024**2,
            top=sorted(top,key=lambda x:x['cpu_one_core'],reverse=True)[:12]))
    cohort_active_check();competing=[]
    for p,name in tracked:
        if p.pid==os.getpid():continue
        values=[next((x['cpu_one_core'] for x in s['top'] if x['pid']==p.pid),0) for s in samples]
        if sum(v>=50 for v in values)>=5:competing.append(dict(pid=p.pid,name=name,cpu_samples=values))
    ready=not competing and sum(s['cpu_percent']>=80 for s in samples)<5 and min(s['available_mib'] for s in samples)>=1024
    decision=dict(**stamp(),ready=ready,samples=samples,competing=competing,
        reason='No sustained competing computation and adequate memory' if ready else 'Sustained competing computation or insufficient memory; checkpoint and stop')
    path=slot/('readiness_'+str(time.time_ns())+'.json');atomic(path,decision)
    print('Readiness',slot.name,ready,'CPU',[s['cpu_percent'] for s in samples],flush=True)
    return ready

def cohort_init():
    import platform,importlib.metadata as meta
    cohort_active_check();rows=cohort_rows();assert len(rows)==32
    COHORT.mkdir(exist_ok=True)
    plan=ROOT/'output/paper_v1/unified_pilot_plan'
    selected={}
    for r in rows:
        candidates=[a for a in json.loads(r['existing_attempts_json']) if a.get('protocol_match') or a.get('effective_match_if_edge_exception_accepted')]
        if candidates:
            chosen=min(candidates,key=lambda a:(a['order_timestamp'],a['config_path']))
            d=ROOT/chosen['attempt_path'];c=read(d/'configuration.json')
            assert sha(d/'configuration.json')==chosen['config_sha256']
            assert sha(ROOT/c['input_file'])==r['input_sha256']
            for key in ['alpha','delta1','delta2','delta3']:assert c[key]==float(r[key])
            for path,digest in c['algorithms'].items():assert sha(ROOT/path)==digest
            assert c['phase_budgets']==dict(startup=60.,input_loading=15.,warmup=120.,solver=60.,checkpoint=5.,validation=30.,optional_measurements=30.,reporting=10.)
            if c['method'].startswith('CPS'):
                assert c['certificate'] is False and c['graph_limits']==dict(seconds=55.,max_states=100000,max_transitions=20000000,max_chain_edges=800000)
            else:assert c['method']=='Independent discrete'
            selected[r['configuration_id']]=chosen['attempt_path']
    assert len(selected)==2,selected
    metadata=dict(manifest_sha256=sha(plan/'run_manifest.csv'),protocol_sha256=sha(plan/'protocol.md'),
        selected_reuse=selected,new_slots=30,policy='Shared external 60 s / 512 MiB budgets; method-specific internal limits explicitly adopted',
        edge_exception='Independent discrete and CPS-3F independent preprocessing retain 2000000 edges per curve; main CPS graphs use 800000.',
        algorithms=json.loads(rows[0]['algorithm_hashes_json']),core_code={n:sha(Path(__file__).with_name(n)) for n in ['common.py','worker.py','supervisor.py']},
        environment=dict(executable=sys.executable,executable_sha256=sha(sys.executable),python=sys.version,platform=platform.platform(),
            packages={n:meta.version(n) for n in ['numpy','numba','llvmlite','psutil','pandas','matplotlib']}),
        cache_policy='Fresh output-local cache for each new slot; warmup unchanged; no shared cohort compilation cache',certificate=False)
    if (COHORT/'cohort.json').exists():
        previous=read(COHORT/'cohort.json')
        for k in ['manifest_sha256','selected_reuse','algorithms','core_code']:assert previous[k]==metadata[k]
    else:atomic(COHORT/'cohort.json',metadata)
    for r in rows:
        slot=COHORT/'slots'/f"{int(r['plan_order']):02d}";slot.mkdir(parents=True,exist_ok=True)
        path=slot/'slot.json'
        if path.exists():continue
        atomic(path,dict(configuration_id=r['configuration_id'],reuse_path=selected.get(r['configuration_id']),plan_order=int(r['plan_order'])))
    print('Verified cohort: 2 reused observations; 30 new slots. No optimization code changed.',flush=True)
    return rows

def cohort_result(r):
    slot=COHORT/'slots'/f"{int(r['plan_order']):02d}";info=read(slot/'slot.json')
    folders=list((slot/'attempts').glob('*/configuration.json'));assert len(folders)<=1
    d=ROOT/info['reuse_path'] if info.get('reuse_path') else folders[0].parent if folders else None
    result={k:r[k] for k in ['configuration_id','plan_order','domain','pair','nA','nB','units','input_file','input_sha256','method','alpha','delta1','delta2','delta3']}
    result.update(origin='reused' if info.get('reuse_path') else 'new' if d else 'unexecuted',attempt_path=str(d.relative_to(ROOT)) if d else '',
        status='unexecuted',validated=False,resource_limited=False,kA=None,kB=None,k=None,retained_A=None,retained_B=None,retained_total=None,
        solver_seconds=None,validation_seconds=None,warmup_seconds=None,total_seconds=None,peak_rss_mib=None,
        states=None,transitions=None,states_dequeued=None,completed_states_lower_bound=None,queue_size=None,chain_graphs_completed=None,
        graph_entered=False,solution_route=None,solver_optimality_claim=False,validation_status='not_run',limiting_phase=None,
        coupling_distance=None,coupling_feasible=None,raw_checkpoint_sha256=None,validation_proves_optimality=False,
        applicable_internal_seconds=55 if r['method'].startswith('CPS') else None,
        effective_chain_cap='2000000 default' if r['method']=='Independent discrete' else '800000 main; 2000000 preprocessing' if r['method']=='CPS-3F' else '800000 main' if r['method']=='CPS-2F' else 'N/A')
    if not d:return result
    def get(name):return read(d/name) if (d/name).exists() else {}
    raw=get('raw_solver.json');ex=get('solver_exit.json');sup=get('supervisor.json');val=get('validation.json');err=get('worker_error.json')
    phases=[json.loads(s) for s in (d/'phases.jsonl').read_text().splitlines()] if (d/'phases.jsonl').exists() else []
    obs=next((e for e in reversed(phases) if 'solver_subphase' in e),{})
    obs={**obs,**raw.get('observation',{}),**ex.get('counters',{})}
    payload=raw.get('raw_result',{});status=raw.get('solver_status') or ex.get('status') or err.get('status') or sup.get('status','attempt_incomplete')
    valid=val.get('status')=='validated'
    if valid:
        assert raw and val['raw_checkpoint_sha256']==sha(d/'raw_solver.json')
        c=read(d/'configuration.json');assert sha(ROOT/c['input_file'])==c['input_sha256']
        if r['method'].startswith('CPS'):
            assert c['certificate'] is False and raw.get('certificate_enabled',False) is False
            assert raw['solver_status']=='optimal_configuration_graph' and raw['solution_route']=='configuration_graph'
            assert obs['configuration_graph_entered'] and obs['chain_graphs_completed']==2
    if not sup:status='attempt_incomplete';valid=False
    outcome='validated' if valid else val.get('status') if val.get('status') in ['validation_failed','validation_incomplete'] else status
    if not valid and err.get('status')=='software_error':outcome='software_error'
    resource=any(s in outcome for s in ['budget_limit','timeout','memory_limit','startup_failure'])
    result.update(status=outcome,raw_status=status,validated=valid,resource_limited=resource,
        solver_seconds=raw.get('solver_seconds',ex.get('solver_seconds',sup.get('phase_elapsed_seconds') if sup.get('limiting_phase')=='solver' else None)),
        validation_seconds=val.get('validation_seconds'),warmup_seconds=get('warmup.json').get('seconds'),total_seconds=sup.get('process_seconds'),
        peak_rss_mib=sup.get('approximate_sampled_tree_rss_peak_mib'),states=payload.get('states',obs.get('states_discovered')),
        transitions=payload.get('transitions',obs.get('product_transitions')),states_dequeued=obs.get('states_dequeued'),
        completed_states_lower_bound=obs.get('completed_states_lower_bound'),queue_size=obs.get('queue_size'),
        chain_graphs_completed=obs.get('chain_graphs_completed'),graph_entered=obs.get('configuration_graph_entered',False),
        solution_route=raw.get('solution_route',obs.get('solver_subphase')),solver_optimality_claim=bool(payload.get('status','').startswith('optimal')),
        validation_status=val.get('status','not_run_no_solution'),limiting_phase=sup.get('limiting_phase') or ('solver / '+obs.get('solver_subphase','unknown') if resource else None),
        coupling_distance=val.get('coupling_distance'),coupling_feasible=val.get('coupling_feasible'),
        raw_checkpoint_sha256=sha(d/'raw_solver.json') if raw else None,supervisor_status=sup.get('status'),
        failure_reason=ex.get('failure_reason',err.get('reason')),solver_time_kind='solution_return' if raw and payload.get('status','').startswith('optimal') else 'time_to_termination')
    if raw:
        result.update(kA=payload.get('kA',len(payload.get('indices_A',[]))),kB=payload.get('kB',len(payload.get('indices_B',[]))),k=payload.get('k'))
    if valid:
        result.update(kA=val['kA'],kB=val['kB'],k=val['k'],retained_A=val['kA']/int(r['nA']),retained_B=val['kB']/int(r['nB']),
            retained_total=(val['kA']+val['kB'])/(int(r['nA'])+int(r['nB'])))
    return result

def cohort_save():
    import pandas as pd
    frame=pd.DataFrame([cohort_result(r) for r in cohort_rows()]);frame.to_csv(COHORT/'results.csv',index=False)
    counts=dict(new_executed=int(frame.origin.eq('new').sum()),reused=int(frame.origin.eq('reused').sum()),validated=int(frame.validated.sum()),
        resource_limited=int(frame.resource_limited.sum()),unexecuted=int(frame.origin.eq('unexecuted').sum()))
    atomic(COHORT/'counts.json',counts);return frame

def cohort_execute():
    from supervisor import attempt
    rows=cohort_init();cohort_save()
    for r in rows:
        slot=COHORT/'slots'/f"{int(r['plan_order']):02d}";info=read(slot/'slot.json')
        if info.get('reuse_path'):continue
        existing=list((slot/'attempts').glob('*/configuration.json'))
        if existing:
            outcome=cohort_result(r)
            if outcome['status'] in ['attempt_incomplete','software_error','validation_failed']:
                atomic(COHORT/'pause.json',dict(reason='Existing attempt needs review; no duplicate launch',configuration_id=r['configuration_id'],**stamp()));break
            continue
        if (slot/'launch.json').exists():
            atomic(COHORT/'pause.json',dict(reason='Launch marker without finalized attempt: inspect; never retry automatically',configuration_id=r['configuration_id'],**stamp()));break
        if not cohort_readiness(slot):
            atomic(COHORT/'pause.json',dict(reason='host_readiness',configuration_id=r['configuration_id'],**stamp()));break
        source=json.loads(r['existing_attempts_json'])[0]['parameter_values']
        config={**source,'method':r['method'],'diagnostic':f"slot{int(r['plan_order']):02d}",'certificate':False,
            'input_sha256':r['input_sha256'],'algorithms':json.loads(r['algorithm_hashes_json']),
            'graph_limits':dict(seconds=55.,max_states=100000,max_transitions=20000000,max_chain_edges=800000),
            'phase_budgets':json.loads(r['phase_budgets_json']),'memory_mib':512,'configuration_id':r['configuration_id']}
        assert sha(ROOT/config['input_file'])==config['input_sha256']
        for n,h in config['algorithms'].items():assert sha(ROOT/n)==h
        atomic(slot/'configuration.json',config)
        atomic(slot/'launch.json',dict(**stamp(),configuration_sha256=sha(slot/'configuration.json'),automatic_retries=False))
        print('START',r['plan_order'],r['configuration_id'],flush=True)
        folder=attempt(config,output_directory=slot)
        atomic(slot/'attempt_index.json',dict(attempt_directory=str(folder.relative_to(ROOT))))
        outcome=cohort_result(r);cohort_save();print('END',r['plan_order'],outcome['status'],outcome['k'],flush=True)
        if outcome['status'] in ['software_error','validation_failed','attempt_incomplete'] or str(outcome['raw_status']).startswith('infeasible'):
            atomic(COHORT/'pause.json',dict(reason='Software/validation/infeasibility review required',configuration_id=r['configuration_id'],status=outcome['status'],**stamp()));break
    cohort_save()
    print('COHORT CHECKPOINT',read(COHORT/'counts.json'),flush=True)

def cohort_report():
    os.environ['MPLCONFIGDIR']=str(COHORT/'report_cache')
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import nbformat
    frame=cohort_save();counts=read(COHORT/'counts.json');figdir=COHORT/'figures';figdir.mkdir(exist_ok=True)
    for column in ['plan_order','nA','nB','alpha','delta1','delta2','delta3']:
        frame[column]=pd.to_numeric(frame[column])
    methods=['Independent continuous','Independent discrete','CPS-2F','CPS-3F'];colors=dict(zip(methods,['#0072B2','#777777','#009E73','#E69F00']))
    plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','text.color':'black','axes.labelcolor':'black','font.size':9,'pdf.fonttype':42})
    def save(fig,name):
        fig.tight_layout();fig.savefig(figdir/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    # Every slot is included in the denominator, including unexecuted rows.
    tab=pd.crosstab(frame.method,frame.status).reindex(methods,fill_value=0);tab.to_csv(COHORT/'completion_by_method.csv')
    fig,ax=plt.subplots(figsize=(10,4));bottom=np.zeros(4)
    for status in tab.columns:
        ax.bar(range(4),tab[status],bottom=bottom,label=status.replace('_',' '));bottom+=tab[status].to_numpy()
    ax.set_xticks(range(4),methods);ax.set_ylim(0,8.6);ax.set_ylabel('Configurations (8 planned per method)');ax.legend(bbox_to_anchor=(1.02,1),loc='upper left',fontsize=8)
    ax.set_title('Validated completion and termination reasons');save(fig,'completion')
    cases=list(dict.fromkeys(frame.configuration_id.str.split('|').str[0]));labels=[c.split(':')[1]+' / '+c.split(':')[2] for c in cases]
    fig,ax=plt.subplots(figsize=(12,5));width=.18
    for j,m in enumerate(methods):
        subset=frame[frame.method.eq(m) & frame.validated]
        x=[cases.index(v.split('|')[0])+(j-1.5)*width for v in subset.configuration_id]
        ax.scatter(x,subset.k,label=m,color=colors[m],s=40)
    ax.set_xticks(range(len(cases)),labels,rotation=90);ax.set_ylabel(r'$k=\max(|A\prime|,|B\prime|)$');ax.set_title('Validated solution sizes; missing outcomes are omitted, never zero');ax.legend(ncol=2);ax.grid(axis='y',alpha=.2);save(fig,'solution_size')
    for field,name,label,cap in [('solver_seconds','solver_time','Solver seconds (return or termination)',60),('peak_rss_mib','memory','Peak sampled process-tree RSS (MiB)',512)]:
        fig,ax=plt.subplots(figsize=(9,4))
        for j,m in enumerate(methods):
            subset=frame[frame.method.eq(m)]
            for t,r in enumerate(subset.itertuples()):
                value=getattr(r,field)
                if pd.isna(value):continue
                marker='o' if r.validated else 'X';face='none' if r.origin=='reused' else colors[m]
                ax.scatter(j+(t-3.5)*.035,value,marker=marker,facecolors=face,edgecolors=colors[m],s=55)
        ax.axhline(cap,ls='--',color='black',lw=1);ax.set_xticks(range(4),methods);ax.set_ylabel(label)
        if field=='solver_seconds':
            ax.set_yscale('log');ax.set_ylabel(label+'; log scale')
        ax.set_title('Circle: validated; X: no validated solution; hollow: reused observation');ax.grid(axis='y',alpha=.2);save(fig,name)
    comparisons=[]
    for case in cases:
        block=frame[frame.configuration_id.str.startswith(case+'|')]
        for joint,baseline in [('CPS-2F','Independent continuous'),('CPS-3F','Independent discrete')]:
            a=block[block.method.eq(joint)].iloc[0];b=block[block.method.eq(baseline)].iloc[0];ok=bool(a.validated and b.validated)
            comparisons.append(dict(case=case,pair=a.pair,alpha=a.alpha,delta1=a.delta1,delta2=a.delta2,delta3=a.delta3,method=joint,baseline=baseline,
                joint_status=a.status,baseline_status=b.status,matched_validated=ok,joint_k=a.k if a.validated else None,baseline_k=b.k if b.validated else None,
                vertex_cost=a.k-b.k if ok else None,k_ratio=a.k/b.k if ok else None))
    costs=pd.DataFrame(comparisons);costs.to_csv(COHORT/'matched_comparisons.csv',index=False)
    examples=[]
    for domain in ['hurricane','protein']:
        choices=frame[frame.domain.eq(domain)&frame.validated&frame.method.str.startswith('CPS')].sort_values('plan_order')
        if choices.empty:continue
        r=choices.iloc[0];d=ROOT/r.attempt_path;raw=read(d/'raw_solver.json')['raw_result'];data=np.load(ROOT/r.input_file);A,B=data['A'],data['B']
        fig=plt.figure(figsize=(7,5));ax=fig.add_subplot(111,projection='3d') if A.shape[1]==3 else fig.add_subplot(111)
        for P,ids,c,label in [(A,raw['indices_A'],'#0072B2','A'),(B,raw['indices_B'],'#D55E00','B')]:
            ax.plot(*P.T,color=c,alpha=.3,lw=1,label=label+' input');ax.plot(*P[ids].T,color=c,lw=2,marker='o',ms=4,label=label+' simplified')
        ax.set_xlabel('x ('+r.units+')');ax.set_ylabel('y ('+r.units+')')
        if A.shape[1]==3:
            ax.set_zlabel('z ('+r.units+')')
            for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_facecolor('white')
        else:ax.set_aspect('equal',adjustable='datalim')
        ax.set_title(f'{r.pair} | alpha={r.alpha} | {r.method}\n'+f'delta1={r.delta1:.6g}, delta2={r.delta2:.6g}, delta3={r.delta3:.6g}')
        ax.legend(fontsize=8);name='geometry_'+domain;save(fig,name);examples.append(dict(configuration_id=r.configuration_id,figure=name,attempt_path=r.attempt_path))
    atomic(COHORT/'geometric_examples.json',examples)
    pairs=costs[costs.matched_validated];findings=[]
    for m in ['CPS-2F','CPS-3F']:
        q=pairs[pairs.method.eq(m)]
        findings.append(f"- {m}: {len(q)}/8 matched validated comparisons; vertex-cost range "+(f"{q.vertex_cost.min():g} to {q.vertex_cost.max():g}." if len(q) else 'unavailable.'))
    table='| Method | Validated / 8 | Resource-limited | Unexecuted |\n|---|---:|---:|---:|\n'
    for m in methods:
        q=frame[frame.method.eq(m)];table+=f'| {m} | {int(q.validated.sum())}/8 | {int(q.resource_limited.sum())} | {int(q.origin.eq("unexecuted").sum())} |\n'
    pause=read(COHORT/'pause.json') if (COHORT/'pause.json').exists() and counts['unexecuted'] else None
    report=f'''# Unified bounded pilot results

## Execution and coverage
New executions: {counts['new_executed']}; reused observations: {counts['reused']}; validated: {counts['validated']}/32 intended configurations; resource-limited: {counts['resource_limited']}; unexecuted: {counts['unexecuted']}. Validated among observed slots: {counts['validated']}/{32-counts['unexecuted']}. All 32 rows, including missing outcomes, are in results.csv. Reused observations retain their original paths and are not newly measured timings. Historical mixed-budget 24/32 coverage is separate and is not this cohort's completion rate.

Pause: {pause if pause else 'none; all planned slots observed'}.

{table}
## Adopted policy
Shared 60 s external solver allowance and 512 MiB sampled worker-tree RSS; method-specific internal limits disclosed. CPS: 55 s soft limit, 100,000 product states, 20,000,000 product transitions, 800,000 edges per main single-curve graph. Independent discrete retains its existing 2,000,000-edge default per curve; CPS-3F retains that default in independent preprocessing. Independent methods ignore CPS-only time/state/transition configuration fields. These are not identical computational-work limits. Startup 60, input loading 15, warm-up 120, checkpoint 5, external validation 30, optional measurements 30, reporting 10 seconds remain unchanged. Each new slot has a fresh local compilation cache, unchanged warm-up and no concurrent worker. CPU/memory readiness observations precede every new attempt. Caps are practical choices, not sufficiency guarantees; no retries or escalation.

Four fixed input pairs, alpha=0.5 and 1.0, frozen full-precision thresholds and arrays. Protein inputs are w=16 (22+22), not full-resolution backbones; hurricane tracks are 25+25 with no additional subsampling/alignment. The protein pairs share 1o7j.a. Different alphas and repeated attempts are not independent input pairs. No threshold tuning, greedy processing or algorithm modifications. All CPS calls explicitly use certificate=False; references never guide search.

## Quality comparisons
{chr(10).join(findings)}

See matched_comparisons.csv for each pair, alpha and exact deltas. CPS-2F is compared with independent continuous; CPS-3F with independent discrete. Independent baselines do not require output coupling <=delta3. Both CPS methods require discrete output coupling, but continuous versus discrete fidelity gives different feasible sets; a size difference is not evidence that methods solve the same constrained problem.

Only validated outputs enter quality figures and retained-vertex fractions. External validation establishes feasibility, not optimality independently; solver optimality claims and routes are separate columns. CPS internal finish() and independent-continuous internal validation remain in solver time. Raw output is atomically saved before external validation. Missing objectives are not zero or infeasible instances; resource stops are completed bounded observations. Optional measurements failing do not invalidate solver checkpoints.

## Figures and interpretation
- figures/completion.pdf/.png: all 8 intended slots per method, including unexecuted and distinct termination reasons.
- figures/solution_size.pdf/.png: validated k on matched inputs/thresholds; omissions are visible through the coverage table.
- figures/solver_time.pdf/.png: circles for validated outcomes, X for other outcomes; hollow points are reused. Failed-run values are time to termination, never successful runtimes; no mean termination time is presented as a speed comparison.
- figures/memory.pdf/.png: approximate sampled process-tree RSS, not exact algorithm memory. Sampling can miss peaks and external time limits may overshoot.
- Geometric illustrations: choose the earliest manifest-order validated CPS result in each domain, at most two. No selection by compression or visual appeal. Input lines are faint, simplified lines bold with selected vertices. They illustrate saved outputs only, not topological preservation. Selections are in geometric_examples.json.

These are single-attempt pilot observations, not a repeated timing benchmark. Historical runtime provenance is incomplete for some reused observations; inspect origin and attempt_path. Current host readiness does not prove past scheduling conditions. Do not infer completion percentages from queue sizes or general speedups from this cohort.

## Implementation and provenance
Existing common.py, worker.py, supervisor.py and optimization algorithms are unchanged. Only controlled_followup.py gained manifest orchestration and saved-result reporting functions. The runner persists a per-slot launch marker, checks for prior attempts, executes sequentially in manifest order, and checkpoints results after each slot. Incomplete launches or systemic software/validation failures stop the batch without retry. Metadata and reuse selection are in cohort.json; each attempt saves code/configuration/input hashes and phase/resource logs. Initial expected reuse was verified: 2 observations and 30 new slots. The original plan and historical results are preserved.

No further experiments or timing repetitions were launched beyond authorized missing slots.
'''
    (COHORT/'report.md').write_text(report,encoding='utf-8')
    notebook=Path(__file__).with_name('00_controlled_transition_budget.ipynb');nb=nbformat.read(notebook,as_version=4)
    marker='## Unified bounded pilot: saved results'
    nb.cells=[c for c in nb.cells if not c.get('metadata',{}).get('unified_pilot_display')]
    cells=[nbformat.v4.new_markdown_cell(marker+'\n\n'+table+'\nReads saved outputs only; never launches optimization. Full constraints, missing results and provenance: [report](../../output/paper_v1/unified_pilot_results/report.md).'),
        nbformat.v4.new_code_cell("from pathlib import Path\nimport pandas as pd\nfrom IPython.display import display, Image\nroot = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p/'output/paper_v1/unified_pilot_results/results.csv').exists())\nsaved = root/'output/paper_v1/unified_pilot_results'\ndisplay(pd.read_csv(saved/'completion_by_method.csv'))\ndisplay(pd.read_csv(saved/'matched_comparisons.csv')[['pair','alpha','method','joint_status','baseline_status','vertex_cost']])\nfor name in ['completion','solution_size','solver_time','memory']:\n    display(Image(filename=str(saved/'figures'/f'{name}.png')))" )]
    for c in cells:c.metadata['unified_pilot_display']=True
    nb.cells.extend(cells);nbformat.write(nb,notebook)
    print('REPORT',counts,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','report','cohort_execute','cohort_report']);p.add_argument('--output',required=True);p.add_argument('--config');p.add_argument('--previous')
    args=p.parse_args();out=Path(args.output).resolve()
    if args.action=='cohort_execute':cohort_execute()
    elif args.action=='cohort_report':cohort_report()
    elif args.action=='prepare':
        if args.previous:prepare_states(out,Path(args.previous).resolve())
        else:prepare(out)
    elif args.action=='run':run(Path(args.config).resolve(),out)
    else:report(out)
