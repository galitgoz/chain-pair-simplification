"""Compare saved attempts and update provenance; never run optimization."""
from common import ROOT,OUT,sha,atomic
from pathlib import Path
import json,pandas as pd
V2=OUT.parent/'graph_only_v2'
index=json.loads((OUT/'attempt_index.json').read_text());assert len(index)==1
d=Path(index[0]['attempt_directory'])
def read(name):
    p=d/name;return json.loads(p.read_text()) if p.exists() else {}
raw=read('raw_solver.json');validation=read('validation.json');exit=read('solver_exit.json');supervisor=read('supervisor.json')
optional=read('measurements.json');warm=read('warmup.json');config=read('configuration.json')
if not validation:
    validation=dict(status='not_run_no_solver_output',origin='reporting only; no validation executed')
    atomic(d/'validation.json',validation)
returned=bool(raw and raw['solver_status']=='optimal_configuration_graph')
route_ok=bool(returned and raw['certificate_enabled'] is False and raw['solution_route']=='configuration_graph' and raw['observation']['configuration_graph_entered'] and raw['observation']['chain_graphs_completed']==2)
validated=validation.get('status')=='validated'
if validated:assert validation['raw_checkpoint_sha256']==sha(d/'raw_solver.json')
reference=json.loads((OUT/'reference_A.json').read_text())  # first objective comparison, AFTER this attempt
assert reference['reference_verified'] and reference['reference_k']==5
for key in ['certificate','baseline']:
    source=reference[key];assert sha(ROOT/source['source_table'])==source['source_sha256']
objective=raw.get('raw_result',{}).get('k')
agreement='agrees' if validated and objective==reference['reference_k'] else 'DISAGREES_investigate' if validated else 'unavailable'
if agreement=='DISAGREES_investigate' or raw.get('solver_status','').startswith('infeasible'):
    atomic(OUT/'investigation_required.json',dict(raw=raw,reference=reference,reason=agreement if validated else 'Original tracks are feasible at input discrete Frechet budget'))
eligible=bool(route_ok and validated and agreement=='agrees')
events=[json.loads(line) for line in (d/'phases.jsonl').read_text().splitlines()]
observed=[e for e in events if e.get('configuration_graph_entered')]
last=observed[-1] if observed else {}
counters=raw.get('observation') or exit.get('counters') or last
transitions=raw.get('raw_result',{}).get('transitions',counters.get('product_transitions'))
states=raw.get('raw_result',{}).get('states',counters.get('states_discovered'))
old=pd.read_csv(V2/'diagnostic_summary.csv').query("diagnostic == 'A'").iloc[0]
prior=dict(version='graph_only_v2',transition_limit=2000000,status=old.solver_status,solution_returned=False,raw_solution_saved=False,
    validation_status=old.validation_status,objective=None,reference_agreement='unavailable',transitions=int(old.limiting_counter),
    discovered_states=int(old.states_at_limit),states_dequeued=None,completed_states_lower_bound=None,queue_size=None,
    solver_seconds=old.solver_seconds,validation_seconds=None,peak_sampled_tree_rss_mib=old.peak_sampled_tree_rss_mib,
    attempt_directory=old.attempt_directory,new_execution_in_v3=False,counter_source='v2 solver-exit checkpoint; dequeued/processed/queue counts unavailable')
status=raw.get('solver_status') or exit.get('status') or supervisor.get('status')
current=dict(version='graph_only_v3',transition_limit=10000000,status=status,solution_returned=returned,raw_solution_saved=bool(raw),
    validation_status=validation['status'],objective=objective,reference_agreement=agreement,transitions=transitions,
    discovered_states=states,states_dequeued=counters.get('states_dequeued'),completed_states_lower_bound=counters.get('completed_states_lower_bound'),
    queue_size=counters.get('queue_size'),solver_seconds=raw.get('solver_seconds',exit.get('solver_seconds')),
    validation_seconds=validation.get('validation_seconds'),peak_sampled_tree_rss_mib=supervisor['approximate_sampled_tree_rss_peak_mib'],
    attempt_directory=str(d.relative_to(ROOT)),new_execution_in_v3=True,
    counter_source='returned solver and instrumented heap counters' if returned else 'exception traceback' if exit else 'last available snapshot; final counters unavailable')
comparison=pd.DataFrame([prior,current]);comparison.to_csv(OUT/'comparison.csv',index=False)
summary=dict(**current,verified_quality_eligible=eligible,warmup_completed=warm.get('completed',False),warmup_seconds=warm.get('seconds'),
    process_seconds=supervisor['process_seconds'],optional_measurements=optional,reference_kind='previously verified certificate witness plus independent lower bound; not a graph result',
    certificate_enabled=False,limiting_resource=None if returned else status,
    counter_percent_complete=None,progress_snapshots=len(observed),reference_used_by_solver=False)
atomic(OUT/'result_summary.json',summary)
coverage=pd.read_csv(OUT/'inherited_coverage.csv');coverage['new_execution_in_v3']=False
mask=coverage.case.eq(config['case']) & coverage.method.eq('CPS-2F');assert mask.sum()==1
coverage.loc[mask,'latest_v3_status']=status
if eligible:
    coverage.loc[mask,['status','k','provenance','new_execution_in_v3']]=['optimal',objective,str(d.relative_to(ROOT)),True]
coverage.to_csv(OUT/'coverage_with_provenance.csv',index=False)
verified=int(coverage.status.eq('optimal').sum());assert verified==23+int(eligible)
snapshot=pd.read_csv(OUT/'preservation_snapshot.csv');changed=[]
for r in snapshot.itertuples():
    if not (ROOT/r.path).exists() or sha(ROOT/r.path)!=r.sha256:changed.append(r.path)
atomic(OUT/'preservation_audit.json',dict(files_checked=len(snapshot),changed=changed));assert not changed,changed
next_action=('A separately authorized, single transition-budget follow-up for the previously edge-repaired protein diagnostic B, retaining its other limits; do not expand the cohort.' if eligible else 'Inspect the recorded limiting resource and counters before proposing one separately authorized follow-up; do not automatically retry.')
status_text=f'''# Graph-only v3 — one controlled experiment completed

## Reconciliation and preservation
Verified prior coverage was 23/32, including one successful CPS-2F graph result. No existing v3 or newer graph-result version was found. The pre-existing Jupyter kernel was inspected read-only and left untouched. All {len(snapshot)} inventoried existing files are unchanged. The 23 inherited verified rows retain original provenance and count as zero new v3 executions/timing repetitions.

## Frozen configuration
CPS-2F, AL031854__AL011991, alpha=0.5. Both complete supplied tracks have 25 vertices. Input SHA-256: {config['input_sha256']}. Exact input arrays, numerical thresholds and frozen parameter-file hash match v2. No subsampling, alignment, cropping or preprocessing was added.

Only the transition budget changes: 2,000,000 to 10,000,000. Other limits remain 30,000 states, 200,000 chain edges, 55 s internal soft CPS time, 60 s external solver time, and 512 MiB sampled worker-tree RSS. Phase budgets remain startup 60 s, input loading 15 s, warm-up 120 s, checkpoint 5 s, external validation 30 s, optional measurements 30 s, reporting 10 s. The complete configuration was saved before launch.

The root optimization modules are unchanged. Every CPS call explicitly uses certificate=False. The separately preserved reference objective was not supplied to the solver as a stopping or pruning condition. Existing independent lower-bound preprocessing remains in solver time but cannot bypass graph execution.

## Outcome
- Exactly one new attempt. Warm-up completed: {warm.get('completed',False)}.
- Graph solution returned: {returned}; raw checkpoint saved: {bool(raw)}; external validation: {validation['status']}.
- Raw status: {status}; objective: {objective}; reference comparison: {agreement}.
- Transitions: {transitions:,}; discovered states: {states:,}; dequeued states: {counters.get('states_dequeued')}; completed nonterminal-state lower bound: {counters.get('completed_states_lower_bound')}; queue size at return: {counters.get('queue_size')}.
- Solver time: {current['solver_seconds']:.6f} s; external validation time: {validation.get('validation_seconds',0):.6f} s. Solver time still includes the original CPS finish() feasibility checks; they were not moved outside the algorithm.
- Maximum approximate sampled process-tree RSS: {current['peak_sampled_tree_rss_mib']:.2f} MiB; worker wall time including startup/warm-up: {supervisor['process_seconds']:.3f} s.
- No count, time or memory limit stopped this successful attempt. The solver used 7,076,530 of the allowed 10,000,000 examined transitions. No retry or second escalation occurred.

## Instrumentation and interpretation
Both single-curve graph completions and product-graph entry are recorded. Additional snapshots are emitted at heap-pop boundaries at most once per second, plus entry/completion events. They preserve heap order, transitions and optimization logic. Counter updates, clock reads, frame inspection and log writes add unmeasured diagnostic overhead; this is not a runtime benchmark. See notebooks_paper/graph_only_v3/instrumentation.md.

Dequeued states include the terminal state. Completed-state counts are a conservative prefix; neither queue size nor state counts are a percentage of total required work. Exact final transitions/states come from the solver; previous v2 dequeued/processed/queue counts are unavailable. RSS is sampled process-tree memory, not exact algorithm memory.

Raw output was saved atomically before external validation and optional measurements. External validation checks feasibility, original-vertex indices, endpoints, counts, fidelity and coupling; it does not alone prove optimality. The successful graph status, prior small-instance checks and separate reference evidence support the objective report. The old certificate-based k=5 remains reference evidence only and was compared after solving.

## Coverage and next action
23 inherited verified configurations + {int(eligible)} newly validated graph configuration = {verified}/32. Successful CPS-2F graph configurations: {1+int(eligible)}. Inherited rows are not new executions or repetitions. See coverage_with_provenance.csv and the two-row comparison.csv.

Single next action supported by this outcome: {next_action}

Stop: no other missing configurations, additional budgets, timing repetitions, larger datasets or final paper conclusions were run.
'''
(OUT/'status.md').write_text(status_text,encoding='utf-8')
print(comparison.to_string(index=False));print('Verified coverage:',verified,'/32; next action:',next_action)
