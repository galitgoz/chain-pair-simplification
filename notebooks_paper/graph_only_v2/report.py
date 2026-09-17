"""Saved-result diagnostic report; never runs optimization."""
from common import ROOT,OUT,atomic,sha
import json,pandas as pd,numpy as np
index=json.loads((OUT/'attempt_index.json').read_text());assert len(index)==3
refs={x['diagnostic']:x for x in json.loads((OUT/'references.json').read_text())}
rows=[]
for entry in index:
    from pathlib import Path
    d=Path(entry['attempt_directory']);config=json.loads((d/'configuration.json').read_text())
    def read(name):
        p=d/name;return json.loads(p.read_text()) if p.exists() else {}
    raw=read('raw_solver.json');exit=read('solver_exit.json');supervisor=read('supervisor.json');warm=read('warmup.json')
    progress=read('progress.json');validation=read('validation.json');optional=read('measurements.json')
    if not validation:
        validation=dict(status='not_run_no_solver_output',validation_seconds=None,feasibility_passed=False,origin='post-run reporting; no validation was executed')
        assert not (d/'validation.json').exists();atomic(d/'validation.json',validation)
    returned=bool(raw and raw['solver_status'].startswith('optimal'))
    passed=validation.get('status')=='validated'
    verified=returned and passed
    raw_result=raw.get('raw_result',{})
    claim=bool(raw_result.get('status','').startswith('optimal'))
    limiting=exit.get('status') or (supervisor['status'] if supervisor['status']!='worker_exited' else None)
    if not limiting and read('worker_error.json'):limiting=read('worker_error.json')['status']
    counters=exit.get('counters',{})
    counter=counters.get('product_transitions',counters.get('single_chain_edges_generated',counters.get('states_discovered')))
    k=raw_result.get('k');reference=refs[config['diagnostic']].get('reference_k')
    agreement='not_comparable_no_solution' if not returned else 'reference_unavailable' if reference is None else 'agrees' if k==reference else 'DISAGREES_investigate'
    if agreement=='DISAGREES_investigate' or raw.get('solver_status','').startswith('infeasible'):
        atomic(OUT/f"investigation_{config['diagnostic']}.json",dict(raw=raw,reference=refs[config['diagnostic']],reason=agreement if returned else 'Original inputs provide a feasible pair at the frozen delta3; infeasibility requires investigation'))
    row=dict(diagnostic=config['diagnostic'],pair=config['pair'],alpha=config['alpha'],method=config['method'],
        attempt_directory=str(d.relative_to(ROOT)),new_execution=True,new_timing_repetition=False,
        warmup_completed=warm.get('completed',False),warmup_seconds=warm.get('seconds'),
        solver_returned_solution=returned,raw_solution_saved=bool(raw),solver_status=raw.get('solver_status',exit.get('status',limiting)),
        solution_route=raw.get('solution_route',progress.get('solver_subphase','unavailable')),
        validation_completed=validation['status'] in ['validated','validation_failed'],validation_passed=passed,validation_status=validation['status'],
        objective=k,kA=validation.get('kA'),kB=validation.get('kB'),solver_optimality_claim=claim,
        validation_proves_optimality=False,verified_quality_eligible=verified,
        solver_seconds=raw.get('solver_seconds',exit.get('solver_seconds')),
        validation_seconds=validation.get('validation_seconds'),optional_status=optional.get('status','not_run'),optional_seconds=optional.get('seconds'),
        limiting_phase='solver / '+progress.get('solver_subphase','unknown') if exit else supervisor.get('limiting_phase'),
        limiting_resource=limiting,limiting_counter=counter,states_at_limit=counters.get('states_discovered'),
        single_chain_graphs_completed=raw.get('observation',progress).get('chain_graphs_completed'),
        reference_optimum=reference,reference_agreement=agreement,
        peak_sampled_tree_rss_mib=supervisor['approximate_sampled_tree_rss_peak_mib'],process_seconds=supervisor['process_seconds'],
        input_sha256=config['input_sha256'],configuration_sha256=sha(d/'configuration.json'))
    rows.append(row)
summary=pd.DataFrame(rows);summary.to_csv(OUT/'diagnostic_summary.csv',index=False)
inherited=pd.read_csv(OUT/'inherited_verified_results.csv');assert len(inherited)==22
old=pd.read_csv(OUT.parent/'graph_only_v1/pilot_results.csv')
coverage=old[['case','pair','alpha','method','status','k']].copy()
coverage['provenance']='graph_only_v1';coverage['new_execution']=False
for r in summary.itertuples():
    mask=coverage.pair.eq(r.pair)&coverage.alpha.eq(r.alpha)&coverage.method.eq(r.method)
    coverage.loc[mask,'latest_diagnostic_status']=r.solver_status
    if r.verified_quality_eligible:
        coverage.loc[mask,['status','k','provenance','new_execution']]=['optimal',r.objective,r.attempt_directory,True]
coverage.to_csv(OUT/'coverage_with_provenance.csv',index=False)
changed=[]
snapshot=pd.read_csv(OUT/'preservation_snapshot.csv')
for r in snapshot.itertuples():
    if not (ROOT/r.path).exists() or sha(ROOT/r.path)!=r.sha256:changed.append(r.path)
atomic(OUT/'preservation_audit.json',dict(files_checked=len(snapshot),changed=changed));assert not changed,changed
new=int(summary.verified_quality_eligible.sum());remaining=32-22-new
status=f'''# Graph-only v2 — repair and diagnostic stage completed

## Reconciled and preserved
The current v1 files confirm 100 passing correctness checks, 22 verified solutions out of 32, and 8 successful CPS graph solutions. The 10 missing results comprise 5 product-transition limits, 2 single-chain-edge limits, and 3 external timeouts. No newer graph-result version was found. The v1 timeouts covered solving plus external validation/measurements; their original limiting phase cannot be recovered from those files.

The 22 verified rows are preserved in inherited_verified_results.csv with source hashes and original timings, and are **not new executions or timing repetitions**. Existing code, notebooks and results passed a {len(snapshot)}-file preservation audit. The pre-existing Jupyter kernel was inspected read-only and not stopped. No optimization implementation was changed.

## Repair
Each attempt has a UUID directory. Limits and code/input/configuration hashes are written before launch. Phases are input loading, warm-up, solving, atomic raw checkpoint, external validation, optional measurements, and reporting. The raw checkpoint uses a temporary file, flush/fsync, then atomic replace, before external validation or distance estimates.

The existing CPS `finish()` performs endpoint/order checks, fidelity tests, and discrete coupling reconstruction/validation before returning. This remains inside solver time; it has not been removed or relabelled as external validation. Independent discrete does not call that finish function. External validation establishes feasibility and count consistency, not optimality by itself; the solver optimality claim is retained separately. The small-instance exhaustive tests check objectives.

All CPS calls explicitly use certificate=False. Compilation warm-up has no internal time cap, is performed on synthetic curves, and has its own 120 s external budget. Startup: 60 s; input loading: 15 s; solver: 60 s external, 55 s internal for CPS (soft, checked by the existing algorithm); checkpoint: 5 s; validation: 30 s; optional measurements: 30 s. Independent discrete has no internal time parameter; its existing single-curve edge default of 2,000,000 is unchanged. RSS cap: 512 MiB over the worker tree, sampled approximately, not exact algorithm memory. Watchdogs can overshoot due to OS scheduling.

Regression outcomes: 21/21 passed, including 16 small exhaustive comparisons, checkpoint survival after validation failure, validation timeout and optional measurement failure, plus an independent baseline whose coupling exceeds delta3 without failing its own constraints. Synthetic tests are not pilot repetitions.

## Exactly three new attempts
All three warm-ups completed. See diagnostic_summary.csv for the three-row table and attempts/ for committed configurations, phase/resource logs, raw/exit checkpoints and separate validation records.

- **A**, AL031854__AL011991, alpha=0.5, CPS-2F: still stopped at the increased product-transition cap, after 2,000,001 examined transitions; 5,452 discovered states. Solver {summary.iloc[0].solver_seconds:.6f} s. No solution returned or raw solution checkpoint; validation not run. Solver exit checkpoint preserved. Other graph-work caps unchanged.
- **B**, 1o7j.a__1hfj.c, alpha=1, CPS-2F: the 800,000 single-chain-edge allowance allowed both chain graphs to complete; the next bottleneck was the unchanged 500,000 product-transition cap (stopped at 500,001, with 17,567 discovered states). Solver {summary.iloc[1].solver_seconds:.6f} s. No solution returned; exit checkpoint preserved; validation not run. No automatic second escalation.
- **C**, 1o7j.a__1d9q.d, alpha=0.5, Independent discrete: returned kA=20, kB=19, k=20. Raw output was saved before external validation. Validation passed. Solver {summary.iloc[2].solver_seconds:.6f} s; validation {summary.iloc[2].validation_seconds:.6f} s; optional measurements {summary.iloc[2].optional_seconds:.6f} s. This does not retrospectively identify the phase responsible for the v1 timeout.

Total diagnostic worker wall time: {summary.process_seconds.sum():.3f} s. Maximum sampled tree RSS: {summary.peak_sampled_tree_rss_mib.max():.2f} MiB. No diagnostic hit the time or memory cap. Limits ended A and B; neither outcome is mathematical infeasibility.

## Reference comparison
A's saved certificate-based k=5 witness was revalidated against the identical input hash and thresholds; the independent continuous k=5 lower bound was checked as saved reference evidence. It is a reference optimum, not a graph-run result. B's saved independently optimal continuous k=5 pair was revalidated and satisfies delta3, supporting the same reference optimum. Neither diagnostic returned a graph solution, so no new graph objective can yet be compared with 5. No identical-configuration verified reference was located for C in the prior frozen-protocol versions. No discrepancy or infeasibility report occurred. Feasibility validation alone is not the argument for either reference's optimality.

## Coverage and recommendation — stop here
Inherited verified results: 22. New diagnostic attempts: 3. Newly validated outputs: {new}. Combined coverage: {22+new}/32; remaining missing configurations: {remaining}. See coverage_with_provenance.csv; inherited rows retain their original provenance and are not timing repetitions.

Use this phase-separated pipeline for any future timeout repair; C now completes and preserves its output, but the old timeout cause remains unidentified. A needs a transition-work investigation or a separately authorized count-only increase; more wall time alone would not bypass its count cap. B's first construction bottleneck was cleared, but a future controlled retry would need to address its transition cap, while monitoring states and memory. Do not generalize these two runs into a runtime or complexity conclusion. For the other missing configurations, identify the resource category before proposing a targeted retry; do not rerun all nine automatically.

No further attempts, timing repetitions, subsampling, threshold changes, larger cohort, greedy preprocessing or paper conclusions were run. This stage is stopped after A, B and C.
'''
(OUT/'status.md').write_text(status,encoding='utf-8')
print(summary[['diagnostic','warmup_completed','solver_returned_solution','raw_solution_saved','validation_status','objective','solver_seconds','validation_seconds','limiting_resource','limiting_counter','reference_agreement']].to_string(index=False))
