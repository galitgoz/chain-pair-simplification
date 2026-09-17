# Graph-only v3 — one controlled experiment completed

## Reconciliation and preservation
Verified prior coverage was 23/32, including one successful CPS-2F graph result. No existing v3 or newer graph-result version was found. The pre-existing Jupyter kernel was inspected read-only and left untouched. All 1603 inventoried existing files are unchanged. The 23 inherited verified rows retain original provenance and count as zero new v3 executions/timing repetitions.

## Frozen configuration
CPS-2F, AL031854__AL011991, alpha=0.5. Both complete supplied tracks have 25 vertices. Input SHA-256: a2ab35b783cce3d389472b1841e85b92acc3171b29993a2566d04b6759b0b287. Exact input arrays, numerical thresholds and frozen parameter-file hash match v2. No subsampling, alignment, cropping or preprocessing was added.

Only the transition budget changes: 2,000,000 to 10,000,000. Other limits remain 30,000 states, 200,000 chain edges, 55 s internal soft CPS time, 60 s external solver time, and 512 MiB sampled worker-tree RSS. Phase budgets remain startup 60 s, input loading 15 s, warm-up 120 s, checkpoint 5 s, external validation 30 s, optional measurements 30 s, reporting 10 s. The complete configuration was saved before launch.

The root optimization modules are unchanged. Every CPS call explicitly uses certificate=False. The separately preserved reference objective was not supplied to the solver as a stopping or pruning condition. Existing independent lower-bound preprocessing remains in solver time but cannot bypass graph execution.

## Outcome
- Exactly one new attempt. Warm-up completed: True.
- Graph solution returned: True; raw checkpoint saved: True; external validation: validated.
- Raw status: optimal_configuration_graph; objective: 5; reference comparison: agrees.
- Transitions: 7,076,530; discovered states: 5,678; dequeued states: 5678; completed nonterminal-state lower bound: 5677; queue size at return: 0.
- Solver time: 8.986753 s; external validation time: 0.001161 s. Solver time still includes the original CPS finish() feasibility checks; they were not moved outside the algorithm.
- Maximum approximate sampled process-tree RSS: 223.24 MiB; worker wall time including startup/warm-up: 18.988 s.
- No count, time or memory limit stopped this successful attempt. The solver used 7,076,530 of the allowed 10,000,000 examined transitions. No retry or second escalation occurred.

## Instrumentation and interpretation
Both single-curve graph completions and product-graph entry are recorded. Additional snapshots are emitted at heap-pop boundaries at most once per second, plus entry/completion events. They preserve heap order, transitions and optimization logic. Counter updates, clock reads, frame inspection and log writes add unmeasured diagnostic overhead; this is not a runtime benchmark. See notebooks_paper/graph_only_v3/instrumentation.md.

Dequeued states include the terminal state. Completed-state counts are a conservative prefix; neither queue size nor state counts are a percentage of total required work. Exact final transitions/states come from the solver; previous v2 dequeued/processed/queue counts are unavailable. RSS is sampled process-tree memory, not exact algorithm memory.

Raw output was saved atomically before external validation and optional measurements. External validation checks feasibility, original-vertex indices, endpoints, counts, fidelity and coupling; it does not alone prove optimality. The successful graph status, prior small-instance checks and separate reference evidence support the objective report. The old certificate-based k=5 remains reference evidence only and was compared after solving.

## Coverage and next action
23 inherited verified configurations + 1 newly validated graph configuration = 24/32. Successful CPS-2F graph configurations: 2. Inherited rows are not new executions or repetitions. See coverage_with_provenance.csv and the two-row comparison.csv.

Single next action supported by this outcome: A separately authorized, single transition-budget follow-up for the previously edge-repaired protein diagnostic B, retaining its other limits; do not expand the cohort.

Stop: no other missing configurations, additional budgets, timing repetitions, larger datasets or final paper conclusions were run.
