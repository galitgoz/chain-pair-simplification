# Controlled protein CPS-2F follow-up

Pair 1o7j.a__1hfj.c, alpha=1, existing w=16 inputs (22 vertices each), angstrom units. Frozen hashes match diagnostic B; no new subsampling/alignment/cropping. Delta1=18.95299827936545, delta2=19.07576340534633, delta3=0.9993252723713152.

Exactly one new attempt. Only product-transition cap changed: 500,000 to 10,000,000. Limits: 800,000 single-chain edges, 30,000 states, 55 s internal soft limit, 60 s external solver allowance, 512 MiB sampled worker-tree RSS. Separate budgets unchanged: startup 60, loading 15, warmup 120, checkpoint 5, validation 30, optional measurements 30, reporting 10 seconds.

## Outcome
Warmup completed: True. Both single-curve graphs completed: 2. Product graph entered: True. Graph execution completed with an optimal solution: False.

Status: state_budget_limit. Raw solution saved: False. External validation: not_run_no_solver_output. kA=None, kB=None, k=None. Reference agreement: unavailable.

Transitions: 3985753; discovered states: 30000; dequeued: 307; completed-state lower bound: 306; queue: 29693. These are not completion percentages.

Solver: 13.339564999987488 s; external validation: None s; approximate peak sampled process-tree RSS: 303.5703125 MiB. Limiting phase: solver / configuration_graph; resource: state_budget_limit. Resource limits do not establish infeasibility.

## Interpretation and provenance
certificate=False at every CPS call; unchanged algorithms. Existing heap-boundary instrumentation retained with unmeasured diagnostic overhead. No new optimization instrumentation. Solver finish() internally validates fidelity, endpoints/order and coupling before returning; this remains included in solver time. External validation is a separate phase after atomic raw checkpointing, and establishes feasibility rather than independently proving optimality. Independent optimality is an inherited solver claim supported by prior correctness tests.

Reference k=5 is the saved independent continuous optimum whose saved vertex indices, source hashes and feasibility/coupling validation were verified before launch. It was not supplied to the worker as a bound, pruning rule or stop condition. Objective comparison occurs only after graph validation. No missing graph objective is replaced by the reference.

Existing supervisor/worker reused. Only output-location plumbing was added to common.py and supervisor.py; no modules were copied. Code hashes are stored in the attempt configuration. All historical output hashes are unchanged. See provenance.json, phase/resource logs, comparison.csv, validation.json inside the attempt, and coverage_with_provenance.csv.

Coverage: 24 inherited verified configurations + 0 new = 24/32. Inherited rows are not new executions or timing repetitions. No automatic retry, other missing case, dataset expansion or repeated measurement was launched. Stop after this report.
