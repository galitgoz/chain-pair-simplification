# Controlled protein CPS-2F follow-up

Pair 1o7j.a__1hfj.c, alpha=1, existing w=16 inputs (22 vertices each), angstrom units. Frozen hashes match diagnostic B; no new subsampling/alignment/cropping. Delta1=18.95299827936545, delta2=19.07576340534633, delta3=0.9993252723713152.

Exactly one new attempt. Only max_transitions changed: 10,000,000 to 20,000,000; configuration equality after restoring this single field was asserted. Limits: 20,000,000 transitions, 800,000 single-chain edges, 100,000 states, 55 s internal soft limit, 60 s external solver allowance, 512 MiB sampled worker-tree RSS. Separate budgets unchanged: startup 60, loading 15, warmup 120, checkpoint 5, validation 30, optional measurements 30, reporting 10 seconds.

## Outcome
Warmup completed: True. Both single-curve graphs completed: 2. Product graph entered: True. Graph execution completed with an optimal solution: False.

Status: transition_budget_limit. Raw solution saved: False. External validation: not_run_no_solver_output. kA=None, kB=None, k=None. Reference agreement: unavailable.

Transitions: 20000001; discovered states: 57332; dequeued: 980; completed-state lower bound: 979; queue: 56352. These are not completion percentages.

Solver: 50.764502300007734 s; external validation: None s; approximate peak sampled process-tree RSS: 325.8828125 MiB. Limiting phase: solver / configuration_graph; resource: transition_budget_limit. Resource limits do not establish infeasibility.

## Interpretation and provenance
certificate=False at every CPS call; unchanged algorithms. Existing heap-boundary instrumentation retained with unmeasured diagnostic overhead. No new optimization instrumentation. Solver finish() internally validates fidelity, endpoints/order and coupling before returning; this remains included in solver time. External validation is a separate phase after atomic raw checkpointing, and establishes feasibility rather than independently proving optimality. Independent optimality is an inherited solver claim supported by prior correctness tests.

Reference k=5 is the saved independent continuous optimum whose saved vertex indices, source hashes and feasibility/coupling validation were verified before launch. It was not supplied to the worker as a bound, pruning rule or stop condition. Objective comparison occurs only after graph validation. No missing graph objective is replaced by the reference.

Existing supervisor/worker reused. Only outer preparation/readiness and reporting code was adapted for this transition-budget configuration; worker, supervisor, common and optimization modules are unchanged in this follow-up. No modules were copied. Code hashes are stored in the attempt configuration. All historical output hashes are unchanged. See provenance.json, phase/resource logs, comparison.csv, validation.json inside the attempt, and coverage_with_provenance.csv.

Coverage: 24 inherited verified configurations + 0 new = 24/32. Inherited rows are not new executions or timing repetitions. No automatic retry, other missing case, dataset expansion or repeated measurement was launched. Stop after this report.

Total process time: 70.51338620000752 s. Failed-attempt solver time is time to termination, not time to solution.

Recommendation: review the new limiting resource and memory headroom before authorizing one further bounded diagnostic; do not retry automatically.

## Resumed launch and preserved history

Readiness_20260916T110754.json passed: CPU 8.2?21.2%, available memory at least 2377.9 MiB, no active worker or sustained CPU-intensive competitor. The earlier decision.json and readiness_20260916T110153.json remain the prior no-launch records. The old status is saved as status_before_launch_20260916T110754.md. Exactly one attempt was executed and no worker remains active. The stop was the product-transition count cap, not the time, state or memory limit; it does not establish infeasibility. All objective values remain missing and no external validation was possible.
