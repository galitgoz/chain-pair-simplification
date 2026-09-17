# Controlled protein CPS-2F follow-up

Pair 1o7j.a__1hfj.c, alpha=1, existing w=16 inputs (22 vertices each), angstrom units. Frozen hashes match diagnostic B; no new subsampling/alignment/cropping. Delta1=18.95299827936545, delta2=19.07576340534633, delta3=0.9993252723713152.

Exactly one new attempt. Only max_states changed: 30,000 to 100,000; configuration equality after restoring this single field was asserted. Limits: 10,000,000 transitions, 800,000 single-chain edges, 100,000 states, 55 s internal soft limit, 60 s external solver allowance, 512 MiB sampled worker-tree RSS. Separate budgets unchanged: startup 60, loading 15, warmup 120, checkpoint 5, validation 30, optional measurements 30, reporting 10 seconds.

## Outcome
Warmup completed: True. Both single-curve graphs completed: 2. Product graph entered: True. Graph execution completed with an optimal solution: False.

Status: solver_timeout. Raw solution saved: False. External validation: not_run_no_solver_output. kA=None, kB=None, k=None. Reference agreement: unavailable.

Transitions: 1441943; discovered states: 17567; dequeued: 148; completed-state lower bound: 147; queue: 17419. These are not completion percentages.

Solver: 56.10024400000111 s; external validation: None s; approximate peak sampled process-tree RSS: 299.13671875 MiB. Limiting phase: solver; resource: solver_timeout. Resource limits do not establish infeasibility.

## Interpretation and provenance
certificate=False at every CPS call; unchanged algorithms. Existing heap-boundary instrumentation retained with unmeasured diagnostic overhead. No new optimization instrumentation. Solver finish() internally validates fidelity, endpoints/order and coupling before returning; this remains included in solver time. External validation is a separate phase after atomic raw checkpointing, and establishes feasibility rather than independently proving optimality. Independent optimality is an inherited solver claim supported by prior correctness tests.

Reference k=5 is the saved independent continuous optimum whose saved vertex indices, source hashes and feasibility/coupling validation were verified before launch. It was not supplied to the worker as a bound, pruning rule or stop condition. Objective comparison occurs only after graph validation. No missing graph objective is replaced by the reference.

Existing supervisor/worker reused. Only controlled_followup.py was adapted to prepare and report this state-budget configuration; worker, supervisor, common and optimization modules are unchanged in this follow-up. No modules were copied. Code hashes are stored in the attempt configuration. All historical output hashes are unchanged. See provenance.json, phase/resource logs, comparison.csv, validation.json inside the attempt, and coverage_with_provenance.csv.

Coverage: 24 inherited verified configurations + 0 new = 24/32. Inherited rows are not new executions or timing repetitions. No automatic retry, other missing case, dataset expansion or repeated measurement was launched. Stop after this report.

Total process time: 154.7321290000109 s. Failed-attempt solver time is time to termination, not time to solution.

## Timeout details and next recommendation

The solver raised its internal wall-time ResourceLimit in configuration-graph search after 56.100244 s (55 s soft budget), and saved solver_exit.json. The supervisor subsequently recorded its external 60 s solver allowance being exceeded while the worker was still exiting; sampled phase elapsed time was 62.711111 s and the owned worker tree was terminated. These are two recorded timeout events, not a state-limit exit. Final counters are from the saved exception traceback, not the earlier progress snapshot. Neither count cap nor the sampled memory cap was reached. The 147 completed-state count is a conservative lower bound; 148 states were dequeued and 17,419 remained queued.

Warm-up took 91.761453 s versus 20.898204 s previously; total process time was 154.732129 s versus 35.759099 s previously. Host scheduling or contention is a possible explanation for the slowdown, not an established diagnosis. No state-budget scaling or timing-performance conclusion follows from these two diagnostic runs.

Recommendation: inspect host load and scheduling before authorizing any further run; if a repeat is later authorized, keep this exact configuration rather than escalating budgets based on this timeout. No repeat or further experiment was executed.
