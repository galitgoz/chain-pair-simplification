# Controlled protein CPS-2F follow-up

## Latest readiness decision and executed repeat

Readiness recheck `readiness_20260916T104907.json` found no active experiment, CPU samples 7.8–13.0%, at least 2414 MiB available memory, and no sustained CPU-intensive competitor. Exactly one repeat was launched after that check. The earlier `decision.json` and `investigation.json` remain unchanged and describe the earlier no-launch decision; old summary files are preserved in `pre_launch_history_20260916T104907/`. The new readiness record also contains the launch-time interpreter/package/cache provenance.

The earlier slowdown was not reproduced in this single repeat: warm-up was 9.859 s (previously 91.761 s), and the solver reached 10,000,001 transition checks in 16.365 s (previously 1,441,943 in 56.100 s). This does not establish a cause or a general timing conclusion. The repeat stopped solely at the transition budget; state, time and sampled-memory caps were not reached. No solution or validation was produced, so kA, kB and k remain missing and coverage remains 24/32. Both single-curve graphs completed. Final counters: 57,332 discovered states, 594 dequeued, conservative completed-state lower bound 593, and queue size 56,738. Total process time was 26.930 s, peak sampled worker-tree RSS 325.11 MiB. See comparison.csv for all three attempts.

Pair 1o7j.a__1hfj.c, alpha=1, existing w=16 inputs (22 vertices each), angstrom units. Frozen hashes match diagnostic B; no new subsampling/alignment/cropping. Delta1=18.95299827936545, delta2=19.07576340534633, delta3=0.9993252723713152.

Exactly one new attempt. Exact repeat: configuration equals the first 100,000-state configuration in every field; no resource or experimental parameter changed. Limits: 10,000,000 transitions, 800,000 single-chain edges, 100,000 states, 55 s internal soft limit, 60 s external solver allowance, 512 MiB sampled worker-tree RSS. Separate budgets unchanged: startup 60, loading 15, warmup 120, checkpoint 5, validation 30, optional measurements 30, reporting 10 seconds.

## Outcome
Warmup completed: True. Both single-curve graphs completed: 2. Product graph entered: True. Graph execution completed with an optimal solution: False.

Status: transition_budget_limit. Raw solution saved: False. External validation: not_run_no_solver_output. kA=None, kB=None, k=None. Reference agreement: unavailable.

Transitions: 10000001; discovered states: 57332; dequeued: 594; completed-state lower bound: 593; queue: 56738. These are not completion percentages.

Solver: 16.36467770001036 s; external validation: None s; approximate peak sampled process-tree RSS: 325.10546875 MiB. Limiting phase: solver / configuration_graph; resource: transition_budget_limit. Resource limits do not establish infeasibility.

## Interpretation and provenance
certificate=False at every CPS call; unchanged algorithms. Existing heap-boundary instrumentation retained with unmeasured diagnostic overhead. No new optimization instrumentation. Solver finish() internally validates fidelity, endpoints/order and coupling before returning; this remains included in solver time. External validation is a separate phase after atomic raw checkpointing, and establishes feasibility rather than independently proving optimality. Independent optimality is an inherited solver claim supported by prior correctness tests.

Reference k=5 is the saved independent continuous optimum whose saved vertex indices, source hashes and feasibility/coupling validation were verified before launch. It was not supplied to the worker as a bound, pruning rule or stop condition. Objective comparison occurs only after graph validation. No missing graph objective is replaced by the reference.

Existing supervisor/worker reused. The outer runner was adapted to accept and report an identical repeat; investigate_repeat.py records pre-run observations without solver-loop instrumentation; worker, supervisor, common and optimization modules are unchanged in this follow-up. No modules were copied. Code hashes are stored in the attempt configuration. All historical output hashes are unchanged. See provenance.json, phase/resource logs, comparison.csv, validation.json inside the attempt, and coverage_with_provenance.csv.

Coverage: 24 inherited verified configurations + 0 new = 24/32. Inherited rows are not new executions or timing repetitions. Exactly one authorized diagnostic repeat was launched. No automatic retry, other missing case or dataset expansion was launched. This is not a new experimental configuration or a runtime study repetition. Stop after this report.

Total process time: 26.929635000007693 s. Failed-attempt solver time is time to termination, not time to solution.

Recommendation: review the new limiting resource and memory headroom before authorizing one further bounded diagnostic; do not retry automatically.

## Slowdown investigation

Worker, supervisor, warm-up, validation, instrumentation, and algorithm hashes match across the original attempts and repeat. The original configurations differ only in max_states. Outer preparation/reporting changes do not execute within the solver. All attempts use the same output-local cache policy and warm-up procedure; each new output starts without a cache. No caches were cleared.

The historical records do not include host CPU load, CPU time, paging or disk activity. The full historical environment and executable hashes were also not recorded, so identical package/runtime state cannot be proven retrospectively. Current interpreter, package versions, executable hash and cache settings are saved in environment.json. No observed core-code or parameter mismatch blocks this diagnostic repeat.

Host load was sampled before the repeat only; no new monitoring was inserted in solver loops. Current light load does not prove the cause of the earlier slowdown. Longer resource-sampling gaps and warm-up time indicate a broader wall-clock slowdown, consistent with scheduling, system load or I/O effects, but do not identify a cause. The raised state cap was never reached in the earlier timed-out run and does not explain the unchanged warm-up procedure taking longer.

This is one diagnostic repeat of the same configuration, not a new configuration or a runtime benchmark. Failed solver times mean time to termination. Completed-state counts are conservative lower bounds, not completion percentages. Reference comparison occurs only after saved output passes validation.

Recommended next step: Inspect host CPU scheduling and disk/paging activity in a separately authorized diagnostic before making algorithmic or runtime-performance claims. No further run or diagnostic was performed.
