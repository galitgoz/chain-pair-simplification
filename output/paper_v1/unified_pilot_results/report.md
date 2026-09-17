# Unified bounded pilot results

## Execution and coverage
New executions: 30; reused observations: 2; validated: 29/32 intended configurations; resource-limited: 3; unexecuted: 0. Validated among observed slots: 29/32. All 32 rows, including missing outcomes, are in results.csv. Reused observations retain their original paths and are not newly measured timings. Historical mixed-budget 24/32 coverage is separate and is not this cohort's completion rate.

Pause: none; all planned slots observed.

| Method | Validated / 8 | Resource-limited | Unexecuted |
|---|---:|---:|---:|
| Independent continuous | 8/8 | 0 | 0 |
| Independent discrete | 8/8 | 0 | 0 |
| CPS-2F | 5/8 | 3 | 0 |
| CPS-3F | 8/8 | 0 | 0 |

## Adopted policy
Shared 60 s external solver allowance and 512 MiB sampled worker-tree RSS; method-specific internal limits disclosed. CPS: 55 s soft limit, 100,000 product states, 20,000,000 product transitions, 800,000 edges per main single-curve graph. Independent discrete retains its existing 2,000,000-edge default per curve; CPS-3F retains that default in independent preprocessing. Independent methods ignore CPS-only time/state/transition configuration fields. These are not identical computational-work limits. Startup 60, input loading 15, warm-up 120, checkpoint 5, external validation 30, optional measurements 30, reporting 10 seconds remain unchanged. Each new slot has a fresh local compilation cache, unchanged warm-up and no concurrent worker. CPU/memory readiness observations precede every new attempt. Caps are practical choices, not sufficiency guarantees; no retries or escalation.

Four fixed input pairs, alpha=0.5 and 1.0, frozen full-precision thresholds and arrays. Protein inputs are w=16 (22+22), not full-resolution backbones; hurricane tracks are 25+25 with no additional subsampling/alignment. The protein pairs share 1o7j.a. Different alphas and repeated attempts are not independent input pairs. No threshold tuning, greedy processing or algorithm modifications. All CPS calls explicitly use certificate=False; references never guide search.

## Quality comparisons
- CPS-2F: 5/8 matched validated comparisons; vertex-cost range 0 to 1.
- CPS-3F: 8/8 matched validated comparisons; vertex-cost range 0 to 0.

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

## Descriptive runtime and memory observations

These ranges summarize distinct configurations, not timing repetitions. They exclude reused observations and resource exits from successful-runtime ranges.

| Method | Newly validated n | Solver seconds, min-max | New-attempt peak RSS, min-max MiB |
|---|---:|---:|---:|
| Independent continuous | 8 | 0.009002 to 0.012559 | 210.78 to 227.71 |
| Independent discrete | 7 | 0.005237 to 0.012077 | 193.85 to 214.66 |
| CPS-2F | 5 | 0.375706 to 20.406551 | 222.39 to 346.34 |
| CPS-3F | 8 | 0.029120 to 0.346828 | 194.26 to 216.56 |

Resource-limited observations (all objectives missing; times are time to termination):

| Pair | Alpha | Method | Origin | Limit | Solver seconds | States | Transitions |
|---|---:|---|---|---|---:|---:|---:|
| AL031854__AL011991 | 1 | CPS-2F | new | transition_budget_limit | 26.114 | 16327 | 20000001 |
| 1o7j.a__1d9q.d | 1 | CPS-2F | new | state_budget_limit | 11.773 | 100000 | 9667290 |
| 1o7j.a__1hfj.c | 1 | CPS-2F | reused | transition_budget_limit | 50.765 | 57332 | 20000001 |

CPS-2F used one additional vertex in k only for EP021949__EP031949 at alpha=0.5 (6 versus 5); the other four matched validated continuous comparisons had zero vertex cost. CPS-3F had zero cost in all eight matched discrete comparisons. CPS-2F produced smaller k than CPS-3F on all five settings where both validated, but their fidelity constraints differ and CPS-2F had three resource exits. These observations do not establish general speed or compression guarantees.

