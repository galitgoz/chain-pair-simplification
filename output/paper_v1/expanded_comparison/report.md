# Primary expanded comparison: alpha=0.5

## Coverage and provenance
Ten frozen pairs x four methods = **40 primary configurations**. alpha=0.5 was
chosen using earlier pilot evidence; this is an exploratory, evidence-informed
comparison, not a prespecified independent confirmatory study.

| Method | Validated / planned | Resource exits | Unexecuted |
|---|---:|---:|---:|
| Independent continuous | 10/10 | 0 | 0 |
| Independent discrete | 10/10 | 0 | 0 |
| CPS-2F | 9/10 | 1 | 0 |
| CPS-3F | 10/10 | 0 | 0 |

- Inherited observations: 16 (original solver/validation paths retained).
- New expansion slots executed: 24/24, including the previously executed slot 33.
- Solvers returned an optimal-status output: 39/40.
- Successfully validated: 39/40.
- Resource-limited observations without validated output: 1.
- Unexecuted: 0; termination categories: {"transition_budget_limit": 1, "validated": 39}.

All alpha=1 artifacts are preserved separately in alpha1_pilot_evidence.csv and
all_alpha_results.csv. No pending alpha=1 configuration was launched. Original
configuration_manifest.json is unchanged; execution_manifest_alpha05.json is the
revised execution plan. presentation_before_alpha05 preserves the prior presentation.
results.csv is the primary 40-row table, including failures and missing outputs.

## Validation recovery and targeted fix
Slot 33 returned in 0.0478647 s and saved 16 A indices and 20 B indices. The original
validation exceeded 30 s before optional measurements. Its frozen protein arrays
are Fortran-contiguous; warmup compiled the C/C layout of continuous_decision, and
the saved cache lacked the F/C specialization required by the original validator.
This exposed compilation inside validation. Contention was also observed later;
the original log does not identify the exact time spent compiling versus scheduled
CPU work, so a unique historical cause cannot be proved.

The targeted fix copies validation inputs to C-contiguous layout without changing
coordinate values, dtype, fidelity predicate, tolerances, endpoints, ordering or
membership checks. Five focused checks compare layouts and reject invalid fidelity,
endpoints, ordering and membership. Independent coupling is descriptive, not required.
Achieved continuous-distance searches remain in optional measurements outside validation.
No validation budget was increased. Validation-only replay passed in
0.003733 s; kA=16,
kB=20, k=20. No optimization was rerun for this output.
The original timeout remains unchanged, and validation_recovery contains the separate
replay, source checkpoint hash, worker code provenance, phase and resource logs.
The row retains original solver/warmup/total-process measurements; recovery timing
is a separate measurement, not a second solver repetition.

Slot34 subsequently saved an optimal-status independent-discrete output but hit
WinError5 while atomically replacing progress.json immediately before validation.
The orchestrator was paused after the current worker finished; no worker or unrelated
process was killed for this repair. A bounded 0.7-second retry of the I/O rename
(not a solver retry) handles transient Windows sharing denials within unchanged
phase budgets. The exact competing reader is unknown. A focused injected-denial
test passed. Slot34 validation-only recovery passed in approximately 0.0030 seconds;
its original solver output and error remain preserved. Outcome aggregation now
exposes worker software errors even when a raw optimal-status checkpoint exists.
See orchestration_repair.json, atomic_replacement_checks.json and validation_recovery_slot34/.

## Inputs, methods and resource policy
All ten input pairs and exact thresholds were frozen before expansion. Every protein
pair shares 1o7j and has 22 vertices per curve after existing w=16 sampling of supplied
prepared backbones. These pairs establish neither protein input-size scaling nor
independent biological replication. Newly selected hurricane tracks retain every
supplied point (including 25/82, 62/81 and 30/29 lengths), without alignment or cropping.
Shared basin coordinates are preserved; coordinates alone do not imply synchronized time.
Pair provenance and input hashes are in input_manifest.json. Original pilot pairs and
newly selected pairs are identified in results.csv and the size figure.

delta1=0.5*mean edge length(A); delta2=0.5*mean edge length(B);
delta3=unrounded discrete Frechet distance(A,B). Identical inputs, endpoints and
thresholds apply across methods. Objective: min max(kA,kB). Independent continuous
and discrete baselines enforce their respective global fidelity only. CPS-2F enforces
continuous fidelity plus discrete output coupling; CPS-3F enforces discrete fidelity
plus discrete coupling. All outputs use original vertices; unequal lengths and repeated
coupling indices are allowed. certificate=False disables both shortcuts. No greedy
processing, vertex replacement, threshold tuning or reference-guided stopping was added.

External solver60s, internal CPS55s, sampled tree RSS512MiB; CPS100000 states,
20000000 candidate checks and 800000 main single-graph edges. Independent discrete
and CPS-3F independent preprocessing retain the 2000000-edge default. Independent
continuous has no equivalent graph-work cap. These are shared external limits with
disclosed method-specific internal limits, not equal computational work guarantees.
Separate unchanged phase budgets: startup60, input15, warmup120, checkpoint5,
validation30, optional measurements30, reporting10 seconds. Internal soft time is checked
after state expansion; external monitoring is the backstop. No automatic retries.

## Observed size and coupling comparisons
- CPS-2F: 9/10 validated matched comparisons; vertex cost range 0 to 1; 8 equal-size outcomes.
- CPS-3F: 10/10 validated matched comparisons; vertex cost range 0 to 0; 10 equal-size outcomes.

Only matched validated outputs enter vertex-cost comparisons. retained_A, retained_B
and retained_total record retained fractions. CPS-2F is compared against independent
continuous; CPS-3F against independent discrete. Direct CPS-2F/CPS-3F comparisons
change fidelity constraints as well as the algorithm. Missing outputs are not zero.
In the fixed AL031854/AL011991 illustration both independent continuous and CPS-2F
have k=5, but coupling is 725.571465 versus 426.184218 km for a 426.184218 km budget.
This demonstrates equal size with different coupling feasibility, not universal
zero vertex cost or topological preservation.

## Graph bounds and computational cost
- CPS-2F: 6/10 configurations have the tighter active-state product bound; observed occupancy ranges from 13.55% to 79.07%, with 1 interrupted search(es) included explicitly.
- CPS-3F: 6/10 configurations have the tighter active-state product bound; observed occupancy ranges from 14.97% to 77.34%, with 0 interrupted search(es) included explicitly.
- AL072006__AL011900, CPS-2F: transition_budget_limit; 20,000,001 transition checks, 14,181 discovered states, 27.975 s to termination, 248.96 MiB sampled peak RSS; no returned objective.

graph_bounds.md cites thesis Observation4 (printed p18/PDF26) and section5.1
(printed p26/PDF34), and manuscript sections3.1/4.1/4.4 (PDF7-8/11-12/15-16).
For CPS-2F, auxiliary bound U(n)=2n(n-1); finite potential-state bound
n*m*(n+U(n))*(m+U(m)). CPS-3F has no auxiliary points and bound n^2*m^2.
Actual extended-size bound is n*m*(n+aA)*(m+aB); completed active single-graph
counts give H_A*H_B. Neither product enumerates coupling-valid states. Source and
terminal are included. Auxiliary matching locations are not selectable output vertices.
Big-O expressions are not assigned a hidden constant of one.

graph_bounds.csv separates discovered states, dense single-graph slots, active states,
single-graph edges including waits, auxiliary locations and candidate transition checks.
Single-graph edge bound H(H+1)/2 counts the same objects (unique forward edges plus
H waits). Candidate-check bound E_A*E_B-H_A*H_B is derived separately; the implicit
product graph is not materialized, and candidate checks are not relabelled as unique
stored edges. Unsupported product-edge and DP-label observations remain missing.
Percentages compare matching counting quantities and are never completion estimates.
Interrupted counts are observed before termination. Historical missing counters are
not filled by rerunning solvers. No packedness or general scaling law is inferred.

Solver time includes existing internal fidelity/coupling checks. External validation
has its own timer; recovery is separately identified. Successful times, solver resource
exits and returned-but-unvalidated outputs use distinct markers. Inherited observations
are hollow. RSS is sampled whole-process-tree memory, including warmup, not exact
algorithm memory. Boundary graph counters and five-second supervisor status updates
add small unquantified overhead, so these single attempts are not a timing benchmark.
Feasibility validation alone does not establish optimality; solver optimality claims
remain separate, supported by prior small-instance exhaustive checks.

## Live progress and files
progress.json is atomically refreshed approximately every five seconds outside solver
loops; progress_history.jsonl retains updates. New slots finished /24 counts recorded
terminal outcomes, not just successes. Reused observations are separate. Graph counters
carry timestamps; budget-used percentages never mean search completion. The saved-only
presentation notebook displays progress and results without starting optimization.

[Results](results.csv), [completion](completion.csv), [matched comparisons](matched_comparisons.csv),
[graph counts and formulas](graph_bounds.csv), [compact state table](graph_state_summary.csv),
[captions](figure_captions.md), figures/*.pdf and *.png, and original per-attempt checkpoints.
Algorithms and scientific parameters are unchanged. common.py has the layout-only
validation fix; worker.py supports validation replay; supervisor.py has an optional
five-second observation callback; expanded_comparison.py manages the revised manifest
and reporting. No duplicate infrastructure or new versioned code directory was created.

The batch was resumed after an orchestration pause. Historical progress entries named
total_batch_elapsed_seconds measured their current execution segment. The final progress
record separates that segment from total wall time since manifest freeze (including
pauses and recovery); timing_definitions.json records this clarification. No solver
timers or historical progress records were rewritten.
