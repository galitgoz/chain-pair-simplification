# Graph-only v2 — repair and diagnostic stage completed

## Reconciled and preserved
The current v1 files confirm 100 passing correctness checks, 22 verified solutions out of 32, and 8 successful CPS graph solutions. The 10 missing results comprise 5 product-transition limits, 2 single-chain-edge limits, and 3 external timeouts. No newer graph-result version was found. The v1 timeouts covered solving plus external validation/measurements; their original limiting phase cannot be recovered from those files.

The 22 verified rows are preserved in inherited_verified_results.csv with source hashes and original timings, and are **not new executions or timing repetitions**. Existing code, notebooks and results passed a 1433-file preservation audit. The pre-existing Jupyter kernel was inspected read-only and not stopped. No optimization implementation was changed.

## Repair
Each attempt has a UUID directory. Limits and code/input/configuration hashes are written before launch. Phases are input loading, warm-up, solving, atomic raw checkpoint, external validation, optional measurements, and reporting. The raw checkpoint uses a temporary file, flush/fsync, then atomic replace, before external validation or distance estimates.

The existing CPS `finish()` performs endpoint/order checks, fidelity tests, and discrete coupling reconstruction/validation before returning. This remains inside solver time; it has not been removed or relabelled as external validation. Independent discrete does not call that finish function. External validation establishes feasibility and count consistency, not optimality by itself; the solver optimality claim is retained separately. The small-instance exhaustive tests check objectives.

All CPS calls explicitly use certificate=False. Compilation warm-up has no internal time cap, is performed on synthetic curves, and has its own 120 s external budget. Startup: 60 s; input loading: 15 s; solver: 60 s external, 55 s internal for CPS (soft, checked by the existing algorithm); checkpoint: 5 s; validation: 30 s; optional measurements: 30 s. Independent discrete has no internal time parameter; its existing single-curve edge default of 2,000,000 is unchanged. RSS cap: 512 MiB over the worker tree, sampled approximately, not exact algorithm memory. Watchdogs can overshoot due to OS scheduling.

Regression outcomes: 21/21 passed, including 16 small exhaustive comparisons, checkpoint survival after validation failure, validation timeout and optional measurement failure, plus an independent baseline whose coupling exceeds delta3 without failing its own constraints. Synthetic tests are not pilot repetitions.

## Exactly three new attempts
All three warm-ups completed. See diagnostic_summary.csv for the three-row table and attempts/ for committed configurations, phase/resource logs, raw/exit checkpoints and separate validation records.

- **A**, AL031854__AL011991, alpha=0.5, CPS-2F: still stopped at the increased product-transition cap, after 2,000,001 examined transitions; 5,452 discovered states. Solver 2.109743 s. No solution returned or raw solution checkpoint; validation not run. Solver exit checkpoint preserved. Other graph-work caps unchanged.
- **B**, 1o7j.a__1hfj.c, alpha=1, CPS-2F: the 800,000 single-chain-edge allowance allowed both chain graphs to complete; the next bottleneck was the unchanged 500,000 product-transition cap (stopped at 500,001, with 17,567 discovered states). Solver 1.838192 s. No solution returned; exit checkpoint preserved; validation not run. No automatic second escalation.
- **C**, 1o7j.a__1d9q.d, alpha=0.5, Independent discrete: returned kA=20, kB=19, k=20. Raw output was saved before external validation. Validation passed. Solver 0.006537 s; validation 0.874883 s; optional measurements 0.037013 s. This does not retrospectively identify the phase responsible for the v1 timeout.

Total diagnostic worker wall time: 8.837 s. Maximum sampled tree RSS: 221.68 MiB. No diagnostic hit the time or memory cap. Limits ended A and B; neither outcome is mathematical infeasibility.

## Reference comparison
A's saved certificate-based k=5 witness was revalidated against the identical input hash and thresholds; the independent continuous k=5 lower bound was checked as saved reference evidence. It is a reference optimum, not a graph-run result. B's saved independently optimal continuous k=5 pair was revalidated and satisfies delta3, supporting the same reference optimum. Neither diagnostic returned a graph solution, so no new graph objective can yet be compared with 5. No identical-configuration verified reference was located for C in the prior frozen-protocol versions. No discrepancy or infeasibility report occurred. Feasibility validation alone is not the argument for either reference's optimality.

## Coverage and recommendation — stop here
Inherited verified results: 22. New diagnostic attempts: 3. Newly validated outputs: 1. Combined coverage: 23/32; remaining missing configurations: 9. See coverage_with_provenance.csv; inherited rows retain their original provenance and are not timing repetitions.

Use this phase-separated pipeline for any future timeout repair; C now completes and preserves its output, but the old timeout cause remains unidentified. A needs a transition-work investigation or a separately authorized count-only increase; more wall time alone would not bypass its count cap. B's first construction bottleneck was cleared, but a future controlled retry would need to address its transition cap, while monitoring states and memory. Do not generalize these two runs into a runtime or complexity conclusion. For the other missing configurations, identify the resource category before proposing a targeted retry; do not rerun all nine automatically.

No further attempts, timing repetitions, subsampling, threshold changes, larger cohort, greedy preprocessing or paper conclusions were run. This stage is stopped after A, B and C.
