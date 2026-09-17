# Experimental results for thesis and manuscript

## Status
Original cohort: 39/40 validated configurations. Fresh timing study: **120/120 recorded**, 117 validated, 3 resource-limited, 0 pending. Three rounds are scheduled; inherited observations are not repetitions. This draft is provisional while timing slots remain pending.

## Supported quality findings (original cohort)
- Continuous-fidelity comparison: 8/9 matched cases have equal k; 1 has a larger joint k.
- Discrete-fidelity comparison: all 10 matched cases have equal k.
- 3 continuous-fidelity cases repair coupling without increasing k. Independent coupling violations are permitted by their problem definition.
- The Atlantic illustrative example has equal k=5 with different coupling feasibility. EP111954/EP021984 illustrates CPS-2F k=16 versus CPS-3F k=79 under different fidelity constraints.
- 0 unexpected objective differences among recorded validated timing outputs. Pending slots support no timing conclusion.

## Interpretation and limits
The state figure distinguishes input-only finite bounds, active single-graph products and discovered states. Missing historical counters remain missing. Occupancy is not search completion. Transition checks are not unique edges. Alpha was chosen using pilot evidence; all protein inputs have 22 vertices and share 1o7j. These data do not establish general scaling, independent biological replication, packedness, topological preservation or statistical significance.

## Files and traceability
- `main.tex` includes the new `experiments.tex`; the user's original file is untouched.
- `tables/timing_observations.csv`: all 120 planned repetitions, individual phases, raw-result paths and outcomes.
- `tables/timing_summary.csv`: 40 configurations, all three statuses, success-only median/range/count, resource-exit times separately.
- `tables/claim_traceability.csv`: generated numerical chapter claims and source filters; table entries derive directly from CSVs.
- `tables/original_cohort.csv` and `tables/graph_bounds.csv`: preserved primary evidence, copied for this report.
- `figures/`: PDF and PNG; `figure_captions.md`: self-contained captions.
- `provenance/`: schedule/source hashes and build audit.

## Execution policy
The frozen schedule interleaves methods across three rounds. One worker at a time, unchanged warm-up/cache/thread policy, certificate=False, no algorithm changes, no replacement repetitions. Shared external time/memory caps have disclosed method-specific work limits. The schedule is resumable and completed slots are never duplicated.

## Recorded batch state
Phase: completed. Last progress update: 2026-09-16T17:37:13.918837+00:00. None

The readiness records preserve temporary PyCharm contention and all pauses. No unrelated processes were terminated and no system settings were changed.

## LaTeX checks
File inclusions, references and environment structure passed. No LaTeX compiler was found on PATH; the draft has not been compiled. PDF figures are rendered separately for visual inspection.

## Coverage extension update

178/180 recorded combined single observations validated; 2 resource exits; 8 preparation-blocked slots. See [extension report](../coverage_extension/report.md) for all 188 intended slots and preparation caveats.

## LaTeX checks
File inclusions, references and environment structure passed. No LaTeX compiler was found on PATH; the draft has not been compiled. PDF figures are rendered separately for visual inspection.

## LaTeX checks
File inclusions, references and environment structure passed. No LaTeX compiler was found on PATH; the draft has not been compiled. PDF figures are rendered separately for visual inspection.
