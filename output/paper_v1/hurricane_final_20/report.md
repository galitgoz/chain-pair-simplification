# Stratified final hurricane cohort

**COMPLETED: 80/80 scheduled configurations observed; 79 validated, one transition-budget exit, zero never launched.**

The resumed execution retained four previous observations and launched exactly 76 new attempts. Execution-session wall time was 2774.748 seconds (46.25 minutes), including readiness, warm-up, solver and validation phases. All readiness checks passed in this session: continuous and cumulative host waiting were zero. No resource-limited attempt was retried. See [bounded-wait implementation and checks](provenance/bounded_wait_notes.md).

## Scope and selection
Seed 20260916; 2209 eligible supplied unordered pairs; full track sizes [25, 98]. Quartile cutpoints: [32.0, 37.0, 44.0]; ratio split at 1.5. Allocation: {'Q1_balanced': 4, 'Q1_unequal': 0, 'Q2_balanced': 4, 'Q2_unequal': 0, 'Q3_balanced': 3, 'Q3_unequal': 3, 'Q4_balanced': 3, 'Q4_unequal': 3}. All candidates, exclusions, empty strata and deterministic redistribution are saved. No upper-size cap. The selected tracks span 25–69 vertices. Prior storm reuses: 0; within-cohort reuses: 0.

Preference-constrained stratified random sampling is not simple random sampling. Unweighted percentages describe this cohort, not the hurricane archive. Alpha=.5 was chosen using earlier pilot evidence. The sampling frame is the supplied pairs.csv, not every combinatorial pairing of all storms.

## Execution
80/80 observed; 79 validated; 1 resource-limited; 0 validation-incomplete; 0 pending. Current phase: completed. Pause reason: none.
Resource exits are completed observations and never retried. Raw checkpoints precede external validation. Independent coupling is descriptive. Validation checks feasibility, not independent optimality. Common external 60-second/512-MiB budgets coexist with disclosed method-specific internal work limits. No optimization algorithm, endpoint rule, tolerance, warmup, cache policy or thread setting was changed.

## Matched quality
Continuous: 17/19 matched objectives equal; 2 higher coupled objectives; 13 equal-k coupling repairs.
Discrete: 20/20 matched objectives equal; 0 higher coupled objectives.
Continuous/discrete comparisons have different fidelity feasible sets. Unvalidated quality stays missing. Timing and memory are single-attempt observations; no scaling law or significance claim is supported.

## Evidence layers and files
This new cohort is reported separately in results_with_strata.csv, tables/, graph_bounds.csv, and figures/. Earlier cohorts remain in ../coverage_extension/combined_results.csv with explicit cohort labels. The separate 120-repetition study is unchanged. Original geometric illustrations are retained; the gallery covers every saved attempt rather than selecting only attractive outcomes.
[Gallery](../coverage_extension/attempt_gallery/index.html), [metrics](../coverage_extension/attempt_gallery/per_attempt_metrics.csv). Auxiliary locations are matching positions, not selectable output vertices. Transition checks are not unique configuration edges; unique edges and DP labels remain unavailable. Finite bounds count like objects; interrupted counts are observed before termination, not progress percentages.

## Resume / stopping point
If pending slots remain: `.venv/Scripts/python.exe notebooks_paper/coverage_extension.py execute --cohort hurricane_final_20` after readiness passes. Recorded attempts are skipped. No further datasets, thresholds, budget increases or timing rounds are authorized.
