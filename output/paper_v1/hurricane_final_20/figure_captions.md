## final20_completion
All 80 scheduled configurations, including pending and resource outcomes. Left: 20 per method. Right: recorded allocation by length quartile and ratio; balanced means ratio <=1.5. Empty strata are shown. Completion is not validation success; resource exits are not infeasibility.

## final20_quality
Matched input pairs and thresholds, only externally validated outputs. Black crosses mark independent coupling violations, not independent-method failures. Missing solutions are omitted, never assigned zero. S aliases follow frozen selection order.

## coupling_vertex_cost
Does enforcing coupling require more vertices? Output discrete coupling divided by the shared threshold (top), and matched-fidelity objective difference (bottom). Dashed line is the coupling bound. Independent violations are permitted. Stars identify equal-objective coupling repairs. Missing comparisons are not zero cost; their outcome or pending status is in the result table. Pair aliases and exact thresholds are in the input manifest.

## graph_size_levels
How large is the explored state set relative to finite construction bounds? D counts discovered configurations, including the source. The input-only bound follows sphere–edge intersection counting; H_A H_B is the product of active single-curve vertices after anchored pruning. Missing counters remain missing. Crosses mark counts observed before termination. Neither completed search nor bound occupancy implies enumeration of all valid states; percentages do not measure progress. Formulas and assumptions are in graph_bounds.csv. Input-only bounds can be shown before execution; they are not observed graphs.

## final20_compression_cost
Continuous and discrete fidelity define different feasible sets. Left: retained fraction versus solver time for validated solutions only. Right: resource-exit times, with no invented output quality. These are single attempts, not repeated timing estimates.

## final20_work
Descriptive work indicators only: input size, discovered states, and candidate transition checks. Circles are validated solution times; crosses are unsuccessful termination times. Missing counters are omitted. Transition checks are not unique edges or all operations. No complexity law is fitted.

## final20_memory
Supplementary sampled process-tree RSS, including runtime overhead; not exact algorithm memory. Black crosses indicate resource exits. Missing measurements are not zero.