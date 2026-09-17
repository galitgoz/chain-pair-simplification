## coupling_vertex_cost

Does enforcing coupling require more vertices? Original alpha=0.5 cohort: output discrete coupling divided by the shared threshold (top), and matched-fidelity objective difference (bottom). Dashed line is the coupling bound. A violation is permitted for independent baselines. Stars identify equal-objective coupling repairs. The missing CPS-2F output is a resource exit, not zero cost. Pair aliases and exact thresholds are in the input table.

## graph_size_levels

How large is the explored state set relative to finite construction bounds? Original cohort only. D counts discovered configurations, including the source. The input-only bound follows sphere–edge intersection counting; H_A H_B is the product of active single-curve vertices after anchored pruning. Missing historical H counters remain missing. Crosses mark counts observed before termination. Neither completed search nor bound occupancy implies enumeration of all valid states; percentages do not measure progress. Formulas and assumptions are given in the graph-size section and graph_bounds.csv.

## completion

What completed under the fixed resource policy? Original configurations and fresh scheduled repetitions have separate denominators. Green means a saved output passed external feasibility validation; orange denotes resource termination, never mathematical infeasibility. Pending repetitions have no measured outcome. Exact termination reasons are retained in timing_observations.csv and the completion table.

## timing_stability

How stable are repeated outcomes and solver times? Each point is one scheduled fresh repetition, offset by round. Bars show median and min–max only among validated successes, with counts and all outcomes in timing_summary.csv. Crosses are time to resource termination and are excluded from success statistics. Blank locations are unexecuted. Solver time retains internal checks; warm-up and external validation are separate. Repetitions are not independent input instances.

## compression_cost

What compression is obtained at the observed computational cost? Only validated fresh repetitions appear in the quality panel. Protein and hurricane markers differ; colors identify methods. Continuous fidelity and discrete fidelity define different feasible sets, so shorter CPS-2F outputs do not establish superiority for an identical optimization problem. Resource exits appear separately because no solution quality is available. No historical runtime is pooled into this timing plot.

## work_indicators

Which recorded work indicators accompany solver time? Fresh repetitions only; circles are validated successes and crosses are recorded unsuccessful terminations. Historical missing counters are not filled. Discovered states are not all valid states; transition checks are not unique graph edges or all machine operations. These are descriptive associations with no complexity-law fit. All protein inputs have the same sizes and share one structure.

## supplement_memory_phases

Supplement: sampled peak process-tree RSS and separate recorded times for fresh repetitions. Memory includes interpreter, compilation, allocator and monitoring-visible descendants, and is an approximate sampled peak, not exact algorithm memory. Phase values and unsuccessful outcomes remain individually available in the observation table; total process time is not solver time.

## geometric_examples

Illustrative selections, not representative averages. Top: Atlantic AL031854/AL011991 has equal k=5 for independent continuous and CPS-2F, but only the coupled output satisfies the coupling bound. Bottom: EP111954/EP021984 has CPS-2F k=16 and CPS-3F k=79, illustrating the change in feasible set between continuous and discrete fidelity. Pale curves are supplied inputs; strong curves and markers are selected original vertices. The dashed segment marks a bottleneck pair in the saved discrete coupling. Distances are rounded for display only; full-precision saved thresholds are unchanged. No independent storm alignment or temporal synchronization is inferred.
# Coverage extension figures

## coverage_completion

Combined single-observation coverage with denominators of 14 planned protein pairs and 10 hurricane pairs per method. Hatched slots are preparation-blocked, not experiments or infeasible instances. Resource limits retain their exact reasons. The original 40 observations retain their provenance; the 120 timing repetitions are not pooled into this figure.

## coverage_quality

Matched simplification sizes from validated outputs, with missing outcomes omitted rather than zero. Black crosses identify independent outputs that violate the coupling bound; this is permitted for those baselines. Hatched protein slots lack an eligible prepared input. CPS-2F and CPS-3F have different fidelity feasible sets. Pair aliases map to the combined input manifest.

## coverage_coupling_cost

Does coupling improve feasibility at a vertex cost? Normalized output coupling (top) and objective change against the independent baseline with matching fidelity (bottom). Stars mark equal-k coupling repairs. Missing points are missing comparisons, not zero cost. Independent coupling violations are not baseline failures. Both curves use identical saved thresholds within each pair.

## coverage_runtime_memory

Time and sampled process-tree memory versus input vertex count. Hollow circles preserve inherited single observations; filled circles are new validated observations; crosses are resource exits, whose times are time to termination. Repeated timing measurements remain separate. Protein preparation regimes differ, so these descriptive points do not identify a causal size effect or complexity law. RSS includes runtime overhead.

## graph_size_levels

Combined eligible single observations: input-only finite bound, active single-graph product H_A H_B, and discovered states D. Percentages name the denominator; neither ratio is progress. Crosses mark interrupted searches. Missing historical counters stay missing. Blocked inputs have no constructed graph and are excluded from this graph-count comparison. Finite bounds follow the saved derivation, not a hidden constant of one in asymptotic notation.

## coverage_single_graphs

Supplement: each point describes one single-curve graph (circles A, triangles B), not a pooled configuration count. Auxiliary matching locations are not selectable output vertices. CPS-3F has no auxiliaries; zero is shown on a symmetric-log axis. Edges include waits after anchored pruning; the construction budget applies earlier during enumeration. Missing counters are omitted; crosses identify counts from interrupted construction where available.