# Graph-only Stage 1 — final report

## Protocol correction
All CPS calls, including warm-up and correctness checks, explicitly pass certificate=False. The adapter rejects True. The independent early return and alternative-subsequence search are disabled; no shortcut functions were deleted. Independent continuous and discrete baselines are separate runs. Existing lower-bound preprocessing remains inside solve_cps and is included in total time, never used to bypass graph execution.

The same frozen cohort and numerical thresholds were reused: two protein pairs at existing w=16, two hurricane pairs, alpha=0.5 and 1, four methods, 32 attempts. No larger or later-stage experiment was started.

## Verification
- 100/100 final small-input checks PASS, including exhaustive enumeration, independently feasible examples, and the 3-vertex example eligible for alternative certification.
- The guarded alternative-subsequence function was called zero times.
- An initial observer revision failed 16 assertions because it counted baseline single-curve graphs as joint construction. The observer was corrected before the pilot. The initial check log is preserved; this was an instrumentation error, not an optimization mismatch.
- The observer now counts only direct chain-graph calls from solve_cps and marks product-graph entry at its first heap pop. Successful joint rows additionally require raw status optimal_configuration_graph. Baseline graph instrumentation is not a performance metric.
- 8/16 CPS attempts returned verified optimal graph solutions. All 8 successful CPS rows actually entered the product graph and completed both chain graphs. No certificate success is included.
- 22/32 total verified solutions, 10/32 missing. Status counts: {'optimal': 22, 'other_resource_limit': 7, 'time_limit': 3}.
- Joint route outcomes: {('optimal', 'configuration_graph'): 8, ('other_resource_limit', 'configuration_graph'): 5, ('other_resource_limit', 'graph_construction'): 2, ('time_limit', 'configuration_graph'): 1}.

## Preservation and exclusion
Six previous successful CPS certificate rows are listed in excluded_certificate_results.csv. All 32 previous pilot rows are tracked in prior_result_compatibility.csv and none is reused as a new timing measurement. Original notebooks, results, and root algorithm modules are unchanged; preservation audit passed for 1225 files. Read-only process inspection found the pre-existing Jupyter kernel, which was neither attached to nor stopped.

## Call sites and thesis text
See callsite_audit.csv for every executable solve/solve_cps call in the revised modules. The new call path replaces the previous certificate-default pilot adapter: graph_core.solve -> solve_cps(certificate=False), pilot_worker measured call, graph_core warm-up, graph_checks tests, and shortcut_probe all use explicit False. Prior experiment modules are preserved as historical software.

The methods/protocol insert is notebooks_paper/graph_only_v1/thesis_methods_protocol.md. No editable thesis draft was found, so reference PDFs were preserved. The alternative-subsequence search is not a studied thesis method, contribution, or local replacement algorithm. No paper conclusions are drafted.

## Resources and limitations
One attempt per configuration, diagnostic timings only. Worker wall-time sum 209.701 s; maximum sampled worker-tree RSS 214.95 MiB. Limits: 20 s execution watchdog after warm-up, 90 s startup, 512 MiB RSS; graph caps 8 s soft time, 30,000 states, 500,000 product transitions and 200,000 chain edges. Watchdog polling can overshoot under OS scheduling delays. Warm-up, preprocessing, graph construction, and graph resource exits are distinguished by route/phase; missing outputs are not infeasible. No per-phase timing benchmark is claimed.

Protein curves are already subsampled, share 1o7j.a, and retain supplied gap bridges; 1d9q.d has the existing documented count discrepancy (325 supplied versus paper length 297). Hurricane inputs retain the common basin frame and lack timestamps. These remain pipeline examples, not evidence of topological preservation or broad performance conclusions.

## Outputs and next step
New notebooks 00–03 are under notebooks_paper/graph_only_v1. Tables, figures, logs, exclusions and validation are under output/paper_v1/graph_only_v1. See pilot_summary.csv, pilot_results.csv, graph_route_summary.csv and validation.csv.

Stage 1 stops here. Any retry of the 10 missing configurations, change of resource caps, or larger experiment requires a further request. Notebooks 04/05 and later stages remain unexecuted.
