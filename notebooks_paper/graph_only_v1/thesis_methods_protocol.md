# Thesis methods insert — configuration-graph evaluation

This evaluation compares optimal independent continuous simplification, optimal independent discrete simplification, CPS-2F, and CPS-3F. All methods use identical input curves, thresholds, and anchored endpoints. The objective is min max(|A′|, |B′|); output lengths may differ and a discrete coupling may repeat indices. Auxiliary positions are not selectable output vertices. Continuous fidelity is global.

The current software contains an independent-solution certificate shortcut and an alternative-subsequence certificate search. **Both are disabled in this evaluation by explicitly passing `certificate=False` to every CPS-2F and CPS-3F call, including correctness checks and warm-up.** The alternative-subsequence search is neither a method studied in this thesis nor a claimed contribution or local vertex-replacement method.

Successful CPS results are accepted only when their raw status is `optimal_configuration_graph`, both chain graphs were constructed, and execution entered the configuration-graph recurrence. Independently feasible solutions cannot bypass this execution. The existing implementation still computes an independent lower bound during preprocessing; this work is included in total solver time and is not mislabelled as graph-phase time. Separate independent baseline runs are reported independently.

Endpoint rejection, preprocessing resource exits, graph-construction resource exits, and recurrence resource exits are recorded separately. A resource limit is not mathematical infeasibility. Timing is diagnostic and includes the existing preprocessing; phase benchmarking is deferred. Warm-up and process startup are recorded separately.

Previous certificate-enabled experiments are preserved but are incompatible with this graph-only protocol. Certificate-labelled successes are explicitly excluded and are never relabelled as graph results. The new pilot reuses the frozen two protein and two hurricane pairs, with α∈{0.5,1}, δ₁=αs_A, δ₂=αs_B, δ₃=d_dF(A,B), and no additional subsampling or alignment. Existing protein w=16 curves are not full-resolution backbones. Single executions are not a runtime study.

This is a methods/protocol insert only. No editable thesis draft was found in the project; reference PDFs and previous notebooks remain unchanged. No experiments-section conclusions or later-stage benchmarks are introduced.
