## completion

Validated / intended configurations by domain and method. Colored stacks retain each distinct termination reason; white-grey denotes not executed.

## solution_size

Validated outputs only on all 10 matched pair/alpha slots. Black crosses mark independent outputs exceeding delta3 (not baseline failures). Missing outputs are omitted, never zero; see completion table. Shading identifies newly selected pairs; alpha was selected using pilot evidence.

## solver_time

Circles: validated solutions; X: solver resource exits (time to termination); square: solver returned but validation incomplete (time to raw return, not validation termination). Hollow markers: inherited observations. One observation/configuration; added boundary instrumentation precludes strict timing equivalence.

## configuration_states

Input-only finite bound: CPS-2F n*m*(n+2n(n-1))*(m+2m(m-1)); CPS-3F n^2*m^2. Includes all index combinations and endpoints, not just valid states. X indicates interrupted observations. Dashed line is equality, not a fitted scaling law.

## state_utilization

Left: S_A*S_B from actual extended lengths (CPS-2F only; CPS-3F repeats its input bound and is omitted). Right: H_A*H_B from saved active counts. Missing historical counters stay missing. Percentages measure occupancy of bounds, never progress.

## single_curve_sizes

Each point represents curve A or B for a CPS configuration; thresholds are not independent repetitions. No auxiliary locations in CPS-3F. Graph edges include self/wait entries. Only saved counters are shown; no implied linear growth or packedness claim.

## single_curve_edge_bounds

Completed single-curve graphs only. H is the active vertex count; unique forward edges plus one wait at each active vertex obey H(H+1)/2. Circles: completed product search; X: subsequently interrupted product search. These are single-graph edges, not product-transition checks. Historical missing H counts remain missing.

## algorithm_diagram

Anchored global-fidelity configuration construction and coupled output. Independent baselines solve the single-curve problems without enforcing output coupling.

## atlantic_coupling_example

Pre-specified illustration AL031854 / AL011991, alpha=0.5, full supplied 25-point tracks. Faint lines: inputs; bold lines: simplified curves; dotted lines: saved discrete coupling. Equal k does not imply the same coupling feasibility. Coordinates do not imply temporal synchronization.

## memory

Sampled whole-worker-tree RSS, including preparation and warm-up; not exact algorithm memory and not guaranteed peak capture. Dashed line: 512 MiB cap. Circles validated, X unsuccessful.