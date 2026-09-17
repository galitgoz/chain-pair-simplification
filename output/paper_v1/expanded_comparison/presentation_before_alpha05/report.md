# Expanded bounded CPS comparison

## Scope and completion
Ten fixed pairs, two alpha values, four methods: 80 configurations.
Inherited observations: 32; new attempts: 1;
validated: 29; resource-limited: 4; unexecuted: 47.
No automatic retries, budget increases, new thresholds or timing repetitions.
Inherited observations retain their original attempt paths and are not new measurements.

| Method | Validated / planned | Resource exits | Unexecuted |
|---|---:|---:|---:|
| Independent continuous | 8/20 | 1 | 11 |
| Independent discrete | 8/20 | 0 | 12 |
| CPS-2F | 5/20 | 3 | 12 |
| CPS-3F | 8/20 | 0 | 12 |

## Inputs and protocol
[Input manifest](input_manifest.csv), [exact configurations](configuration_manifest.json),
[selection rule](selection_rule.md), [candidate ranking evidence](candidate_selection.json).
Three new pairs per domain were selected before solving, using input geometry and
reuse penalties only. Protein sizes are all 22/22 after supplied w=16 sampling of
already prepared 325-point curves; these are not full-resolution backbones.
The common 1o7j.a anchor is unavoidable in the supplied aligned pair list. Structure
reuse means pairs are not fully independent. Hurricanes preserve complete supplied
tracks and basin coordinates, including unequal lengths. No temporal synchronization
is inferred from coordinates. Gaps/provenance remain those of the supplied datasets.

alpha is 0.5 or 1; delta1=alpha*mean edge length of A, delta2 analogously for B;
delta3 is the unrounded discrete Frechet distance of saved inputs. Identical thresholds,
input arrays, tolerances and anchored endpoints apply to all methods in a configuration.
All methods minimize max(kA,kB). Independent continuous/discrete methods impose their
respective GLOBAL fidelity constraints separately and do NOT enforce delta3. CPS-2F
has continuous fidelity and discrete output coupling; CPS-3F has discrete fidelity
and discrete output coupling. Original vertices alone are output candidates; unequal
output lengths and repeated coupling indices are allowed. Certificates are disabled.

Shared external solver allowance 60s and sampled tree RSS 512MiB. CPS soft time 55s,
100000 states, 20000000 transition checks and 800000 main single-graph edges.
Independent discrete and CPS-3F independent preprocessing retain their 2000000-edge
defaults. Independent continuous has no corresponding graph-work cap. Internal soft
time is checked after expanded states, not continuously; external monitoring remains
the backstop. Separate budgets: startup60, loading15, warmup120, checkpoint5,
validation30, optional30, reporting10 seconds. Warmup is unchanged and output-local
caches are fresh per new slot. These practical bounds do not guarantee completion.

## Observed compression and coupling
- CPS-2F: 5/20 validated matched comparisons; vertex cost range 0 to 1; 4 equal-size outcomes.
- CPS-3F: 8/20 validated matched comparisons; vertex cost range 0 to 0; 8 equal-size outcomes.
Only paired validated outputs support these comparisons; failures stay in the
denominator. Full retained fractions, selected sizes and coupling values are in
[results.csv](results.csv); [matched comparisons](matched_comparisons.csv) uses the
baseline with the same fidelity metric. Comparing CPS-2F with CPS-3F also changes
the fidelity constraint; size differences are not pure implementation speedups.
The pre-specified Atlantic alpha=0.5 illustration has independent continuous and
CPS-2F k=5: coupling 725.571465 km versus 426.184218 km with delta3=426.184218 km.
Thus joint feasibility can change even when optimum output size does not increase.
This illustration is not a general zero-cost claim or evidence of topology preservation.

## Computational cost and graph size
[Graph-bound table](graph_bounds.csv), [derivation and source references](graph_bounds.md).
[Compact state comparison](graph_state_summary.csv); [source files and hashes](source_references.json).
Observed discovered states are compared with explicit finite products of index sets,
not a big-O expression with constant one. Active states, dense slots, auxiliary
locations, edges including waits and candidate checks are distinct. Missing historical
counters stay missing. Interrupted counts are observed before termination; occupancy
percentages do not estimate search completion. DP labels and all coupling-valid states
were not separately counted. No packedness claim or general complexity law is inferred.
See figures/configuration_states and figures/state_utilization for conservative
input-only and tighter data-dependent bounds respectively.

## Validation and limitations
Raw outputs are atomically saved before external validation. CPS finish() already
performs internal fidelity/coupling checks; this work remains in solver time.
External validation checks original indices, order, endpoints, counts, applicable
fidelity, and discrete coupling; baseline coupling violations are informational.
Independent continuous also checks each returned simplification internally; this
check is included in solver time. Warmup compiles C-contiguous continuous-decision
inputs, whereas the frozen protein arrays are Fortran-contiguous. The stopped
slot's cache index contains only the C/C continuous-decision signature: validation
may compile an additional layout specialization. The existing warmup procedure
was preserved, not silently extended. Validation time is therefore not guaranteed
to exclude compilation. This and concurrent host activity are plausible contributors,
not a proven causal attribution for the timeout.
Validation establishes feasibility, not optimality; optimality is the solver's
claim supported by existing exhaustive small-instance checks, not independently
proved for each large input. Resource exits are not infeasibility. Original inputs
are feasible joint outputs because delta3 equals their discrete Frechet distance.
Runtime observations are single attempts with separate phase clocks; resource exits
are times to termination, never successful solution runtimes. RSS includes the worker
tree and startup/warmup, is sampled, and is not exact algorithm memory.

## Provenance and changes
Algorithms, worker, supervisor and warmup are unchanged. common.py adds only scalar
metadata at main single-graph completion; it does not instrument inner transition loops.
The added boundary observation may add small unquantified timing overhead relative to
inherited observations. expanded_comparison.py supplies selection, orchestration and
presentation using the existing runner/worker/validation. Hashes are in provenance.json
and each attempt configuration. instrumentation_checks.json records two small graph
regressions. This is not a repeated timing benchmark; additional thresholds or repeated
runs of a pair do not create independent input pairs. Manuscript experiments are not
used as verified evidence. No further experiments are authorized or launched here.

## Files
[Complete results](results.csv), [completion](completion.csv), [captions](figure_captions.md),
figures/*.pdf and *.png, raw per-attempt checkpoints/validation/phase/resource logs in slots/.

## Execution pause and limits on conclusions

The expansion is **not complete**. One new independent-continuous solver returned
and saved its raw output (protein 1o7j.a / 1qd1.b, alpha=0.5), but external validation
exceeded 30 seconds. The raw objective is unverified and excluded from all quality
comparisons. The next readiness check recorded sustained PyCharm activity:
603.8, 621.4, 243.7, 153.2, 157.3, 148.1, 135.3, 90.9, 42.0, 37.5 percent of one CPU
core. The runner stopped before launching slot 34. Host activity at that check does
not prove the cause of the earlier validation timeout. No process was terminated
other than the attempt worker tree at its validation budget; user kernels remain.
No new CPS graph was run, so graph figures currently describe inherited evidence
only, and tighter active-graph bounds lack observed counters. The 29 verified
outputs are all inherited. No expanded-cohort compression, timing, or graph-scaling
conclusion is justified yet. Resume only the 47 untouched slots after host readiness;
do not retry the completed resource-limited observation automatically.
{
  "configuration_id": "protein:1o7j.a__1qd1.b:a0.5|Independent discrete",
  "reason": "host_readiness",
  "utc": "2026-09-16T12:58:30.860845+00:00",
  "unix": 1789563510.8608575,
  "monotonic": 141992.9085607
}
