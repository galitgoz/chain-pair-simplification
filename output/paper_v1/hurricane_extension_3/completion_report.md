# Completed continuation and three-pair extension

The previously pending 27 eligible coverage-extension configurations completed once each. The coverage extension now has 48/48 validated eligible observations; eight preparation-blocked slots remain unexecuted. Historical resource exits were not retried.

The three-pair extension has 12 recorded single observations, 12 validated, 0 resource exits. No timing repetitions were added.

| Pair | nA/nB | Independent continuous k | CPS-2F k | Independent discrete k | CPS-3F k |
|---|---:|---:|---:|---:|---:|
| EP172000__EP061970 | 32/29 | 5 | 5 | 30 | 30 |
| EP111995__EP032014 | 41/41 | 8 | 8 | 39 | 39 |
| EP061986__EP131991 | 38/70 | 9 | 9 | 67 | 67 |

Combined primary coverage: 100/108 intended configurations observed, 99 validated, 1 resource exits, 8 preparation-blocked. Repetitions and historical diagnostics are excluded from this denominator.

## Scientific comparison
- Continuous fidelity: 21/24 validated matched objectives are equal; 3 coupled objectives are larger.
- Discrete fidelity: 25/25 validated matched objectives are equal; 0 coupled objectives are larger.
- Independent coupling violations are permitted. Continuous/discrete methods have different feasible sets. Validation establishes feasibility, not independent proof of optimality.
- Shared protein structure, subsampling, differing historical preparation regimes, pilot-informed alpha, and single new observations limit generalization. Resource exits do not establish infeasibility.

## Evidence and presentation
- [Gallery](../coverage_extension/attempt_gallery/index.html): 306 attempt records, 169 canonical PNG/PDF figures. Includes legacy integrated validation, phase-separated validation and timing repetitions with distinct provenance; these counts are not unique-configuration coverage.
- [Metrics](../coverage_extension/attempt_gallery/per_attempt_metrics.csv), [availability](../coverage_extension/attempt_gallery/metric_availability.csv), and [finite bounds](../coverage_extension/attempt_gallery/per_attempt_bounds.csv). Exact thresholds, selected indices, hashes, recovery links, graph counters and phase/resource records are retained where available.
- Unique configuration-graph edge and DP-label counts are unavailable. Missing historical measurements remain missing; zeros and not-applicable fields are distinguished. Transition checks are not unique edges. Bound occupancy is not search progress.
- Figures were generated after timed execution. Protein views are static 3D; hurricanes retain the common projected frame. Unvalidated attempts show inputs only. Identical timing-repeat geometries may share a figure while their measurements remain separate.
- [Combined results](../coverage_extension/combined_results.csv), [chapter](../experimental_chapter/experiments.tex), and existing saved-results notebook were updated. The LaTeX draft passed structural checks; no compiler was available.

## Implementation and stopping point
Only orchestration output/schedule parameterization and saved-evidence reporting were extended. Worker, validator and optimization algorithms are unchanged. The pre-extension orchestrator source is retained as provenance text. One attempt per new slot, no budget increases, no certificate shortcuts, no new timing rounds. No experiment worker remains active.
