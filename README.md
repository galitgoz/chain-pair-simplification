# Chain pair simplification

## Paper results: start here

This repository compares **CPS-2F** (continuous input fidelity, discrete output coupling), **CPS-3F** (discrete fidelity and coupling), and independent continuous/discrete baselines. The objective is `k = max(kA, kB)`.

| Goal | Entry point |
| --- | --- |
| Inspect the algorithms | [cps_paper_algorithms.py](cps_paper_algorithms.py), [curve_algorithms.py](curve_algorithms.py) |
| Recompute the main numerical comparisons | [Quick start](#quick-start-numerical-results-without-solver-execution) |
| Understand reproduction requirements | [Reproduction guide](docs/REPRODUCING.md) |
| Find current versus historical files | [Repository map](docs/REPOSITORY_MAP.md) |
| Generate additional compression views | [Reporting guide](reporting/README.md) |

### Quick start: numerical results without solver execution

Using Python 3.10 or newer, run from the repository root:

```sh
git clone https://github.com/galitgoz/chain-pair-simplification.git
cd chain-pair-simplification
python reporting/summarize_primary.py --output-dir reproduced/primary
```

No third-party packages are needed for this command. It checks the common cohort and writes `common_comparison.csv`, `solver_time_medians.csv`, and `summary.json`, including the source table's SHA-256. It does not run solvers or modify historical results.

The primary dataset contains **45 eligible pairs** (12 protein, 33 hurricane). The matched comparison contains **43 pairs** (12 protein, 31 hurricane) with validated outputs from all four methods. Two additional hurricane pairs have resource-limited CPS-2F observations and remain in the saved evidence. Preparation-blocked records are separate from these 45 eligible pairs.

The comparison uses `alpha = 0.5`, chosen using earlier pilot evidence. Fidelity thresholds are half the respective mean input edge lengths; the coupling threshold is the input discrete Fréchet distance. Compression differences depend on the fidelity definitions and tolerance regime. Primary observations, historical diagnostics and timing repetitions must not be pooled as distinct input pairs.

## Data and methods

Protein inputs are prepared C-alpha backbone curves from RCSB PDB. Hurricane inputs are projected NOAA HURDAT2 tracks. Source versions, preprocessing, shared protein references and interpretation limits are documented in the [reproduction guide](docs/REPRODUCING.md#data-provenance).

Output vertices are selected from the inputs in order, with fixed endpoints. Auxiliary points are internal matching locations, not extra output vertices. Feasibility validation and the algorithm's optimality claim are distinct.

## Reproduction status

The numerical summary command is checked against the saved paper counts and timing medians. The original full experiment environment used Python 3.14 on Windows; full solver execution on a clean installation has not yet been validated. See [environment requirements and limitations](docs/REPRODUCING.md#solver-execution-current-limitations) before running experiments.

Historical workflows and evidence are retained for traceability. Use the [repository map](docs/REPOSITORY_MAP.md) to distinguish them from the current comparison. The [cleanup record](docs/REPOSITORY_HYGIENE.md) documents removal of regenerable files. `snapshot_manifest.json` describes the initial import, not subsequent revisions.

## Citation and reuse

Repository: <https://github.com/galitgoz/chain-pair-simplification>. The [v1.0.0 research snapshot](https://github.com/galitgoz/chain-pair-simplification/tree/v1.0.0) provides a fixed version for citation. Citation metadata are available in [CITATION.cff](CITATION.cff), and changes are recorded in [CHANGELOG.md](CHANGELOG.md).

The original software is licensed under the [MIT License](LICENSE). This license does not relicense upstream datasets or third-party research publications; those retain their respective rights.
