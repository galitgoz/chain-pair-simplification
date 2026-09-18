# Chain pair simplification

## Paper results: start here

This repository compares **CPS-2F** (continuous input fidelity, discrete output coupling), **CPS-3F** (discrete fidelity and coupling), and independent continuous/discrete baselines. The objective is `k = max(kA, kB)`.

| Goal | Entry point |
| --- | --- |
| Inspect the algorithms | [cps_paper_algorithms.py](cps_paper_algorithms.py), [curve_algorithms.py](curve_algorithms.py) |
| Recompute the main numerical comparisons | Standard-library command below |
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

## Historical snapshot and implementation notes

Research code, experiment notebooks, and saved experimental records for continuous and discrete chain pair simplification.

### File inventory

- `cps_paper_algorithms.py` and `curve_algorithms.py`: core algorithms.
- `notebooks_paper/`: experiment orchestration, correctness checks, reporting, and notebooks. Subdirectories retain the historical experiment stages.
- `data/`: supplied hurricane tracks and protein inputs/prepared data.
- `output/paper_v1/`: saved experimental records, configurations, validations, tables, and figures. These are historical observations, not new executions.
- `requirements-analysis.txt`: pinned analysis environment requirements.
- Existing `README-*.md` and stage READMEs document earlier workflows and may describe an earlier stage rather than the latest analysis.

### Full analysis environment

Create a separate Python environment and install the dependencies from `requirements-analysis.txt`. The saved experiment environment used Python 3.14. Additional runner dependencies may be required; consult imports and the saved environment audits. This snapshot has not yet been tested by rerunning the full experiment suite on a clean machine.

Open this directory as the project root in PyCharm. Start by reading saved results and notebook code. Do not execute all notebooks indiscriminately: historical runners can write to fixed output paths. Use a fresh output directory for new experiments. Notebook embedded outputs were cleared for version control; external result records and figures are retained.

Some historical helper scripts and provenance records contain original absolute paths. Those records describe the source machine; adapt operational paths before rerunning an affected script.

### Interpretation

IC is the independent continuous baseline. CPS-2F uses continuous fidelity to each input and discrete coupling between simplified outputs; CPS-3F uses discrete distances for all three constraints. Auxiliary points are internal matching locations, not extra output vertices. The output objective is `k = max(kA, kB)`.

The primary common comparison contains 43 pairs (12 protein and 31 hurricane pairs) with validated outputs from all methods. Historical incomplete attempts are retained for provenance and are not part of that completed comparison. Timing repetitions are separate observations, not additional distinct pairs.

### Snapshot and data provenance

`snapshot_manifest.json` records the initial copied file hashes, not a live integrity manifest for later repository revisions. The original working project was left unchanged. See the [data provenance guide](docs/REPRODUCING.md#data-provenance), data READMEs and per-run manifests for source provenance. Source data and third-party materials retain their respective rights. No redistribution license is asserted by this snapshot.

### Citing and reusing the code

Repository: <https://github.com/galitgoz/chain-pair-simplification>. Record the commit used alongside this URL when citing experiments. No numbered release is implied. A software license has not yet been selected; public visibility alone does not grant a reuse license.
