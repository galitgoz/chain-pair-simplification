# Reproducing results

## Numerical summaries

Run from the repository root with Python 3.10+:

```sh
python reporting/summarize_primary.py --output-dir reproduced/primary
```

This uses only the standard library and reads `output/paper_v1/coverage_extension/combined_results.csv`. It excludes preparation-blocked and nonvalidated records, rejects duplicate validated methods, and checks the common cohort: 12 protein and 31 hurricane pairs, each with all four methods.

Outputs contain the matched observations, CPS solver-time medians, IC/CPS-2F objective and coupling comparisons, and CPS-2F/CPS-3F objective comparisons. Running again replaces generated summaries only. Destinations inside the source `data/` or `output/` directories are rejected.

These are saved-result reproductions, not new optimizer executions or timing measurements.

Run the summary and repository-layout checks with:

```sh
python -m unittest discover -s reporting -p test_summarize_primary.py
python -m unittest discover -s tests
```

The summary checks verify recorded paper counts and rounded medians and reject duplicate observations or mismatched thresholds. The layout checks verify import locations and resolution of relocated provenance paths. Neither suite executes the optimizers.

## Additional figures

Follow [reporting/README.md](../reporting/README.md) to create the two supplementary compression views. These scripts require Matplotlib and write under `reporting/saved/compression_preview/`. This is not a single-command reproduction of every manuscript figure.

The analysis package versions are pinned in `requirements-analysis.txt`, including SciPy (required by Numba's linear-algebra routines) and psutil (used by the process supervisor). Create a separate environment with Python 3.14:

```sh
python -m venv .venv
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
python -m pip install -r requirements-analysis.txt
python -m pip install -e .
```

For saved-result analysis on a POSIX shell:

```sh
source .venv/bin/activate
python -m pip install -r requirements-analysis.txt
python -m pip install -e .
```

The pinned dependencies and editable installation were verified in an isolated Windows/Python 3.14.0 environment with system and user site-packages disabled. `python -m pip check`, the six summary/layout tests, 100 small-instance algorithm checks, the documented plotting commands and a fresh notebook kernel passed. The dependency-free summary command is a separate entry point and does not require this environment.

All 180 observed configurations were also rerun from the saved inputs, using their original thresholds and resource budgets with fresh configurations recording the current source hashes. All 178 previously validated solutions passed feasibility validation and reproduced both output vertex counts. The same two CPS-2F configurations remained resource-limited; one reached its time limit rather than its historical transition limit. These checks do not reproduce historical timings or rebuild inputs from external sources. New run artifacts were kept outside the published repository.

## Solver execution: current limitations

Python modules now live in `src/`. The editable installation above preserves their original import names and requires this repository's data and supporting files. Historical scripts generally expect the repository root as their working directory and may overwrite archived outputs; inspect their destinations before execution. Earlier notebooks are in `notebooks/legacy/`; their initial setup cell locates the repository and sets the working directory. Install the modules in the same environment as the notebook kernel.

Historical provenance paths are mapped through `docs/provenance/relocated_files.json` when checking saved source hashes. Recorded results and hashes are not rewritten. Strict checks of historical orchestration code may still reject a resume after code changes; use the cited `v1.0.0` snapshot when investigating that exact historical layout.

- Use a fresh output location. Historical scripts can write fixed paths, refuse reruns or resume old schedules; do not execute all notebooks as one pipeline.
- `notebooks_paper/graph_only_v3/supervisor.py` uses Windows-specific process flags. That solver runner is not currently portable to POSIX.
- Protein preparation may require US-align. `src/setup_structural_alignment.py` downloads its Windows executable from the upstream provider. It is unnecessary for reading saved inputs or numerical summaries.
- Preparation in `notebooks_paper/coverage_extension.py` also expects `papers/chain_pair_simplification.pdf`, which is not bundled. Supply that research source or refactor the preparation dependency before rebuilding the inputs.
- Original machine paths remain in some historical helpers and provenance records. Distinguish operational paths from archived provenance text.

Small-instance correctness checks are in `src/verify_cps_papers.py` and `notebooks_paper/correctness_checks.py`. Their runners write result files and should be isolated before execution. Checking the saved-result summaries is not a substitute for checking algorithm correctness.

For example, `src/verify_cps_papers.py` writes `output/cps_papers/algorithm_checks.json` relative to the working directory. After the editable installation, its small-instance checks can be run in an empty temporary working directory with `python -m verify_cps_papers`, keeping that generated report separate from archived evidence. This does not apply automatically to other runners, which may derive fixed output paths from their source location.

## Data provenance

Proteins originate from [RCSB PDB](https://www.rcsb.org/). Download records are in `data/pdb_cache/sources.json`; chain/residue mappings and omissions are in `output/paper_v1/coverage_extension/protein_audit/`. Inputs use selected C-alpha positions, alignment and endpoint-preserving sampling. The 12 pairs share the reference chain `1o7j.a` and are not independent biological replicates.

Hurricanes originate from [NOAA HURDAT2](https://www.nhc.noaa.gov/data/), specifically:

- [Atlantic 1851–2015, July 6, 2016 release](https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-1851-2015-070616.txt)
- [Northeast and North Central Pacific 1949–2015, May 9, 2016 release](https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-nepac-1949-2015-050916.txt)

`data/hurricanes/README.md` documents filtering, consecutive-duplicate removal and basin-specific projection. Experiments retain the supplied track vertices; those inputs are processed planar coordinates, not unprocessed latitude/longitude records.

## Interpretation

The common-cohort medians describe one validated primary observation per pair and method. Resource exits are excluded. Solver time excludes warm-up and external validation, but includes validation performed inside the solver. Feasibility validation is distinct from independent evidence of optimality. Resource exits do not prove infeasibility.
