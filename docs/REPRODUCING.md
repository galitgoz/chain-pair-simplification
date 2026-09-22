# Reproducing results

## Numerical summaries

Run from the repository root with Python 3.10+:

```sh
python reporting/summarize_primary.py --output-dir reproduced/primary
```

This uses only the standard library and reads `output/paper_v1/coverage_extension/combined_results.csv`. It excludes preparation-blocked and nonvalidated records, rejects duplicate validated methods, and checks the common cohort: 12 protein and 31 hurricane pairs, each with all four methods.

Outputs contain the matched observations, CPS solver-time medians, IC/CPS-2F objective and coupling comparisons, and CPS-2F/CPS-3F objective comparisons. Running again replaces generated summaries only. Destinations inside the source `data/` or `output/` directories are rejected.

These are saved-result reproductions, not new optimizer executions or timing measurements.

Run the summary regression checks with:

```sh
python -m unittest discover -s reporting -p test_summarize_primary.py
python -m unittest discover -s tests
```

They verify the recorded paper counts and rounded medians and reject duplicate observations or mismatched thresholds. They do not execute the optimizers.

## Additional figures

Follow [reporting/README.md](../reporting/README.md) to create the two supplementary compression views. These scripts require Matplotlib and write under `reporting/saved/compression_preview/`. This is not a single-command reproduction of every manuscript figure.

The historical analysis package versions are in `requirements-analysis.txt`. Create a separate environment with the original Python 3.14 interpreter:

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

A clean installation of this historical environment and a full solver rerun have not yet been validated. The dependency-free summary command is a separate entry point and does not require that environment.

## Solver execution: current limitations

Python modules now live in `src/`. The editable installation above preserves their original import names. Run script commands from the repository root, using the new path (for example, `python src/verify_cps_papers.py`). Earlier notebooks are in `notebooks/legacy/`; their initial setup cell locates the repository and sets the working directory. Install the modules in the same environment as the notebook kernel.

Historical provenance paths are mapped through `docs/provenance/relocated_files.json` when checking saved source hashes. Recorded results and hashes are not rewritten. Strict checks of historical orchestration code may still reject a resume after code changes; use the cited `v1.0.0` snapshot when investigating that exact historical layout.

- Use a fresh output location. Historical scripts can write fixed paths, refuse reruns or resume old schedules; do not execute all notebooks as one pipeline.
- `notebooks_paper/graph_only_v3/supervisor.py` uses Windows-specific process flags. That solver runner is not currently portable to POSIX.
- The supervisor imports `psutil`; a solver environment should include it explicitly rather than relying on transitive notebook dependencies.
- Protein preparation may require US-align. `src/setup_structural_alignment.py` downloads its Windows executable from the upstream provider. It is unnecessary for reading saved inputs or numerical summaries.
- Preparation in `notebooks_paper/coverage_extension.py` also expects `papers/chain_pair_simplification.pdf`, which is not bundled. Supply that research source or refactor the preparation dependency before rebuilding the inputs.
- Original machine paths remain in some historical helpers and provenance records. Distinguish operational paths from archived provenance text.

Small-instance correctness checks are in `src/verify_cps_papers.py` and `notebooks_paper/correctness_checks.py`. Their runners write result files and should be isolated before execution. Checking the saved-result summaries is not a substitute for checking algorithm correctness.

## Data provenance

Proteins originate from [RCSB PDB](https://www.rcsb.org/). Download records are in `data/pdb_cache/sources.json`; chain/residue mappings and omissions are in `output/paper_v1/coverage_extension/protein_audit/`. Inputs use selected C-alpha positions, alignment and endpoint-preserving sampling. The 12 pairs share the reference chain `1o7j.a` and are not independent biological replicates.

Hurricanes originate from [NOAA HURDAT2](https://www.nhc.noaa.gov/data/), specifically:

- [Atlantic 1851–2015, July 6, 2016 release](https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-1851-2015-070616.txt)
- [Northeast and North Central Pacific 1949–2015, May 9, 2016 release](https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-nepac-1949-2015-050916.txt)

`data/hurricanes/README.md` documents filtering, consecutive-duplicate removal and basin-specific projection. Experiments retain the supplied track vertices; those inputs are processed planar coordinates, not unprocessed latitude/longitude records.

## Interpretation

The common-cohort medians describe one validated primary observation per pair and method. Resource exits are excluded. Solver time excludes warm-up and external validation, but includes validation performed inside the solver. Feasibility validation is distinct from independent evidence of optimality. Resource exits do not prove infeasibility.
