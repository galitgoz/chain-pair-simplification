# Reproducing the paper comparison

## Saved-result summaries

Run from the repository root with Python 3.10+:

```sh
python reporting/summarize_primary.py
python -m unittest discover -s reporting -p test_summarize_primary.py
python -m unittest discover -s tests
```

These commands need only the standard library. The summary selects 43 pairs with validated outputs from all four methods, checks shared inputs and thresholds, and writes the comparison, timing medians and quality counts under `reproduced/primary/`. Resource exits are excluded from these medians. This reads saved evidence; it does not measure new solver times.

## Install the solver and plotting dependencies

The pinned environment uses Python 3.14. Create an isolated environment without system site-packages:

```sh
python -m venv .venv
```

Activate it with `.venv\Scripts\Activate.ps1` in Windows PowerShell, or `source .venv/bin/activate` in a POSIX shell, then run:

```sh
python -m pip install -r requirements-analysis.txt
python -m pip install -e .
python -m pip check
```

SciPy is required by Numba's linear-algebra routines. psutil supports the Windows process supervisor. The editable installation is intended for use with this repository's prepared data and supporting files.

## Run the algorithms

The 100 small-instance checks compare CPS outputs with exhaustive independent oracles:

```sh
python -m verify_cps_papers
```

They write under `reproduced/algorithm_checks/`. For the saved experiments, first inspect available configurations (no solver dependencies needed for listing):

```sh
python experiments/run.py --list
```

On Windows, rerun the hurricane example or the entire 180-configuration cohort:

```sh
python experiments/run.py --pair EP111954__EP021984 --output-dir reproduced/example
python experiments/run.py --output-dir reproduced/full
```

Filters `--domain`, `--method` and `--limit` allow smaller selections. Each output directory must be new and inside `reproduced/`. The runner uses original thresholds and resource budgets, records current source hashes, validates feasibility and compares output sizes with the reference. A nonzero exit means a startup failure or an outcome outside the expected size/resource comparison. Inspect `comparison.json` and the attempt records before interpreting an incomplete run.

The process supervisor is Windows-specific. Numerical reporting and core algorithms are separate from that supervisor. Two reference CPS-2F configurations are resource-limited; they are not evidence of infeasibility. Timing and the first resource limit reached can vary by machine.

All 180 configurations were previously rerun in an isolated Windows/Python 3.14 environment: all 178 reference solutions passed feasibility validation and reproduced both output sizes, while the same two configurations remained resource-limited. The compact package preserves the core algorithm code and input bytes. Rebuilding inputs from external sources, historical timing replication and every archived notebook are outside this prepared-input workflow.

## Figures and provenance

See the [reporting guide](../reporting/README.md) for supplementary plots and [data provenance](../data/README.md) for source versions and preparation context. Selected final auxiliary-point and coupling/size figures, with plotted values, are retained under `results/figures/`. Their original construction scripts and the full research gallery remain in the [archive](https://github.com/galitgoz/chain-pair-simplification/tree/v1.0.0).

The 12 protein pairs share reference chain `1o7j.a` and are not independent biological replicates. The fidelity parameter was chosen using earlier pilot evidence. Do not pool archived repetitions or diagnostic cohorts with the primary comparison. Feasibility validation is distinct from independently proving the algorithms' optimality claims.
