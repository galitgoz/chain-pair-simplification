# Repository map

## Current comparison

| Location | Purpose |
| --- | --- |
| `src/cps_paper_algorithms.py`, `src/curve_algorithms.py` | Core algorithms |
| `reporting/summarize_primary.py` | Dependency-free numerical reproduction |
| `output/paper_v1/coverage_extension/combined_results.csv` | Combined primary records, including blocked and resource-limited outcomes |
| `output/paper_v1/coverage_extension/attempt_gallery/` | Saved attempt metrics and figures |
| `output/paper_v1/plots_auxiliary_bound_percentages/` | Auxiliary-point and IC/CPS-2F comparison figures |
| `reporting/saved/experiment_tables/` | Frozen tables for additional plotting scripts |

## Data and evidence

`data/` contains supplied/prepared curves and cached PDB files. Cohort `inputs/` directories contain frozen solver inputs. Attempt directories contain configurations, raw results, validation and resource records. Missing quantities must be interpreted alongside status fields.

## Historical files

Earlier `0.*` and `1.*` notebooks live in `notebooks/legacy/`. Their `build_*` and `execute_*` scripts live with the analysis modules in `src/`; these are earlier workflows, not an ordered current reproduction pipeline. Twelve one-off maintenance helpers live in [scripts/legacy](../scripts/legacy/README.md). Historical guides are under `docs/legacy/`, and the initial import manifest is under `docs/provenance/`.

The `graph_only_v1`, `v2` and `v3` subdirectories document diagnostic stages. Other cohorts include `unified_pilot_results`, `expanded_comparison`, `coverage_extension`, `hurricane_extension_3` and `hurricane_final_20`. Do not concatenate all their results: that would duplicate observations and mix primary results with repetitions.

`output/paper_v1/experimental_chapter/` is an earlier chapter draft and may describe an older data snapshot. Use the combined primary records for current numerical comparisons.

The paper workflows in `notebooks_paper/` and all input data and saved experimental evidence retain their locations. Install the source modules with `python -m pip install -e .` before running solver or notebook workflows; run script commands from the repository root.
