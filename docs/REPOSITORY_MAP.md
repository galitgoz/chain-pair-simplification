# Repository map

## Current comparison

| Location | Purpose |
| --- | --- |
| `cps_paper_algorithms.py`, `curve_algorithms.py` | Core algorithms |
| `reporting/summarize_primary.py` | Dependency-free numerical reproduction |
| `output/paper_v1/coverage_extension/combined_results.csv` | Combined primary records, including blocked and resource-limited outcomes |
| `output/paper_v1/coverage_extension/attempt_gallery/` | Saved attempt metrics and figures |
| `output/paper_v1/plots_auxiliary_bound_percentages/` | Auxiliary-point and IC/CPS-2F comparison figures |
| `reporting/saved/experiment_tables/` | Frozen tables for additional plotting scripts |

## Data and evidence

`data/` contains supplied/prepared curves and cached PDB files. Cohort `inputs/` directories contain frozen solver inputs. Attempt directories contain configurations, raw results, validation and resource records. Missing quantities must be interpreted alongside status fields.

## Historical files

Root-level `0.*` and `1.*` notebooks, `build_*`, `execute_*` and reporting-update scripts are earlier workflows, not an ordered current reproduction pipeline.

The `graph_only_v1`, `v2` and `v3` subdirectories document diagnostic stages. Other cohorts include `unified_pilot_results`, `expanded_comparison`, `coverage_extension`, `hurricane_extension_3` and `hurricane_final_20`. Do not concatenate all their results: that would duplicate observations and mix primary results with repetitions.

`output/paper_v1/experimental_chapter/` is an earlier chapter draft and may describe an older data snapshot. Use the combined primary records for current numerical comparisons.

Historical files retain their locations because saved evidence and notebooks may depend on relative paths.
