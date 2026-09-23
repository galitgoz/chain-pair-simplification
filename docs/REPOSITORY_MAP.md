# Repository map

The default branch contains the paper reproduction package. Start with the [main README](../README.md).

- `src/curve_algorithms.py` and `src/cps_paper_algorithms.py`: core algorithms.
- `data/inputs.json`: manifest of all 45 prepared pairs, their original paths, vertex counts and SHA-256 hashes.
- `experiments/configurations.json`: 180 configurations, with current input paths, original numerical thresholds and budgets, and reference output sizes.
- `results/primary.csv`: unmodified primary table, including resource-limited and preparation-blocked rows.
- `results/evidence.json`: original configuration, raw-result, validation and process-summary JSON payloads for the 180 primary attempts. Historical path fields remain unchanged.
- `results/plot_data.json`: the 180 primary observations from the frozen plotting table; timing repetitions are excluded.
- `results/figures/`: selected final figures and their plotted values.
- `reporting/`: saved-result summaries and supplementary figure generation.

Historical experiments, raw-source preparation, notebooks, per-process logs and earlier figure versions are available in the [v1.0.0 archive](https://github.com/galitgoz/chain-pair-simplification/tree/v1.0.0). They are not required for the documented prepared-input workflow.
