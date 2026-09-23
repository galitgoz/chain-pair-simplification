# Reference results

- `primary.csv` is the original combined table, preserved byte-for-byte. It contains 180 observed method configurations on 45 eligible pairs, plus eight preparation-blocked method records. The matched comparison uses 172 observations on 43 pairs.
- `evidence.json` consolidates the available original JSON records for each of the 180 primary attempts: configuration, raw solver checkpoint, validation, process summary, errors, progress and optional measurements. File hashes refer to original bytes; payloads have been reserialized into this container.
- `plot_data.json` contains the 180 primary rows from the frozen plotting data; timing repetitions remain in the archive.
- `figures/` retains the selected auxiliary-point and coupling/size figures with their plotted CSV values. Supplementary compression figures can be regenerated using [reporting scripts](../reporting/README.md).

Paths embedded in reference tables and evidence describe the original experiment layout. Resolve them against the [v1.0.0 archive](https://github.com/galitgoz/chain-pair-simplification/tree/v1.0.0), not this compact checkout. Current prepared-input paths are in `data/inputs.json`; current runnable configurations are in `experiments/configurations.json`.

New calculations belong under `reproduced/` and must not replace these reference records. Historical times exclude warm-up and external validation, but include checks performed inside the solver. Resource exits are not proofs of infeasibility.
