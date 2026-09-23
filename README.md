# Chain pair simplification

Research code and saved primary results for **CPS-2F** (continuous input fidelity, discrete output coupling), **CPS-3F** (discrete fidelity and coupling), and independent baselines. The objective is `k = max(kA, kB)`; output vertices preserve input order and fixed endpoints.

## Reproduce the numerical results

Python 3.10+; no third-party packages required:

```sh
git clone https://github.com/galitgoz/chain-pair-simplification.git
cd chain-pair-simplification
python reporting/summarize_primary.py
```

Generated files go under `reproduced/`. The matched comparison contains **43 pairs** (12 proteins, 31 hurricanes). The retained inputs include all **45 eligible pairs**; two hurricane CPS-2F configurations are resource-limited. Preparation-blocked records remain in the reference table.

The experiment uses `alpha = 0.5`: fidelity thresholds are half each input's mean edge length, and the coupling threshold is the input discrete Fréchet distance. Compression differences depend on this tolerance regime as well as the fidelity definitions.

## Repository contents

| Directory | Contents |
| --- | --- |
| [src](src/) | Two core algorithm modules and small-instance correctness checks |
| [data](data/README.md) | 45 prepared input pairs, hashes and source information |
| [experiments](experiments/run.py) | Primary configurations and isolated solver runner |
| [results](results/README.md) | Original primary table, consolidated evidence and selected figures |
| [reporting](reporting/README.md) | Numerical summaries and plotting scripts |
| [tests](tests/) | Input and configuration consistency checks |
| [docs](docs/REPRODUCING.md) | Installation and reproduction instructions |

For solver execution, install the dependencies and local modules as described in the [reproduction guide](docs/REPRODUCING.md). The process supervisor currently requires Windows. New runs never overwrite the reference results.

## Reproduction status

An isolated Windows/Python 3.14 environment was used to rerun all 180 observed configurations from saved inputs: 178 returned validated solutions with the original output sizes; the same two configurations remained resource-limited. This verifies execution from prepared inputs, not reconstruction from raw external sources or reproduction of historical timings.

## Archive and citation

The [v1.0.0 research archive](https://github.com/galitgoz/chain-pair-simplification/tree/v1.0.0) retains historical notebooks, intermediate experiments, source preparation and detailed logs. The `main` branch is the compact reproduction package. Historical path fields in saved results refer to that archive.

[CITATION.cff](CITATION.cff) describes `v1.0.0`; record the exact commit when using newer code. See [CHANGELOG.md](CHANGELOG.md) for changes. Original software uses the [MIT License](LICENSE); upstream datasets retain their own rights.
