# Summaries and supplementary plots

From the repository root, run `python reporting/summarize_primary.py` to reproduce the primary numerical comparison without solver execution or third-party dependencies.

After installing the analysis dependencies, run:

```sh
python reporting/plot_retention.py
python reporting/plot_compression_scatter.py
```

The first command reads `results/plot_data.json` and generates exact plotted values and a retention chart under `reproduced/figures/`. The scatter script reads those generated values, so run the commands in order. `python reporting/plot_tradeoff.py` produces an exploratory coupling/compression view under `reproduced/tradeoff/`. Generated previews are not versioned.

Compression is `100 * (1 - (kA+kB)/(nA+nB))`, distinct from the objective `max(kA,kB)`. These views use 43 pairs with validated saved outputs from all four methods. Historical machine paths in the frozen table are provenance fields and are not used to load plot inputs. Selected final figures are available under `results/figures/`; the full historical gallery remains in the research archive.
