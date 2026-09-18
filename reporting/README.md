# Saved-result summaries and additional figures

For the current paper's numerical comparisons, run `python reporting/summarize_primary.py --output-dir reproduced/primary` from the repository root. This requires only Python 3.10+ and writes the common cohort, solver-time medians, and quality counts without changing saved evidence. See the [reproduction guide](../docs/REPRODUCING.md).

## Additional figure scripts

Run `python reporting/plot_retention.py`, followed by `python reporting/plot_compression_scatter.py` from the repository root. They read the supplied frozen table and write figures under `reporting/saved/compression_preview`. No solvers are invoked. `plot_tradeoff.py` is an exploratory alternative, not a selected paper figure.

Compression is 100 * (1 - (kA+kB)/(nA+nB)); this is distinct from the optimization objective max(kA,kB). The figure set includes only pairs with validated saved outputs from all four methods. Metadata may retain historical source-machine paths; these scripts do not use those paths to load results.
