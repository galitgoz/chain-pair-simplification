# Protein chain-pair analysis

Open `0.analyze_CPS.ipynb` from this project directory and run its cells in order.
Use the project Python environment and dependencies in `requirements-analysis.txt`.

The notebook contains executable implementations of global continuous GCS
(van de Kerkhof et al., Algorithm 1, represented by integer cost layers),
Fan et al.'s CPS-3F configuration-graph algorithm, and the continuous-fidelity,
discrete-pairing CPS-2F graph and dynamic program from Galit Gozoltzani's thesis.
All three preserve the input endpoints. The joint objective is
`min max(number of A vertices, number of B vertices)`.

## Inputs and defaults

- `WS = (8,)`: the current requested sampling level; change to any subset of
  `(1, 2, 4, 8, 16)`. Finer levels can require substantially more graph resources.
- `ALPHAS = (1, 2)`: fidelity tolerances are alpha times mean input edge length.
- `INPUT_SCOPE = "supplied"`: retains every supplied aligned CSV vertex, including
  edges across gaps. It does not crop or decimate the supplied curves again.
- The seven Table 3 pairs are supplemented from cached full audited chains,
  sampled independently with endpoints retained. Their rigid transforms come
  from the previous US-align preparation. Unequal chain lengths are preserved.
- `INPUT_SCOPE = "qualified"` selects the previous common continuous intervals.

The pairing bound is the ceiling of the actual input discrete Fréchet distance.
This is an experiment on the available coordinates, not a reproduction of the
historical thresholds or result counts in Tables 1–3. The provided thesis is an
annotated draft; the notebook identifies the dependency on its discretization claim.

## Outputs

`output/cps_papers_w8/` contains the current measured CSV and HTML tables, selected vertex
indices and coupling paths, input/residue provenance, validation report, manifest,
and PNG figures. The executed notebook embeds tables and figures. Earlier stride-16
results remain in `output/cps_papers/`.

The analysis reports `|A|`, `|B|`, `|A′|`, `|B′|`, compression factors and percentage
removed. `auxiliary_counts.csv` lists CPS-2F's added points on each curve and the
extended-curve sizes, including diagnostic counts for certificate solutions.
`cps2f_vs_cps3f.csv` compares the methods on matched pairs and all three tolerances.
Every comparison plot identifies the pair and delta1/delta2/delta3.

Structural fidelity includes continuous Fréchet error and nonlocal contact-map
precision, recall and F1. Contacts use C-alpha distances at most 8 Å and polymer
sequence separation at least 5. Discarded input vertices are mapped to shortcuts
by original subcurve arc fraction, solely for evaluation on a common grid.
Contact preservation is a structural proxy, not a knot or entanglement guarantee.
The full CSV includes true/false-positive and false-negative contact counts.

`refresh_cps_reporting.py` checks solver/input hashes and refreshes measurements
and plots from saved solver outputs without repeating expensive optimization.

Optimal independent-bound certificates and actual configuration-graph execution
have separate status labels. Resource limits remain explicit missing results.
The small-instance validation forces graph execution and compares against
exhaustive anchored subsequences.

The notebook's data and plotting helpers live in `cps_notebook_data.py`,
`cps_notebook_run.py`, and `cps_notebook_plots.py`. The builder embeds algorithm
source from `curve_algorithms.py` and `cps_paper_algorithms.py`.
To regenerate and execute the notebook, run:

```powershell
.venv/Scripts/python.exe execute_cps_w8.py
```

This runner preserves the current notebook cells, selects w=8 and alpha=1,2, and
executes from a fresh kernel. The older `execute_cps_notebook.py` instead regenerates
from `build_cps_notebook.py`, replacing manual notebook edits with the template.
