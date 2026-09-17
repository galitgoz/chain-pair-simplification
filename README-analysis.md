# Protein Curve Simplification

The current stride-decimation experiment is **1.analyze_independent.ipynb**. It uses the exact requested shared-index decimator, w=1/2/4/8/16, bond-scaled alpha=1/2/4 thresholds, and a local-continuous joint shortcut solver. Its CSV, rendered table, F1-F4 figures, residue mappings and qualification audit are under `output/decimation/`. Read its prominent qualification note: unequal or broken chains use explicitly reported common continuous row intervals. Both methods preserve the endpoints of those intervals at every stride. w=8/16 are labelled runtime-only; summaries include the w=1 control.

To regenerate this version, run `execute_decimation_notebook.py`. The older `execute_small_notebook.py` builds the superseded 30-residue experiment and should not be used to rebuild the current notebook.

Open **0.analyze_independent.ipynb** from this directory and use the project's Python environment. The notebook includes executed outputs and defines the continuous global solver in executable cells.

```powershell
.venv/Scripts/python.exe -m pip install -r requirements-analysis.txt
```

The existing `data/raw protein backbones` files are unchanged. Cached primary PDB/mmCIF records augment them for residue identities and missing-region verification. If the cache is absent, run `prepare_sources.py` (network access required).

The notebook uses a verified 69-run cached sweep by default, runs a full-chain example live, and rechecks saved fidelity. Set `RECOMPUTE_ALL = True` to recompute every fragment. Do this after editing inline solver cells. Cache signatures detect changes to shipped modules, structural records, raw inputs, parameters, and preparation policies; they cannot detect arbitrary edits to inline notebook functions.

Outputs:

- `data/prepared/`: separate fragment arrays in angstroms, chain-to-residue CSV mappings, and preparation policy.
- `output/identifier_audit.csv`: printed identifiers, verified interpretations, reported/deposited/prepared lengths, chain namespaces, and provenance.
- `output/chain_breaks.csv` and `output/missing_residues.csv`: explicit exclusions and connectivity findings.
- `output/independent/`: counts, continuous distance brackets, discrete distances, timings, selected indices, and retained-residue mappings.
- `output/paper_experiments.csv` and `output/pair_jobs.json`: all 21 paper experiments and pending pairing/registration choices.
- `output/figures/`: common-view 3D fragments, simplification overlays, and summary figures.

`curve_algorithms.py` supplies the independently testable implementation, also embedded in the notebook. `protein_data.py` handles traceable preparation. `verify_algorithms.py` checks small-instance global optimality, degenerate geometry, and metric/endpoint differences. `benchmark_fragments.py` runs the 69-case independent sweep. `build_notebook.py` generates the notebook source; `execute_notebook.py` validates the stored sweep and executes it with the current interpreter.

The independent solver follows Algorithm 1 in the supplied **full** GCS paper through an equivalent integer-cost-layer representation of its global free-space recurrence, with O(n^4) worst-case time. It does not substitute local shortcut validation for global matching. The joint discrete CPS implementation is an exact, size-guarded exhaustive reference for small windows, not the paper's scalable Algorithm 2. Full joint protein runs remain pending. CPS-2F is an explicit unimplemented interface.

The notebook documents unresolved historical subsets, likely identifier typos, registration, fragment pairing, and endpoint-policy differences. Current PDB counts are not forced to match paper counts, and missing regions are never connected.
