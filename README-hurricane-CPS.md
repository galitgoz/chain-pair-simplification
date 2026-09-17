# Hurricane CPS analysis

Open `1.analyze_CPS.ipynb` from the project root and run its cells in order with
the project Python environment. It reuses `curve_algorithms.py` and
`cps_paper_algorithms.py`; hurricane loading, experiments and plots live in
`hurricane_cps.py`.

The expanded notebooks load all supplied tracks, select 16 pairs per basin across the
supplied ball-packedness rankings, and runs Independent GCS, CPS-2F and CPS-3F
at alpha 1 and 2 (32 pairs, 192 method results). All selected input vertices are retained before simplification.
Set `PAIRS_PER_BASIN=None` to analyze every candidate pair. The time limit is
soft, now 120 seconds, with a 20-million-transition budget; preprocessing and certificate searches can exceed it. Missing joint
solutions remain explicitly marked as resource limits.

Coordinates and all thresholds are in projected kilometres. Each basin uses a
different supplied coordinate frame; pairs never cross basins. Original-pair
distances are recomputed because some supplied distance metadata violates the
continuous/discrete ordering. There are no timestamps in these coordinate files.

Expanded outputs are in `output/hurricane_cps_expanded/` and
`output/hurricane_auxiliary_expanded/`; earlier eight-pair results remain in their
original folders. Outputs include comparison and parameter tables, distance
audit, simplification indices and coupling paths, graph counts, auxiliary
utilization, manifest and PDF/SVG/PNG figures. Protein notebooks and their
results are separate.

Automation:

```powershell
.venv/Scripts/python.exe execute_hurricane_cps.py
.venv/Scripts/python.exe validate_hurricane_cps.py output/hurricane_cps_expanded
.venv/Scripts/python.exe validate_hurricane_auxiliary.py output/hurricane_auxiliary_expanded
```

`build_hurricane_cps_notebook.py` recreates the notebook template and overwrites
the notebook; use it only when intentionally discarding existing notebook edits
and outputs.
