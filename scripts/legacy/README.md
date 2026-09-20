# Historical maintenance scripts

These 12 helpers were used to inspect sources, modify notebook builders, add reporting cells, or run specific earlier experiment settings. They are not the current paper reproduction entry point.

Use [the main README](../../README.md) for the saved-result summary workflow. Preserve the historical output locations if investigating an earlier workflow; some referenced legacy outputs and local kernels are not included in this snapshot.

If intentionally running a helper, run it from the repository root, for example:

```sh
python scripts/legacy/inspect_sources.py
```

Scripts with `add_`, `update_`, `refresh_`, `execute_` or `expand_` in their names may modify notebooks, builders or output files. They were not rerun during repository reorganization. Root-path resolution and the notebook builder's import of `add_graph_reporting` were updated for this location.

Files: `add_dual_bar_reporting.py`, `add_graph_reporting.py`, `add_parameter_reporting.py`, `update_cps_notebook.py`, `refresh_cps_reporting.py`, `refresh_hurricane_tracks.py`, `expand_hurricane_analysis.py`, `finalize_expanded_hurricanes.py`, `execute_cps_w8.py`, `render_paper_checks.py`, `inspect_sources.py`, and `inspect_cps_sources.py`.
