# Graph-only Stage 1 revision

This revision supersedes the certificate-enabled pilot for thesis graph evaluation. Previous files remain historical and unchanged. Read `output/paper_v1/graph_only_v1/status.md` and the notebooks here for the new results.

All actual CPS calls use `graph_core.solve(..., certificate=False)`, which in turn explicitly calls the unchanged `solve_cps(..., certificate=False)`. Passing True is rejected. The two independent methods remain separate baselines. No alternative-subsequence method is evaluated.

Execution order for this recorded run:

1. `setup_revision.py` freezes the same inputs/parameters and inventories prior results. It refuses reinitialization.
2. `graph_checks.py`, then `shortcut_probe.py`, verify small inputs and both disabled shortcuts. The final gate has 100 checks.
3. `run_pilot.py protein hurricane` runs exactly the frozen 32 configurations. It refuses to overwrite result tables.
4. `report_revision.py`, `build_notebooks.py`, and `render_notebooks.py` read saved results and render the Stage 1 report. They do not optimize.

Reopening the notebooks reads saved results, not new timing repetitions. A future benchmark must use this explicit-False adapter and a newly versioned output directory after authorization. Graph observation records entry into the product recurrence; baseline single-curve graphs are excluded. Resource exits remain separate from infeasibility. No later stage is launched.
