# Graph-only v2: orchestration repair and three diagnostics

Start with `00_repair_and_diagnostics.ipynb`. It reads saved results only. The full report and three-row summary are in `output/paper_v1/graph_only_v2`.

The optimization code remains in the unchanged root modules. `common.py` supplies observation, atomic persistence helpers, warm-up and external feasibility validation. `worker.py` separates phases and checkpoints raw output before external validation. `supervisor.py` monitors only its own uniquely identified worker tree with separate phase budgets.

`prepare.py` reconciled v1 and froze inputs, limits and inherited provenance. `regressions.py` ran 21 small checks, including injected validation and optional-measurement failures. `references.py` validated existing reference witnesses without optimization. `run_diagnostics.py` is guarded against restarting; it executed exactly A, B and C once. Do not rerun it in this version.

Each diagnostic directory contains its committed configuration, phase and resource logs, and either a raw solver checkpoint or a solver-exit record. Validation records are separate. Optional measurements cannot erase a raw result or convert a successful solve into a solver failure. The existing CPS `finish()` validation still happens inside solver time and is documented as such.

No new CPS solution was obtained in these three diagnostics: A and B hit transition limits. C is a newly validated independent-discrete solution. The 22 inherited verified rows are not fresh executions or timing repetitions. All CPS calls retain explicit `certificate=False`; no shortcut, optimization or preprocessing algorithm was changed.

This version is complete. Further attempts need a new request and a new versioned directory.
