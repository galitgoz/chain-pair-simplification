# Single controlled graph-only v3 experiment

Start with `00_controlled_transition_budget.ipynb`; it reads saved output and never runs optimization. Results are under `output/paper_v1/graph_only_v3`.

The one authorized configuration is CPS-2F on AL031854__AL011991 at alpha=0.5. Only max_transitions changes from 2 million to 10 million. All inputs, thresholds and other limits match v2. Both 25-vertex tracks are reused unchanged. `run_once.py` refuses a second launch.

The v2 phase worker, atomic checkpointing and worker-tree supervision are copied into this separate version. The original optimization modules are untouched. All CPS calls explicitly disable certificates. The independent lower-bound calculation remains part of the solver's existing preprocessing, but cannot bypass the graph. The reference k=5 is stored separately and is read only by post-run reporting.

`instrumentation.md` defines counter semantics and overhead. Raw outputs, external validation and optional measurements are separate. The solver's original internal validation remains inside solver time. Logs and final counters are diagnostic observations, not runtime repetitions or exact algorithm-memory measurements.

This experiment is completed and must not be rerun. Further work requires a separate request and version.
