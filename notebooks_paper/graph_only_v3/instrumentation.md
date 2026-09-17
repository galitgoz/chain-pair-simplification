# Instrumentation for the single v3 diagnostic

The optimization modules are unchanged. The v2 worker, atomic persistence and phase-budget supervisor are reused in a separate directory. The worker adds an explicit `certificate_enabled=False` field to the raw checkpoint. All CPS calls still pass False, and the reference objective is not supplied to the solver.

The existing chain-graph wrapper records completion of each direct single-curve graph construction in `solve_cps`; independent lower-bound preprocessing is not counted as joint construction. At the existing heap-pop wrapper, v3 records the number of dequeued states, queue length, discovered-state count and transition counter from the calling solver frame. It returns the unmodified result of the original heap pop. No transition, ordering, pruning rule or termination condition changes.

Progress is emitted at first product-graph entry and then at most once per second, only at heap-pop boundaries. The clock and small counter updates occur once per popped state, not per transition. Log flushes and frame inspection add unmeasured diagnostic overhead; times are not runtime benchmark measurements. Single-curve completion events are also retained. Parent RSS sampling and phase timing are unchanged from v2.

`states_dequeued` counts states removed from the queue, including a terminal state if reached. `completed_states_lower_bound` counts earlier dequeued states whose processing finished before a subsequent pop. It is conservative on a timeout checked at the end of a state's processing. It is not a percentage of completion. Queue size is not remaining total work: undiscovered states and their transitions are unknown.

On a Python resource exception, final transition/state/queue counters are read from the existing traceback frame, and the heap-pop counts from the adapter frame. On a hard worker termination, only the last logged observations are available and are not relabelled as exact final counters. On a successful return, the solver's returned transition/state counters are authoritative.

The original CPS `finish()` still validates fidelity and coupling before returning; that remains in solver time. The raw result is then atomically saved before additional external validation and optional distance measurements. No per-phase algorithm benchmark or exact algorithm-memory claim is made.
