# Stage 1 — completed execution, partial pilot coverage

## Completed
- New notebooks 00–03 and shared adapters under notebooks_paper. Existing algorithm implementations reused unchanged.
- Frozen inputs: two protein pairs at existing w=16 (22 vertices per curve), two hurricane pairs (25 per curve); no further sampling or alignment.
- 8 pair/alpha configurations, 4 methods each, 32 attempted configurations. alpha=0.5,1; separate mean-edge fidelity thresholds; exact unrounded input discrete coupling budget.
- Correctness: 99/99 checks PASS; 0 failures; 0 incomplete. Exhaustive small-input checks force certificates off. Finite checks do not prove the draft theorem or eliminate shared distance-oracle errors.
- Pilot optimal outputs: 14/32; missing: 18/32. See [pilot summary](pilot_summary.csv), [full results](pilot_results.csv), [parameters](pilot_parameters.csv), and [validation](validation.csv).
- Preservation audit: no changes to the 219 inventoried old files or 1,054 old output artifacts. Existing Jupyter kernel was not attached to or stopped.

## Failures and limitations
- Protein attempts used a launcher-only Windows memory monitor. Memory values are invalid and excluded; no solution was returned. Original attempts remain in logs/protein_initial_monitoring_raw.csv. They were not repeated. Two warm-ups (one protein, one hurricane) hit the soft resource cap before the pilot solver started and are classified startup_resource_limit after reviewing their tracebacks; these are not algorithmic correctness failures.
- Corrected hurricane monitoring sums RSS over the launcher and descendants and terminates the worker tree on a cap. Polling limits can overshoot while the OS delays scheduling. The pre-existing kernel is excluded.
- {'time_limit': 15, 'optimal': 14, 'startup_resource_limit': 2, 'other_resource_limit': 1}
- Resource exits do not establish mathematical infeasibility. No failed result is used in a quality comparison. Certificate routes are distinct from graph execution.
- Both protein pairs reuse 1o7j.a. Existing long-edge flags are reported as supplied, not newly verified missing-residue diagnoses. The dataset README's older index-wise coupling convention is not used: this pilot computes discrete Fréchet coupling with repeated indices.
- No full-resolution protein experiment, Contact F1/topology claim, phase timing, exact algorithm peak-memory measurement, or population-level conclusion.

## Runtime and resources
- One attempt per configuration; no timing repetitions or cached repetitions.
- Sum of worker wall times: 1087.078 s. Recorded completed solver/soft-limit times: 9.175 s (excludes unavailable killed-worker solve times). Notebook timings are in logs/notebook_execution.jsonl.
- Maximum valid sampled worker-tree RSS: 154.69 MiB. Protein memory unavailable. Limits: 512 MiB, 20 s after warm-up, 90 s startup; joint graph soft limit 8 s, 30,000 states, 500,000 transitions, 200,000 chain edges.
- First attempts used below-normal priority; corrected remaining workers use ordinary priority. Thus these diagnostics are unsuitable for method runtime rankings.

## Planned, not executed
- 04_graph_cost.ipynb and 05_paper_summary.ipynb are not created or run in Stage 1. No paper conclusions drafted.
- Stage 2 proposal: first resolve the 18 missing configurations on the same frozen inputs, with corrected monitoring and a proposed 60 s execution cap: at most 18 worker-minutes plus startup. Do not enlarge the cohort until coverage is complete. This is a proposal only and requires a further request.
- Future runtime study: 10 measured repetitions after warm-up, increased only after observing variability; separate graph phases and certificate conditions. Different alpha settings are not repetitions.

## Additional source caveat
The pre-existing identifier audit records 1d9q.d as 325 supplied coordinates versus paper length 297, with three contiguous fragments. The supplied single curve retains gap-bridging edges. No new fragmentation or biological completeness claim is made.
