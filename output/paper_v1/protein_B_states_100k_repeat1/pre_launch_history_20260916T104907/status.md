# Protein CPS-2F slowdown investigation — repeat not executed

## Decision
No experiment worker was active at the initial process check; only the pre-existing Jupyter kernel was present. Destination had no previous attempt. No experiment was launched. The prepared repeat configuration is exactly equal to protein_B_states_100k/configuration.json; it was saved for review only.

Initial Windows samples were mostly light (4–8% CPU, one 37% sample), but the final pre-run six-sample interval changed the assessment: aggregate CPU 47.7–100%, mean 73.98%. PyCharm PID 11492 consumed 131.3–704.8% of one core across all samples (1.3–7.0 cores). Antivirus, Zoom, browsers and system activity also appeared. Available memory remained 2418–2649 MiB. This is competing computation; the repeat was withheld. A simple automated mean-CPU >80% gate initially passed, but manual inspection of the full samples overrode it before any launch. No processes were terminated, notebooks closed, caches cleared or settings changed.

## Observed historical differences
The two historical configurations differ only in max_states (30,000 vs 100,000). Input SHA-256, thresholds, tolerances (unchanged algorithm code), certificate=False and all other limits match. Worker, supervisor, common orchestration/instrumentation, warm-up and algorithms have identical recorded/current hashes. Only controlled_followup.py, the outer preparation/reporting script, changed to support the state-budget follow-up; it does not execute the optimization in the worker.

Both runs use project .venv/Scripts/python.exe according to saved launch/provenance and Python 3.14 cache suffixes. Historical per-attempt executable hashes and complete package environments were not recorded, so runtime/package identity cannot be proven retrospectively. Current interpreter, executable hash and versions are saved in environment.json. No observed unexplained core-code or parameter mismatch was found.

Both historical runs use fresh output-local NUMBA_CACHE_DIR and MPLCONFIGDIR, with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1, the same synthetic warm-up and unlimited internal warm-up time under its separate 120 s external cap. Cache inventories and modification times are saved in investigation.json. No historical cache was touched. The prepared repeat retains that cache policy but has never warmed or run.

Warm-up rose from 20.898 s to 91.761 s. Resource-monitor gaps grew from mean 0.0738 s / median 0.0698 s / maximum 0.307 s to mean 0.3247 s / median 0.1108 s / maximum 5.709 s, despite the same nominal 0.05 s sampling sleep. Timing boundaries and monitoring code are unchanged. The slow run saved an internal wall-time exit after 56.100 s, then the supervisor reached the external solver allowance at recorded phase elapsed 62.711 s and terminated the exiting worker. Both timeout records are retained.

## Diagnosis and uncertainties
The slowdown extends beyond graph work to warm-up and monitoring responsiveness. Scheduling, competing load, paging or I/O are possible explanations, not established causes. The current competing PyCharm process does not prove it caused the earlier slowdown. Historical host CPU load, worker CPU time, disk/paging activity and clock frequency were not captured. Increasing max_states is not an established cause; the first 100,000-state run stopped with only 17,567 discovered states, and warm-up does not use that experimental state cap.

The only extra monitoring in this investigation was read-only host sampling before any potential execution; no instrumentation was added to solver loops or to the existing worker resource monitor. Current samples include some investigation overhead. No general runtime claim follows from these diagnostic records.

## Results and coverage
See comparison.csv for both prior attempts and the repeat marked not_executed_host_load. Neither prior attempt returned a solution: kA, kB, k and validation times remain missing. Queue and processed-state lower bounds are not completion percentages. Failed-run times are time to termination, not solution time. The repeat has no solver outcome or resource-limit claim because it never started. Reference k=5 was not supplied to any solver or substituted for a missing output.

Coverage remains 24/32, entirely inherited; zero new attempts and zero new verified results. Historical artifacts, including caches, were hash-checked unchanged. Configuration, input, code, cache and environment provenance are saved alongside host samples and decision.json. Only controlled_followup.py was minimally adapted for an identical repeat and the new investigate_repeat.py was added in the existing code directory. No modules were copied and no notebook was changed.

## One next step
Wait until competing work subsides, then request one repeat with the same frozen configuration. No additional sampling, repeat, resource escalation or recommendation was executed after this decision.
