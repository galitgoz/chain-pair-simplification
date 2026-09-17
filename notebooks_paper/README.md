# Paper experiments — Stage 1 only

Start with `00_data_and_protocol.ipynb`, then `01_correctness.ipynb` and the two domain notebooks. Saved tables and the current completion status are under `../output/paper_v1/`. Notebooks 04 and 05 remain planned.

## Shared code

- `paper_core.py`: frozen protocol, input adapters, method dispatch, isolated worker orchestration, measurement and plots. Optimization remains in the unchanged root modules `curve_algorithms.py` and `cps_paper_algorithms.py`.
- `correctness_checks.py`: exhaustive small-input and targeted checks; reuses pure legacy oracle functions without calling their result-writing runners.
- `pilot_worker.py`: warm-up, one solver execution, independent constraint validation, saved output indices/coupling.
- `inspect_stage1.py` / `audit_setup.py`: initial inventory and preservation snapshots. Do not rerun them over an existing audit snapshot.
- `execute_stage1.py`: executes only explicitly named notebooks 00–03 in a separate kernel.
- `finalize_stage1.py`: reads saved results, audits preservation and writes the Stage 1 report; no optimization.

## Reproduction and preservation

Exact input arrays, full-precision thresholds, source hashes, algorithm hashes, environment versions, output indices, validation outcomes and logs are saved under `output/paper_v1`. The pilot cohort was frozen before optimization. Read saved results to inspect this run. `run_pilot` refuses to overwrite an existing domain result table; a new experiment should use a separately versioned output directory after authorization. Cached reads are not new timing repetitions.

The initial protein worker monitor did not include the Windows Python child process. Those failed attempts were preserved, invalid memory values removed, and no optimization retried. The corrected hurricane monitor includes descendants and stops its own worker tree on a resource cap. The existing user Jupyter kernel is excluded. These single executions are diagnostic, not a runtime study.

No additional preprocessing, greedy algorithm, approximation experiment, graph-cost run, or paper conclusion is introduced.
