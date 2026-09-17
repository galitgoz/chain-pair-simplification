# Chain pair simplification

Research code, experiment notebooks, and saved experimental records for continuous and discrete chain pair simplification.

## Start here

- `cps_paper_algorithms.py` and `curve_algorithms.py`: core algorithms.
- `notebooks_paper/`: experiment orchestration, correctness checks, reporting, and notebooks. Subdirectories retain the historical experiment stages.
- `data/`: supplied hurricane tracks and protein inputs/prepared data.
- `output/paper_v1/`: saved experimental records, configurations, validations, tables, and figures. These are historical observations, not new executions.
- `requirements-analysis.txt`: pinned analysis environment requirements.
- Existing `README-*.md` and stage READMEs document earlier workflows and may describe an earlier stage rather than the latest analysis.

## Environment and use

Create a separate Python environment and install the dependencies from `requirements-analysis.txt`. The saved experiment environment used Python 3.14. Additional runner dependencies may be required; consult imports and the saved environment audits. This snapshot has not yet been tested by rerunning the full experiment suite on a clean machine.

Open this directory as the project root in PyCharm. Start by reading saved results and notebook code. Do not execute all notebooks indiscriminately: historical runners can write to fixed output paths. Use a fresh output directory for new experiments. Notebook embedded outputs were cleared for version control; external result records and figures are retained.

Some historical helper scripts and provenance records contain original absolute paths. Those records describe the source machine; adapt operational paths before rerunning an affected script.

## Interpretation

IC is the independent continuous baseline. CPS-2F uses continuous fidelity to each input and discrete coupling between simplified outputs; CPS-3F uses discrete distances for all three constraints. Auxiliary points are internal matching locations, not extra output vertices. The output objective is `k = max(kA, kB)`.

The primary common comparison contains 43 pairs (12 protein and 31 hurricane pairs) with validated outputs from all methods. Historical incomplete attempts are retained for provenance and are not part of that completed comparison. Timing repetitions are separate observations, not additional distinct pairs.

## Snapshot and data provenance

`snapshot_manifest.json` records the initial copied file hashes. The original working project was left unchanged. See data READMEs and per-run manifests for source provenance. Source data and third-party materials retain their respective rights. No redistribution license is asserted by this snapshot.
