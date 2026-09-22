"""Create a short saved-results notebook for this bounded repair stage."""
from common import ROOT,OUT
from pathlib import Path
import nbformat as n
md=n.v4.new_markdown_cell;code=n.v4.new_code_cell
cells=[md('''# Graph-only v2 — execution repair and three diagnostics

**Research question:** Which resource blocks each selected missing configuration, and can raw solutions survive later validation or measurement failure?

This stage preserves 22 verified v1 results and performs only A, B and C once each. It does not expand the cohort or change coordinates, thresholds, sampling or optimization algorithms. All CPS calls explicitly use `certificate=False`.
'''),code('''from pathlib import Path
import json, pandas as pd
from IPython.display import display
ROOT=Path.cwd()
while not (ROOT/'src/cps_paper_algorithms.py').exists(): ROOT=ROOT.parent
OUT=ROOT/'output/paper_v1/graph_only_v2'
# Central parameter cell: the exact limits and thresholds committed before execution.
PARAMETERS=json.loads((OUT/'diagnostics.json').read_text())
assert len(PARAMETERS)==3 and all(c['certificate'] is False for c in PARAMETERS)
display(pd.DataFrame(PARAMETERS)[['diagnostic','pair','method','alpha','delta1','delta2','delta3']])
display(pd.DataFrame([dict(diagnostic=c['diagnostic'],**c['graph_limits'],external_solver=c['phase_budgets']['solver'],warmup=c['phase_budgets']['warmup'],validation=c['phase_budgets']['validation'],memory_mib=c['memory_mib']) for c in PARAMETERS]))
'''),md('''## Preserved inputs and reconciled evidence

The v1 files contain 100 passing checks, 22 verified results, 8 graph-based CPS successes, and 10 missing results: five transition limits, two chain-edge limits, and three unsplit worker timeouts. No newer graph-result version was found. The 22 inherited rows keep their original provenance and are not new executions or timing repetitions.

Proteins remain at their existing w=16 resolution; original gap bridges and shared chains remain unchanged. Hurricane tracks retain their common basin coordinates. [Frozen input manifest](../../output/paper_v1/graph_only_v2/input_manifest.csv) · [Reconciliation](../../output/paper_v1/graph_only_v2/reconciliation.json) · [Inherited results](../../output/paper_v1/graph_only_v2/inherited_verified_results.csv)

## Phase boundaries

Input loading → warm-up → solving → **atomic raw checkpoint** → external constraint validation → optional measurements.

The unchanged CPS `finish()` already checks endpoints, ordering, fidelity and discrete coupling **inside solver time**. The new external validation occurs only after the checkpoint exists. Validation checks feasibility, not optimality by itself; solver optimality claims and small exhaustive objective checks remain distinct.

Warm-up has its own 120 s allowance and no internal CPS time cap. Measured CPS solves use a 55 s internal soft cap and a 60 s external watchdog. Validation has 30 s. Independent discrete exposes no internal time limit, and its existing graph default is unchanged. All workers retain the 512 MiB sampled tree-RSS cap; it is approximate process memory, not exact algorithm memory.
'''),code('''regression=pd.read_csv(OUT/'regressions.csv')
display(regression.groupby('outcome').size().rename('repair checks').to_frame())
assert regression.outcome.eq('PASS').all()
display(regression[regression.check.str.startswith('checkpoint_')][['check','outcome']])
'''),md('''## Three diagnostic outcomes

The 21 repair checks include exhaustive small-input objectives and raw-checkpoint survival after validation failure, validation timeout and optional-measurement failure. These synthetic tests are not pilot repetitions.

The cells below read saved results only. They never launch or retry optimization.
'''),code('''summary=pd.read_csv(OUT/'diagnostic_summary.csv')
display(summary[['diagnostic','warmup_completed','solver_returned_solution','raw_solution_saved','validation_status','objective','reference_agreement']])
display(summary[['diagnostic','solver_seconds','validation_seconds','limiting_phase','limiting_resource','limiting_counter','peak_sampled_tree_rss_mib']])
'''),md('''## Interpretation and stopping point

A still reached its increased transition budget at 2,000,001 transitions. B completed construction under the increased edge allowance and then reached the unchanged product-transition limit at 500,001. Neither returned a solution; their exit checkpoints are preserved and validation was not run.

C returned original-vertex subsequences of sizes 20 and 19, saved the raw output, and passed validation, with objective 20. Its external validation time is separate from its solver time. This does not identify which phase caused its earlier unsplit timeout.

Saved evidence verifies reference optimum 5 for A and B, but no new completed graph objective exists to compare against it. The A reference remains certificate-based evidence, not a graph run. C has no matching verified prior reference in the frozen-protocol versions.

Coverage is **22 inherited + 1 newly validated = 23/32**. Nine configurations remain missing. Use the repaired phase boundaries for future timeout investigations. A needs transition-work diagnosis; B has moved from a construction cap to a transition cap. More time alone does not remove either count limit. Any further count increase or retry requires a separate request.

[Three-row full summary](../../output/paper_v1/graph_only_v2/diagnostic_summary.csv) · [Reference evidence](../../output/paper_v1/graph_only_v2/references.json) · [Attempt index and raw/log directories](../../output/paper_v1/graph_only_v2/attempt_index.json) · [Status and limitations](../../output/paper_v1/graph_only_v2/status.md)

No further attempts, timing repetitions, larger experiment or final paper conclusions are included.
''')]
nb=n.v4.new_notebook(cells=cells,metadata={'kernelspec':{'name':'graph-v2','display_name':'Graph v2 diagnostics','language':'python'},'language_info':{'name':'python'}})
n.write(nb,Path(__file__).with_name('00_repair_and_diagnostics.ipynb'))
