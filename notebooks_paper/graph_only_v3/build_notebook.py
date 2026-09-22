from common import ROOT,OUT
from pathlib import Path
import nbformat as n
md=n.v4.new_markdown_cell;code=n.v4.new_code_cell
nb=n.v4.new_notebook(cells=[
md('''# One controlled CPS-2F graph experiment

**Question:** Does increasing only the product-transition budget from 2 million to 10 million allow this frozen hurricane configuration to finish?

CPS-2F; AL031854__AL011991; α=0.5. Both complete supplied tracks retain all 25 vertices in their existing shared coordinate frame. No preprocessing or threshold change. One new attempt; no timing repetitions. Previous results remain unchanged.
'''),code('''from pathlib import Path
import json,pandas as pd
from IPython.display import display
ROOT=Path.cwd()
while not (ROOT/'src/cps_paper_algorithms.py').exists(): ROOT=ROOT.parent
OUT=ROOT/'output/paper_v1/graph_only_v3'
# Central parameters: frozen before the single attempt.
PARAMETERS=json.loads((OUT/'configuration.json').read_text())
assert PARAMETERS['certificate'] is False
display(pd.DataFrame([{k:PARAMETERS[k] for k in ['pair','method','nA','nB','alpha','delta1','delta2','delta3','units']}]))
print('Graph limits:',PARAMETERS['graph_limits'])
print('Separate phase budgets:',PARAMETERS['phase_budgets'])
print('Sampled worker-tree RSS cap:',PARAMETERS['memory_mib'],'MiB')
'''),md('''## Execution policy and evidence

The unchanged optimization algorithm explicitly receives `certificate=False`. Neither certificate route is allowed. The reference value is read for comparison only after the run; it is not a pruning or early-stop condition.

The v2 phase-separated worker and atomic persistence are reused. Original CPS internal feasibility checks remain in solver time. External validation and optional measurements occur after saving the raw output, under separate budgets. RSS includes the whole worker tree and is approximate process memory, not exact algorithm memory.
'''),code('''comparison=pd.read_csv(OUT/'comparison.csv')
display(comparison[['version','transition_limit','status','raw_solution_saved','validation_status','objective','reference_agreement']])
display(comparison[['version','transitions','discovered_states','states_dequeued','completed_states_lower_bound','queue_size','solver_seconds','validation_seconds','peak_sampled_tree_rss_mib']])
'''),md('''## Observed progress

Coarse snapshots are recorded at heap-pop boundaries, at most once per second. Final counters come from the completed solver. Dequeued states include the terminal state; the completed-state prefix excludes it. Counts and queue size are **not percentages of completion**. Frame inspection, clocks and log writes add unmeasured overhead, so this is a diagnostic trace, not a runtime benchmark.
'''),code('''import matplotlib.pyplot as plt
entry=json.loads((OUT/'attempt_index.json').read_text())[0]
attempt=Path(entry['attempt_directory'])
raw=json.loads((attempt/'raw_solver.json').read_text())
events=[json.loads(line) for line in (attempt/'phases.jsonl').read_text().splitlines()]
points=[dict(seconds=e['observed_solver_seconds'],transitions=e['product_transitions'],states=e['states_discovered']) for e in events if e.get('configuration_graph_entered') and 'observed_solver_seconds' in e]
points.append(dict(seconds=raw['solver_seconds'],transitions=raw['raw_result']['transitions'],states=raw['raw_result']['states']))
progress=pd.DataFrame(points).sort_values('seconds')
with plt.rc_context({'figure.facecolor':'white','axes.facecolor':'white','text.color':'black','axes.labelcolor':'black'}):
    fig,axes=plt.subplots(1,2,figsize=(9,3.3))
    axes[0].plot(progress.seconds,progress.transitions/1e6,'o-',color='#009E73')
    axes[1].plot(progress.seconds,progress.states,'o-',color='#009E73')
    axes[0].set_ylabel('Examined transitions (millions)')
    axes[1].set_ylabel('Discovered states')
    for ax in axes: ax.set_xlabel('Elapsed solver time (s)'); ax.grid(alpha=.2)
    fig.suptitle('AL031854 / AL011991; α=0.5; CPS-2F graph — diagnostic trace')
    fig.tight_layout();fig.savefig(OUT/'progress.png',dpi=180);plt.show();plt.close(fig)
'''),md('''## Supported outcome and next action

The graph returned two 5-vertex outputs and objective **k=5**, with raw status `optimal_configuration_graph`, both chain graphs completed and product-graph execution observed. The raw result was saved before external validation, which passed. It agrees with the previously verified reference optimum; the certificate-based reference remains a separate record.

Final work: **7,076,530 transitions**, **5,678 discovered/dequeued states**, **5,677 completed nonterminal states**, queue size **0**. Solver time was **8.987 s**; external validation **0.00116 s**. Maximum sampled worker-tree RSS was **223.24 MiB**. No limit stopped this attempt.

Coverage is **23 inherited + 1 new = 24/32**, including **two** successful CPS-2F configuration-graph results. This single case does not establish a general runtime or complexity conclusion.

**One next action:** a separately authorized, one-attempt transition-budget follow-up for the previously edge-repaired protein diagnostic B, retaining its other limits. Nothing further is executed here.

[Frozen configuration](../../output/paper_v1/graph_only_v3/configuration.json) · [Comparison table](../../output/paper_v1/graph_only_v3/comparison.csv) · [Attempt and checkpoints](../../output/paper_v1/graph_only_v3/attempt_index.json) · [Coverage and provenance](../../output/paper_v1/graph_only_v3/coverage_with_provenance.csv) · [Status](../../output/paper_v1/graph_only_v3/status.md) · [Instrumentation](instrumentation.md)
''')],metadata={'kernelspec':{'name':'graph-v3','display_name':'Graph v3 diagnostic','language':'python'},'language_info':{'name':'python'}})
n.write(nb,Path(__file__).with_name('00_controlled_transition_budget.ipynb'))
