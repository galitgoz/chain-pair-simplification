"""Create new Stage 1 notebooks; never touch legacy notebooks."""
from pathlib import Path
import nbformat as n
ROOT=Path(__file__).resolve().parents[1]
def md(s):return n.v4.new_markdown_cell(s)
def code(s):return n.v4.new_code_cell(s)
setup="""from pathlib import Path
import sys, json
ROOT = Path.cwd()
if ROOT.name == 'notebooks_paper': ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / 'notebooks_paper'))
from paper_core import OUT, PROTOCOL, METHODS, compact_results, plot_domain
import pandas as pd
from IPython.display import display, Markdown
# Central parameters: frozen before optimization; no notebook-specific overrides.
PARAMETERS = json.loads((OUT / 'protocol.json').read_text())
assert PARAMETERS['alphas'] == [0.5, 1.0]
"""
definitions=r"""## Inputs and definitions
Retain the first and last input vertices. Outputs are original-vertex subsequences and may have unequal lengths. Minimize $k=\max(|A'|,|B'|)$. Discrete couplings may repeat an index. Auxiliary points represent matching positions, never selectable output vertices.

| Method | Input fidelity | Output coupling |
|---|---|---|
| Independent continuous | global continuous Fréchet | measured, unconstrained |
| Independent discrete | discrete Fréchet | measured, unconstrained |
| CPS-2F | global continuous Fréchet | discrete, at most δ₃ |
| CPS-3F | discrete Fréchet | discrete, at most δ₃ |

For each pair, $s_A,s_B$ are the mean input edge lengths; $\alpha\in\{0.5,1\}$, $\delta_1=\alpha s_A$, $\delta_2=\alpha s_B$, and $\delta_3=d_{dF}(A,B)$ without rounding. All four methods receive identical inputs and thresholds. This is a pipeline pilot, not final paper parameter selection.
"""
links="[Full input manifest](../output/paper_v1/input_manifest.csv) · [Numerical thresholds](../output/paper_v1/pilot_parameters.csv) · [Status and gaps](../output/paper_v1/status.md)"
def save(name,cells):
    nb=n.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python (paper Stage 1)','language':'python','name':'paper-stage1'},'language_info':{'name':'python'}})
    n.write(nb,ROOT/'notebooks_paper'/name)
save('00_data_and_protocol.ipynb',[
 md('# 00 — Data and protocol\n\n**Research question:** What are the experiment inputs and definitions?\n\n'+links),
 md(definitions),code(setup),
 code("inventory = pd.read_csv(OUT / 'existing_inventory.csv')\nprint('Existing artifacts inventoried:', len(inventory))\ndisplay(inventory.groupby('kind').size().rename('count').to_frame())\ndisplay(pd.read_csv(OUT / 'prior_experiments.csv'))"),
 md('## Frozen cohort, selected by size before methods ran\nTwo smallest eligible aligned protein pairs at the fixed, already supplied w=16 resolution; the smallest supplied hurricane pair in each basin. Ties use identifiers. This is not all available pairs. Selection uses no method outcomes. [Complete selection audit](../output/paper_v1/selection_audit.csv).'),
 code("manifest = pd.read_csv(OUT / 'input_manifest.csv')\ndisplay(manifest[['pair','domain','nA','nB','units','extent_policy']])\ndisplay(manifest[['pair','subsampling','gaps_A','gaps_B','alignment']])\ndisplay(pd.read_csv(OUT / 'shared_entities.csv'))"),
 md('Both protein pairs reuse 1o7j.a and are statistically dependent. They are entire supplied **subsampled** curves, not full-resolution backbones. Protein gaps and prior alignment remain as supplied. Hurricanes retain the common coordinate frame within each basin; no new fit or subsampling is performed. Coordinates have no timestamps, so no temporal synchronization is inferred.'),
 code("parameters = pd.read_csv(OUT / 'pilot_parameters.csv')\nwith pd.option_context('display.precision', 9):\n    display(parameters[['pair','nA','nB','alpha','s_A','s_B','delta1','delta2','delta3','units']])\nprint('Configurations:', len(parameters) * len(METHODS))\nprint('Limits:', PARAMETERS['solver_limits'], 'hard seconds:', PARAMETERS['hard_seconds'], 'MiB:', PARAMETERS['memory_mib'])"),
 md('## Provenance and interpretation\n[Inventory with hashes](../output/paper_v1/existing_inventory.csv), [reference document hashes](../output/paper_v1/references/manifest.json), and [frozen protocol](../output/paper_v1/protocol.json) provide traceability. The reference definitions and configuration graph are used as specifications. Thesis experiment numbers and conclusions are not reused as evidence. An output contains one more vertex than its number of index advances; validation explicitly checks this bookkeeping.\n\nOne execution per configuration, after worker warm-up; runtimes are diagnostic. A resource exit does not establish infeasibility. Certificate routes are reported separately from configuration-graph execution. Later graph-cost and paper-summary stages are planned only.\n\n## Supported conclusions and remaining gaps\nThe cohort and thresholds are reproducible. Its small size, coarse protein sampling, shared structure, and absent hurricane timestamps limit generalization. See status.md for completed versus missing work.')])
save('01_correctness.ipynb',[
 md('# 01 — Correctness\n\n**Research question:** Do the implementations solve the stated problems?\n\n'+links),md(definitions),code(setup),
 md('## Validation design\nReuse pure existing exhaustive-subsequence oracles without invoking legacy runners that overwrite results. Seed 20260916; 2D and 3D curves of 3–5 vertices, plus targeted singleton, duplicate, zero-threshold, endpoint, unequal-length, repeated-index, infeasibility, auxiliary-membership, and global-versus-local regression cases. Joint comparisons force certificates **off**. Analytic continuous-distance controls and an independently enumerated discrete grid coupling check distance primitives. Exhaustive global-fidelity filtering still shares the existing continuous decision procedure; finite checks do not prove numerical robustness or the reference discretization theorem.'),
 code("from correctness_checks import check_all\nvalidation = check_all()\ndisplay(validation.groupby(['affected_method','outcome']).size().rename('checks').to_frame())\ndisplay(validation.loc[validation.outcome.ne('PASS'), ['check','affected_method','outcome','failure_reason']])\nprint(json.loads((OUT / 'validation_gate.json').read_text()))"),
 md('## Results and remaining gaps\n[Every check, outcome, reason, and duration](../output/paper_v1/validation.csv). PASS means the finite check completed successfully; INCOMPLETE means validation reached a resource limit. A failing or incomplete method is blocked from the pilot. Reference documents contain draft annotations and are not executed as instructions. This notebook verifies solution correctness, not medium-scale graph performance. The later graph-cost study must separately instrument phases and memory.')])
for number,domain,title in [('02','protein','Proteins'),('03','hurricane','Hurricanes')]:
    plural='proteins' if domain=='protein' else 'hurricanes'
    extra='These are supplied w=16 inputs, not full-resolution backbones. No contact score is interpreted as topological preservation.' if domain=='protein' else 'Coordinates remain in their common basin projection. No independent alignment, temporal synchronization, or new preprocessing is performed.'
    save(f'{number}_{plural}.ipynb',[
     md(f'# {number} — {title}\n\n**Research question:** What are the benefits and costs of joint simplification on these curves?\n\n'+links),md(definitions+'\n'+extra),code(setup+f"\nDOMAIN = '{domain}'\n"),
     code("parameters = pd.read_csv(OUT / 'pilot_parameters.csv')\ndisplay(parameters.loc[parameters.domain.eq(DOMAIN), ['pair','nA','nB','alpha','delta1','delta2','delta3','units']])"),
     md('## Execution and measurements\nOne isolated worker per configuration; warm-up excluded from solver time. Continuous fidelity is reported as a numerical lower/upper bracket; discrete fidelity and output coupling are also saved. Compression is the percentage of input vertices removed across both curves. Vertex cost is joint k minus independent k with the **same fidelity metric**; missing baselines remain explicitly missing. Peak RSS includes imports and warm-up and is not phase-specific algorithm memory. Resource limits produce missing outputs, never inferred infeasibility.\n\n`independent_certificate`: independently optimal outputs meet δ₃; `alternate_certificate`: another feasible pair attains the same independent lower bound; `configuration_graph`: graph execution obtains the solution. Certificate time is not graph time.'),
     code("from paper_core import run_pilot\nresults = run_pilot(DOMAIN)\ndisplay(compact_results(results))\ndisplay(results.groupby(['status','route']).size().rename('configurations').to_frame())"),
     md(f'[Full results and all parameters](../output/paper_v1/{domain}_results.csv) · [Output vertex indices and couplings](../output/paper_v1/{domain}_paths.json)\n\nThe paired comparison connects matched parameter configurations. The geometry example uses the first frozen pair and α=0.5. Colors identify methods consistently; faint lines are inputs.'),
     code("plot_domain(DOMAIN)"),
     md('## Supported pilot observations and gaps\nRead the saved table for this run’s coverage and objectives. Differences across α are parameter sensitivity, not independent repetitions. The coupling ratio checks the specified distance budget; it is not evidence of preserved topology. These four pairs cannot support population-level claims. No approximation, packedness study, or graph-cost experiment is run here.')])
print('Created notebooks 00–03; 04–05 remain planned, not created or executed.')
