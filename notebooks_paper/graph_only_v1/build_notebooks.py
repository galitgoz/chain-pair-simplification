from graph_core import ROOT,OUT
from pathlib import Path
import nbformat as n
HERE=Path(__file__).parent
def md(s):return n.v4.new_markdown_cell(s)
def code(s):return n.v4.new_code_cell(s)
setup="""from pathlib import Path
import sys,json
ROOT=Path.cwd()
while not (ROOT/'src/cps_paper_algorithms.py').exists(): ROOT=ROOT.parent
sys.path.insert(0,str(ROOT/'notebooks_paper/graph_only_v1'))
from graph_core import OUT,METHODS,compact_results,plot_domain
import pandas as pd
from IPython.display import display
# Central protocol: all inputs, thresholds and limits are frozen in this file.
PARAMETERS=json.loads((OUT/'protocol.json').read_text())
CERTIFICATE=False
assert PARAMETERS['certificate'] is False
"""
definition="""## Definitions and protocol
Four methods: Independent continuous, Independent discrete, CPS-2F (global continuous fidelity, discrete output coupling), CPS-3F (discrete fidelity and coupling). Minimize max(|A′|,|B′|), retain input endpoints, allow unequal output lengths and repeated coupling indices. Auxiliary positions cannot be selected as output vertices.

**Every CPS call explicitly uses `certificate=False`.** Successful joint results must have raw status `optimal_configuration_graph` and observed graph execution. The software's independent and alternative-subsequence certificate shortcuts are outside this evaluation. Independent lower-bound preprocessing remains in the existing implementation and is included in total solver time; it cannot bypass the graph.

This is the same Stage 1 cohort: two existing w=16 protein pairs (22 vertices per curve) and two hurricane pairs (25 per curve), α=0.5 and 1, δ₁=α mean-edge(A), δ₂=α mean-edge(B), δ₃=unrounded input discrete Fréchet distance. No further sampling or alignment. Shared 1o7j.a and protein gap bridges remain documented limitations. Hurricane coordinates have no timestamps and retain the shared basin frame.
"""
links='[Full protocol](../../output/paper_v1/graph_only_v1/protocol.json) · [Status](../../output/paper_v1/graph_only_v1/status.md) · [Thesis methods insert](thesis_methods_protocol.md)'
def save(name,cells):
    nb=n.v4.new_notebook(cells=cells,metadata={'kernelspec':{'name':'graph-stage1','display_name':'Graph-only Stage 1','language':'python'},'language_info':{'name':'python'}})
    n.write(nb,HERE/name)
save('00_data_and_protocol.ipynb',[md('# 00 — Graph-only data and protocol\n\n**Question:** What inputs and algorithm routes are evaluated?\n\n'+links),md(definition),code(setup),
code("manifest=pd.read_csv(OUT/'input_manifest.csv')\nparameters=pd.read_csv(OUT/'pilot_parameters.csv')\ndisplay(manifest[['pair','domain','nA','nB','units','gaps_A','gaps_B']])\ndisplay(parameters[['pair','alpha','delta1','delta2','delta3','units']])"),
md('## Compatibility with the previous pilot\nPrevious files are retained unchanged. No previous row is imported as a measurement of this run. Six successful certificate-based CPS rows are explicitly excluded. The two independent methods are measured as separate baselines in the new isolated pilot.'),
code("excluded=pd.read_csv(OUT/'excluded_certificate_results.csv')\ndisplay(excluded[['pair','alpha','method','original_raw_status','compatibility']])\nprint('Configurations:',PARAMETERS['configurations'])\nprint('Limits:',{k:PARAMETERS[k] for k in ['hard_seconds','startup_seconds','memory_mib','solver_limits']})"),
md('## Remaining gaps\nThis is a pipeline pilot, not final thesis parameter selection. One execution per configuration; timings are diagnostic. [Input provenance](../../output/paper_v1/graph_only_v1/input_manifest.csv) and [selection audit](../../output/paper_v1/graph_only_v1/selection_audit.csv) were preserved from the preselected cohort. No graph-cost expansion or paper conclusions are run here.')])
save('01_correctness.ipynb',[md('# 01 — Graph-only correctness\n\n**Question:** Do the solvers execute the graph and solve the stated problem on small inputs?\n\n'+links),md(definition),code(setup),
md('## Checks\nThe revised `graph_checks.py` reuses the prior exhaustive and targeted tests. All CPS calls explicitly pass False, including warm-up. Independently feasible examples must still execute the graph. The alternative-subsequence function is replaced by a fail-on-call guard throughout the suite; zero calls are required. Endpoints, global fidelity, discrete coupling, unequal lengths, repeated indices, infeasibility, and original output membership are checked. Continuous-fidelity enumeration shares the existing decision oracle, so finite testing is not a proof.'),
code("validation=pd.read_csv(OUT/'validation.csv')\ndisplay(validation.groupby(['affected_method','outcome']).size().rename('checks').to_frame())\ndisplay(validation.loc[validation.outcome.ne('PASS'),['check','affected_method','outcome','failure_reason']])\nassert validation.outcome.eq('PASS').all()\nprint(json.loads((OUT/'shortcut_guard.json').read_text()))"),
md('## Evidence and remaining gaps\n[Every check](../../output/paper_v1/graph_only_v1/validation.csv). The independent discrete baseline internally constructs single-curve graphs; these are preprocessing and must not be counted as entry into the joint product graph. The observer counts only direct chain-graph construction calls from `solve_cps`. Resource-limited pilots remain missing rather than infeasible.')])
for domain,num,plural in [('protein','02','proteins'),('hurricane','03','hurricanes')]:
    save(f'{num}_{plural}.ipynb',[md(f'# {num} — Graph-only {plural}\n\n**Question:** What quality and resource outcomes does direct joint graph execution produce?\n\n'+links),md(definition),code(setup+f"\nDOMAIN='{domain}'\n"),
code("parameters=pd.read_csv(OUT/'pilot_parameters.csv')\ndisplay(parameters.loc[parameters.domain.eq(DOMAIN),['pair','nA','nB','alpha','delta1','delta2','delta3','units']])"),
md('## Execution and saved results\nThe isolated runner executes one fresh attempt per configuration. This notebook reads its saved results; reopening it does not produce another timing repetition. The baseline for vertex cost uses the matching fidelity metric. Successful CPS rows must satisfy graph-route assertions; preprocessing and graph resource exits are recorded separately.'),
code("results=pd.read_csv(OUT/f'{DOMAIN}_results.csv')\njoint=results[results.method.str.startswith('CPS') & results.status.eq('optimal')]\nassert joint.raw_status.eq('optimal_configuration_graph').all()\nassert joint.configuration_graph_entered.all()\nassert not results.certificate_enabled.any()\ndisplay(compact_results(results))\ndisplay(results.groupby(['status','route']).size().rename('configurations').to_frame())\ndisplay(joint[['pair','alpha','method','chain_graphs_completed','states_discovered','transitions_examined']])"),
md(f'[Full results, fidelity, coupling, route and resource fields](../../output/paper_v1/graph_only_v1/{domain}_results.csv) · [Output indices and coupling paths](../../output/paper_v1/graph_only_v1/{domain}_paths.json)'),code('plot_domain(DOMAIN)'),
md('## Supported observations and limitations\nOnly completed feasible graph outputs are included in joint quality comparisons. A resource exit establishes neither infeasibility nor an optimal vertex count. Small horizontal plotting offsets reveal equal values without changing measurements. Protein sampling/gaps and repeated structures limit interpretation; no topological preservation or thesis contribution is inferred. Phase timing, larger cohorts, and later stages remain out of scope.')])
print('Created separate graph-only notebooks 00–03.')
