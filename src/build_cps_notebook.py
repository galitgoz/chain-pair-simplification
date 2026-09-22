from pathlib import Path
import nbformat as nbf

root=Path(__file__).resolve().parents[1]
nb=nbf.v4.new_notebook();cells=[]
def md(s):cells.append(nbf.v4.new_markdown_cell(s))
def code(s):cells.append(nbf.v4.new_code_cell(s))

md('''# Protein backbones: global GCS, CPS-3F and thesis CPS-2F

This notebook compares three **different** optimization problems on decimated, ordered 3D Cα curves. Coordinates, arc lengths and distance thresholds are in **ångströms (Å)**. It runs every unique pair from Fan et al. §7, Tables 1–3 (14 pairs; Tables 1 and 2 reuse seven pairs).

The executable default is `w=16`, a **coarse-trace/runtime experiment**, with `α ∈ {1,2,4}`. This is not a biological headline result. Change `WS` below to run other supplied levels, including the `w=1` control. Strides 8 and 16 must remain labelled runtime-only. No additional decimation is applied to the supplied CSV curves.

**Inputs:** by default, use every vertex of the supplied **aligned** decimated CSV curves, including their original endpoints and any gap-crossing edges. Tables 1–2 supply seven pairs. The seven Table 3 pairs are supplemented from the full cached audited chains, decimated separately with endpoints retained; unequal lengths are supported. Their cached rigid US-align transforms are held fixed. No historical result counts are treated as new measurements. Set `INPUT_SCOPE="qualified"` to use the previous continuous-interval preparation instead.

Sources: [van de Kerkhof et al., Global Curve Simplification, Algorithm 1](papers/global%20curve%20simplification.pdf); [Fan et al., full version, §4 Algorithm 1 and §7](https://arxiv.org/pdf/1409.2457); [Galit Gozoltzani, local thesis, Chapters 3–5, Algorithms 1–2](papers/M_Sc__Thesis___Galit_Gozoltzani.pdf). The thesis file is an annotated draft with unresolved editorial/proof comments. Its configuration graph is implemented as specified; small-instance tests support implementation correctness, not a new proof of its discretization theorem.''')

md('''## 1. Load data and configure the experiment

Use the project environment in `requirements-analysis.txt`. The default loads the supplied aligned CSV coordinates without cropping or another decimation. Table 3 supplements use full cached residue coordinates and a previously computed rigid alignment. The coordinate-frame provenance is recorded per pair. The curves are mathematical polylines: edges crossing missing residues remain edges in this default experiment.''')
code('''from pathlib import Path
import json, hashlib, inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, HTML, Markdown
ROOT = Path.cwd()
if not (ROOT / 'data').exists():
    raise RuntimeError('Run this notebook from the project root.')
OUT = ROOT / 'output' / 'cps_papers'
OUT.mkdir(parents=True, exist_ok=True)
WS = (16,)  # choices: (1, 2, 4, 8, 16); w >= 8 is runtime/coarse-trace only
ALPHAS = (1, 2, 4)
INCLUDE_TABLE3 = True
INPUT_SCOPE = "supplied"  # or "qualified": previously audited common continuous intervals
BUDGETS = dict(max_states=300_000, max_transitions=15_000_000, seconds=90.)
from cps_notebook_data import load_inputs
data = load_inputs(ROOT, WS, INCLUDE_TABLE3, scope=INPUT_SCOPE)''')
md('''## 2. General properties — three cells

Printed names are preserved beside the verified PDB/author-chain interpretation. A PDB-cache SHA-256 check prevents silently reusing an audit against changed structures. Current deposited lengths can differ from reported lengths; none is trimmed merely to match the paper.''')
code('''display(pd.DataFrame({
    'property': ['CSV vertices', 'CSV chain/frame/stride groups', 'CSV pair/frame/stride groups',
                 'unique selected pairs', 'executed strides', 'coordinate units'],
    'value': [len(data['curves']), len(data['chain_summary']), len(data['supplied_pairs']),
              data['inputs'].pair.nunique(), str(WS), 'Å']}))
display(data['inputs'][['B','tables','w','nA','nB','role']])''')
code('''display(data['audit'][['printed','verified','reported_length','deposited_ca_length',
    'prepared_ca_length','length_matches_paper','author_chain','label_chain','model','fragments','interpretation']])''')
code('''display(data['inputs'][['B','origin','dropped_decimated_vertices_per_chain',
    'source_first','source_last','full_chain_endpoints_preserved','fragments_A','fragments_B']])
display(data['paper'])  # historical thresholds and results, NOT newly computed results''')
md('''## 3. Organize the curves for consumption

`instances["1o7j.a__1hfj.c@16"]` holds `A`, `B_xyz` and the two residue maps. `input_index` is the solver index; `curve_index` is the original source-row index. Author residue number, insertion code, label residue number, model and chosen alternate atom are retained.

The residue maps preserve first-model Cα selections and document missing residues, alternate atom choices and fragment IDs. In `supplied` mode, gaps are reported but not removed; the exact input polyline is simplified. Source row 0 and the final source row are retained. Table 3 curves are sampled independently, so A and B can have different lengths; discrete Fréchet coupling does not require an index-wise correspondence. In `qualified` mode, the earlier common continuous intervals and their exclusions are used instead.''')
code('''instances = data['instances']
example = next(iter(instances.values()))
display(data['residues'][['instance','side','input_index','curve_index','pdb_id','model',
    'author_chain','label_chain','author_seq_id','insertion_code','label_seq_id','altloc','fragment_id']].head(10))
print('Example array shapes:', example['A'].shape, example['B_xyz'].shape)
display(Markdown('Detailed source records: [chain breaks](output/decimation/chain_breaks.csv), '
 '[missing residues](output/decimation/missing_residues.csv), '
 '[qualification exclusions](output/decimation/qualification_excluded_residues.csv).'))''')
md('''## 4. Algorithms and the comparison contract

| Method | Fidelity of A and B to their simplifications | Between the two outputs | Objective |
|---|---|---|---|
| Independent GCS | global continuous Fréchet | no joint constraint | minimize each vertex count separately |
| CPS-3F | global **discrete** Fréchet | discrete Fréchet ≤ δ₃ | minimize max(kA,kB) |
| Thesis CPS-2F | global **continuous** Fréchet | discrete Fréchet ≤ δ₃ | minimize max(kA,kB) |

All outputs are endpoint-preserving vertex subsequences of the prepared input. This is the anchored version used in the algorithms; the optional unanchored variants in the source remarks are not run. A discrete coupling may repeat an A or B index and permits `kA != kB`. Its maximum link is reported along the actual optimal coupling, not by zipping unequal arrays. `indexwise_max_link` is an additional diagnostic only when the sizes agree.

For each actual prepared `(pair,w)`, `e = (mean_edge(A)+mean_edge(B))/2`, `δ₁=δ₂=αe`, and `δ₃=ceil(d_dF(A,B))`. These are bond-scaled comparisons, not a numerical reproduction of the historical tables. The endpoint maximum is a necessary lower bound; a solver may still be infeasible because of fidelity constraints. With the current δ₃ construction the full input curves themselves supply a feasible discrete coupling.

**GCS:** the implementation realizes Algorithm 1's free-space cost envelope with exact integer link-count layers. Shortcut endpoints are input vertices, but their matched positions along the original curve are free to slide globally. A local shortcut path supplies only a feasible upper bound on the search depth. Time is O(n³K), memory O(n²K), K ≤ n−1; it is not a local shortcut-DAG solution.

**CPS-3F:** the explicit configuration graph from Fan et al. §4 Algorithm 1 has owner indices advancing by 0 or 1 and retained-vertex indices advancing by any nonnegative amount. It uses the paper's two-cost DP recurrence. This notebook does **not** claim the faster §5 Algorithm 2 implementation or its runtime bound.

**CPS-2F:** insert each original-vertex sphere's intersections with original edges as owner-only auxiliary points. These do not alter geometry and cannot be selected in the output. Construct the Chapter 5 graph using monotone transitions and continuous subcurve/segment tests. The dynamic program stores the minimum B-hop count at each exact A-hop count; the equivalent exact-count convention avoids conflating an 'at most' budget with initialization. Exclude self-loops, retain the source initialization, and add one to hop counts to report vertices.

For both joint methods we first compute independent minima under the matching fidelity measure. If that pair obeys δ₃, it attains a lower bound on the joint objective and is certified optimal immediately. For CPS-2F a bounded search may also find a different feasible pair at that same lower bound; this is labelled `optimal_alternate_certificate`, not graph execution. Otherwise the explicit configuration graph runs. `resource_limit` is never called infeasible and never assigned fabricated output metrics. The graph solver optimizes vertex count, not the smallest achievable pairing distance.''')
md("### 4.1 Global continuous Fréchet: free-space geometry and GCS\n\nThese executable definitions implement the global cost-layer sweep. Auxiliary functions include a continuous decision predicate and a discrete distance evaluator.")
geometry = (root/'src/curve_algorithms.py').read_text(encoding='utf-8')
geometry = geometry[:geometry.index('def cps_discrete_reference')]
code(geometry.replace('cache=True', 'cache=False'))
md("### 4.2 Joint solvers: discrete and continuous configuration graphs\n\nThe two graph builders share the two-count dynamic program. `certificate=False` forces graph execution; the default also accepts a feasible pair attaining the independent lower bound.")
joint = (root/'src/cps_paper_algorithms.py').read_text(encoding='utf-8')
joint = joint.replace('from curve_algorithms import (ball_segment_interval, segment_frechet_decision,\n    continuous_decision, discrete_frechet, shortest_independent)', '# Geometry functions are defined in the preceding cell.')
code(joint.replace('cache=True', 'cache=False'))

md('''## 5. Validate before measuring

The tests enumerate all anchored subsequences on independent small random 3D curves and compare both joint optima against their respective discrete/continuous distance predicates. They force the graph path, bypassing the independent certificate. A collinear zero-tolerance example requires three vertices under discrete fidelity and only two under continuous fidelity. Compilation in this cell is excluded from experiment timings.''')
code('''import verify_cps_papers as verification
verification.solve_cps = solve_cps
verification.shortest_independent = shortest_independent
verification.independent_discrete = independent_discrete
verification.continuous_decision = continuous_decision
verification.discrete_frechet = discrete_frechet
checks = verification.run_checks()
display(checks)''')
md('''### Exercise both configuration graphs on a protein display instance

This four-vertex window is drawn from the already decimated first pair. Disable both certificate shortcuts to exercise the actual graphs, including the thesis auxiliary points. This is a solver illustration, not an additional full-pair result; no display-instance number enters the comparison CSV.''')
code('''window_A = example['A'][:4].copy()
window_B = example['B_xyz'][:4].copy()
window_delta = float(.5*(np.linalg.norm(np.diff(window_A,axis=0),axis=1).mean()
                         + np.linalg.norm(np.diff(window_B,axis=0),axis=1).mean()))
window_delta3 = float(np.ceil(discrete_frechet(window_A, window_B)))
demonstrations = {}
display(pd.DataFrame([dict(pair=example['pair'], window='first four supplied vertices',
                          delta1=window_delta, delta2=window_delta, delta3=window_delta3)]))
for method in ('CPS-3F','CPS-2F'):
    demonstrations[method] = solve_cps(window_A, window_B, window_delta, window_delta,
                                       window_delta3, method, certificate=False)
display(pd.DataFrame({method: {k:v for k,v in result.items() if not isinstance(v,np.ndarray)}
                      for method,result in demonstrations.items()}).T)''')
md('''## 6. Run all selected pairs

The following cell recomputes results; it does not load a prewritten result table. All requested method rows survive even when a configured resource budget is reached. Certificates and configuration-DAG executions have distinct status labels.''')
code('''import cps_notebook_run as experiment
# Run the solver definitions above; edits to those cells affect this experiment.
experiment.shortest_independent = shortest_independent
experiment.solve_cps = solve_cps
experiment.ResourceLimit = ResourceLimit
experiment.continuous_distance = continuous_distance
experiment.discrete_frechet = discrete_frechet
experiment.discrete_coupling = discrete_coupling
run_comparison = experiment.run_comparison
results, paths = run_comparison(data, OUT, ALPHAS, BUDGETS)
display(results.groupby(['method','status']).size().rename('rows').to_frame())
assert len(results) == len(instances) * len(ALPHAS) * 3
assert results[results.method != 'GCS'].query("status.str.startswith('optimal')", engine='python').pair_bound_verified.all()''')
md('''### 6.1 Compression, auxiliary points and structural topology

**Sizes:** `nA=|A|`, `nB=|B|`, `kA=|A′|`, `kB=|B′|`. Compression factor is input count / output count; percentage removed is `100(1 − output/input)`. The joint objective remains `k=max(kA,kB)`; total compression is a separate diagnostic.

**Auxiliary points:** count the unique additional curve-parameter positions in the thesis extended curves, separately for A and B. Original vertices are excluded. `|A*|=|A|+auxiliary_A` and similarly for B. These points are not vertices of A′ or B′. Counts are computed for every CPS-2F case, including certificate cases, as diagnostic preprocessing outside solver timing. `auxiliary_graph_used` indicates a completed graph solution; resource-limited attempts do not have a completed graph solution. CPS-3F adds no auxiliary points.

**Structural topology proxy:** compare nonlocal Cα contact maps before and after simplification. A contact is a distance ≤8 Å between input vertices whose polymer sequence indices differ by at least 5. An [example of an 8 Å Cα contact definition](https://pmc.ncbi.nlm.nih.gov/articles/PMC6622159/) motivates the cutoff; the sequence exclusion and evaluation correspondence here are explicit experiment choices. To compare on the same vertex grid, place each discarded vertex along its shortcut using its original subcurve arc-length fraction. This geometric interpolation is used only for scoring; it does not reconstruct missing atoms or change solver input/output.

Report contact precision, recall and F1 for each chain, and a pooled F1 from summed true positives, false positives and false negatives. Undefined denominators remain missing, not perfect scores. Contact counts are exported to show sparse evidence at coarse strides. This measures preservation of nonlocal spatial relationships; **it does not certify knot type, entanglement, or absence of self-intersections**. Lower continuous Fréchet error and greater retained arc length give complementary geometric fidelity measures. Compare both methods at the exact same pair and δ1/δ2/δ3; no general topology advantage is assumed.''')
code('''from cps_quality import enrich_results, paired_comparison
CONTACT_CUTOFF = 8.0  # Å
CONTACT_SEQUENCE_SEPARATION = 5
results = enrich_results(results, data, paths, OUT,
                         cutoff=CONTACT_CUTOFF, separation=CONTACT_SEQUENCE_SEPARATION)
paired = paired_comparison(results, OUT)
auxiliary = results[results.method == 'CPS-2F'][[
    'pair','w','alpha','delta1','delta2','delta3','nA','nB','kA','kB',
    'auxiliary_A','auxiliary_B','auxiliary_total','extended_A','extended_B',
    'auxiliary_graph_used','status']]
auxiliary.to_csv(OUT/'auxiliary_counts.csv',index=False)
display(HTML(auxiliary.rename(columns={'nA':'|A|','nB':'|B|','kA':'|A′|','kB':'|B′|',
    'extended_A':'|A*|','extended_B':'|B*|'}).to_html(index=False,float_format=lambda x:f'{x:.3f}')))''')
md('''### 6.2 Per-pair parameter table — sizes, budgets and auxiliary points

This is the main reference table for **every tested protein pair**: |A|, |B|, δ1, δ2, δ3, and the additional CPS-2F auxiliary points on A and B, including their total. The extended-curve sizes |A*| and |B*| are included. Separate blocks identify stride w and α; auxiliary counts depend on the actual fidelity thresholds, so they are not pooled across α values.

**Theoretical upper bounds** are listed beside the measured counts: `U_A = 2|A|(|A|−1)` and `U_B = 2|B|(|B|−1)`, from [Galit Gozoltzani, Chapter 3, Observation 4](papers/M_Sc__Thesis___Galit_Gozoltzani.pdf). Each vertex-centered sphere intersects each edge at most twice, assuming finite intersections. These are worst-case bounds on **added** points; they exclude original vertices and are not predictions. For example, 22 vertices give a bound of 924 added points, so the extended curve has at most 946 vertices. The measured count may be much smaller.''')
code('''from cps_parameter_table import configured_parameter_table
pair_parameters, parameter_html = configured_parameter_table(data, ALPHAS, OUT)
display(HTML(parameter_html))
display(Markdown('[Open the complete parameter table](output/cps_papers/pair_parameters.html) · '
                 '[Download all pair parameters and auxiliary counts (CSV)](output/cps_papers/pair_parameters.csv)'))''')
md('''## 7. Before/after results

The rendered table and CSV report both vertex counts, arc-length ratios, input/output discrete Fréchet distances, actual coupling links, fidelity distances, count overhead over global GCS and elapsed time. Continuous distances are decision brackets with width ≤ 10⁻⁵ Å; output tables display the upper bounds. Joint bound satisfaction verifies the solver. The empirical comparison is how vertex counts and achieved distances differ, including any independent pairing violations.''')
code('''from cps_notebook_plots import render_table, summary_plots, spatial_plots
display(HTML(render_table(results, OUT/'comparison.html')))
display(Markdown('[Full CSV](output/cps_papers/comparison.csv) · '
 '[Standalone rendered table](output/cps_papers/comparison.html) · '
 '[Input/residue map](output/cps_papers/residues.csv) · '
 '[Selected vertices and discrete couplings](output/cps_papers/paths.json)'))''')
md('''### CPS-2F versus CPS-3F: matched comparisons

Each row identifies both proteins, sampling level, all three thresholds, input sizes, and both methods' output sizes. Positive `k_saved_by_2F` means CPS-2F uses a smaller maximum chain count. Positive `extra_removal_pp_2F` means more total vertices removed (percentage points); positive `contact_F1_gain_2F` means better contact preservation. Negative `continuous_error_change_2F` means lower geometric error. These can favor different methods. Missing results are not counted as losses or ties.''')
code('''paired_columns = ['pair','w','alpha','delta1','delta2','delta3','nA','nB',
    'kA_2F','kB_2F','kA_3F','kB_3F','k_saved_by_2F','extra_removal_pp_2F',
    'contact_F1_pooled_2F','contact_F1_pooled_3F','contact_F1_gain_2F',
    'continuous_error_change_2F','status_2F','status_3F']
display(HTML(paired[paired_columns].rename(columns={'nA':'|A|','nB':'|B|',
    'kA_2F':'|A′| CPS-2F','kB_2F':'|B′| CPS-2F',
    'kA_3F':'|A′| CPS-3F','kB_3F':'|B′| CPS-3F'}).to_html(
        index=False,na_rep='—',float_format=lambda x:f'{x:.3f}')))
display(Markdown('[Matched comparisons CSV](output/cps_papers/cps2f_vs_cps3f.csv) · '
                 '[Auxiliary point counts](output/cps_papers/auxiliary_counts.csv)'))''')
code('''valid = results[results.status.str.startswith('optimal')]
gcs = valid[valid.method == 'GCS']
display(pd.DataFrame({'quantity': ['independent pairs exceeding δ₃', 'successful joint verifications',
 'resource-limited rows', 'infeasible rows'], 'value': [int((~gcs.pair_bound_verified.astype(bool)).sum()),
 int(valid[valid.method!='GCS'].pair_bound_verified.sum()),
 int((results.status=='resource_limit').sum()), int(results.status.str.startswith('infeasible').sum())]}))
sampling = results.drop_duplicates(['pair','w'])[['B','w','role','nA','nB','edge_mean_before','edge_mean','edge_growth','edge_min','edge_max']]
display(sampling)
sampling.to_csv(OUT/'sampling_geometry.csv',index=False)''')
md('''## 8. Figures: compression quality and structural fidelity

Each row compares CPS-2F with CPS-3F for one named pair. Every tick gives δ1, δ2 and δ3 in Å. Left: total vertex removal (higher is more compressed). Middle: pooled nonlocal contact F1 (higher preserves more contacts). Right: maximum continuous fidelity error across A and B (lower follows the input curves more closely). Coarse sampling can leave very few input contacts, so inspect the exported contact counts before drawing biological conclusions.''')
md('''### 8.1 Boxplots: coupling budget and vertex cost

Panel (a) shows `d_dF(A′,B′)/δ3` for **all three methods**, including independent GCS. The dashed line at 1 is the pairing budget; values above it violate that budget. Panel (b) shows the requested boxplots of `k=max(|A′|,|B′|)` for each method, with independent GCS as the baseline. A companion boxplot makes the paired extra cost explicit: `k_method − k_GCS` for the **same input pair and tolerances**.

Only cases solved by all three methods enter these plots. An unsolved case is excluded from every method, so the distributions use identical inputs. The cohort table below identifies every protein pair, stride, α and δ1/δ2/δ3, including exclusions. Each point is a parameter case; repeated thresholds for the same protein pair are not independent biological observations. Boxes show medians and quartiles, whiskers extend to the furthest observation within 1.5 IQR, and all observations are overlaid.''')
code('''from cps_quality_plots import budget_vertex_boxplots
boxplot_cohort = budget_vertex_boxplots(results, OUT/'figures')
display(HTML(boxplot_cohort.to_html(index=False,float_format=lambda x:f'{x:.3f}')))
display(Markdown('[Exact pairs and thresholds](output/cps_papers/boxplot_cohort.csv)'))''')
md('''### 8.2 Publication figure: side-by-side grouped bars

For each α and stride, the left panel shows the achieved discrete coupling divided by δ3, with a dashed constraint line at 1. The right panel compares `k=max(|A′|,|B′|)` for Independent, CPS-2F and CPS-3F. Both panels use the same full protein-pair labels, rotated 90°. Colors and hatching agree across panels. Missing solver results are marked **NA**, not zero. The requested titles are descriptive headings; the plotted vertex costs show the actual overhead, including cases where CPS-3F costs more.

Each figure is exported as vector PDF and SVG, plus 300 dpi PNG. Mathematical labels use STIX math rendering and embedded PDF fonts, without requiring a local LaTeX installation. Exact δ1/δ2/δ3 are listed by figure and pair below.''')
code('''from cps_quality_plots import academic_dual_bars
dual_bar_parameters = academic_dual_bars(results, OUT/'figures')
display(HTML(dual_bar_parameters.to_html(index=False,float_format=lambda x:f'{x:.3f}')))
display(Markdown('[Per-pair thresholds](output/cps_papers/dual_bar_parameters.csv) · '
                 '[α=1 vector PDF](output/cps_papers/figures/dual_bars_w16_alpha1.pdf) · '
                 '[α=1 PNG](output/cps_papers/figures/dual_bars_w16_alpha1.png)'))''')
code('''from cps_quality_plots import quality_plots
quality_plots(results, OUT/'figures')''')
md('''## 9. 3D before/after figures for every pair

At α=1, each pair has an input panel and one panel per algorithm. The four panels share the same viewpoint, x/y/z limits, and equal axis scales. Dotted links show the discrete Fréchet coupling, including repeated indices. Faint lines show the full prepared decimated input. Numerical measurements use 3D coordinates, never a projection.''')
code('''for w in WS:
    spatial_plots(data, paths, OUT/'figures', alpha=ALPHAS[0], w=w)''')
md('''## 10. Interpretation, limitations and planned comparisons

- Compare global GCS with CPS-2F to isolate the price of a joint constraint under the same continuous fidelity. Compare CPS-2F with CPS-3F to quantify the effect of replacing discrete fidelity by continuous fidelity. Do not interpret different-sized outputs as an unusable coupling: the papers use discrete Fréchet traversal, not an equal-size zip.
- Report what the measurements show. Independent outputs need not violate δ₃ in every run; a slack threshold and coarse input can make both methods coincide. Joint solvers enforce a threshold but do not minimize their achieved pairing distance.
- `w=16` is an initial executable coarse-trace comparison. `w=1` is available as the undecimated control; `w=2,4` are the next biologically relevant comparisons. No w≤4 headline is inferred from this default run. Explicit graph budgets may limit larger instances; raise them deliberately or implement the source's faster/sparser algorithms without changing the optimization problem.
- The mean chord grows sublinearly with stride because the backbone is folded. `edge_growth`, `edge_min` and `edge_max` make this visible. The appended final vertex creates a short stub when the stride misses the endpoint. At mean edges above roughly 10 Å, regard the curve as a coarse trace that loses helical detail.
- The default retains all supplied vertices and treats gap-crossing edges as straight segments. The optional qualified mode removes regions outside common continuous intervals. Neither mode reconstructs missing residues. Table 3 supplements retain unequal full-chain lengths, but do not recover the unavailable historical extraction.
- The supplied pairs use their provided Kabsch frame; Table 3 supplements use cached US-align transforms fitted on qualified intervals and applied to full curves. These pair-specific frames, present-day coordinates and decimation mean that historical distances and counts need not match.
- The local thesis is an annotated draft. Graph optimality and global continuous feasibility are separately tested; the asserted equivalence of its auxiliary-point graph to the full continuous problem rests on the thesis's discretization claim. Remaining source annotations should be resolved before treating it as an independently established theorem.
- Endpoints are anchored here for all three methods. Fan and the thesis discuss unanchored extensions; results for those variants should be separate. Historical Table 3 also uses unequal δ₁ and δ₂, whereas this bond-scaled comparison intentionally sets them equal.

Reproducibility: configuration, file hashes and all validation counts are recorded below. Supporting scripts are local; no network request is needed to rerun against the cached audited inputs.''')
code('''source_files = [ROOT/'src/cps_paper_algorithms.py', ROOT/'src/curve_algorithms.py',
 ROOT/'src/cps_notebook_data.py', ROOT/'src/cps_notebook_run.py', ROOT/'src/cps_quality.py',
 ROOT/'src/cps_quality_plots.py', ROOT/'src/cps_notebook_plots.py', ROOT/'src/cps_parameter_table.py',
 ROOT/'data/decimated protein backbones/decimated_curves.csv',
 ROOT/'papers/M_Sc__Thesis___Galit_Gozoltzani.pdf', ROOT/'papers/global curve simplification.pdf',
 ROOT/'papers/chain_pair_simplification.pdf']
manifest = dict(ws=list(WS), alphas=list(ALPHAS), include_table3=INCLUDE_TABLE3,
 input_scope=INPUT_SCOPE, units='angstrom', endpoints='anchored to prepared inputs', budgets=BUDGETS,
 rows=len(results), pairs=results.pair.nunique(), checks=checks,
 contact_cutoff=CONTACT_CUTOFF, contact_sequence_separation=CONTACT_SEQUENCE_SEPARATION,
 sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files})
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
display(manifest)''')
import sys
sys.path.insert(0, str(root))
from scripts.legacy.add_graph_reporting import heading as graph_heading, source as graph_source
graph_index = next(i for i,c in enumerate(cells) if c.source.startswith('## 9.'))
cells[graph_index:graph_index] = [nbf.v4.new_markdown_cell(graph_heading), nbf.v4.new_code_cell(graph_source)]
nb.cells=cells
nb.metadata.update(kernelspec=dict(display_name='Protein analysis',language='python',name='protein-analysis'),language_info=dict(name='python'))
nbf.write(nb,root/'notebooks/legacy/0.analyze_CPS.ipynb')
print(f'Created notebook with {len(cells)} cells')
