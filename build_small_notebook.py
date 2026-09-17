from pathlib import Path
import nbformat as nbf
root=Path(__file__).resolve().parent
cells=[]
def md(s): cells.append(nbf.v4.new_markdown_cell(s.strip()))
def code(s): cells.append(nbf.v4.new_code_cell(s.strip()))

md(r"""
# Small protein curves: independent versus joint simplification

**Question:** can simplifying two structurally superimposed backbones independently move the resulting curves farther apart than their originals? We compare 14 protein pairs using **contiguous 30-Cα windows**, then sweep $\delta_1=\delta_2\in\{2,4,8\}$ Å. The complete pipeline has one row per pair and threshold (42 rows); no example is selected based on whether it supports the headline.

**Definitions fixed for this notebook:**

- **Independent:** shortest endpoint-preserving vertex subsequence under the **local per-shortcut continuous Fréchet criterion** requested here. This differs from notebook 0's global criterion.
- **CPS:** exact **synchronous, equal-cardinality, endpoint-preserving discrete CPS**: both fidelity constraints are discrete Fréchet and corresponding retained vertices must be within $\delta_3$. This stricter variant makes the requested index-wise CPS columns meaningful. Fan et al.'s general CPS permits unequal sizes and a non-diagonal coupling; we do not attribute this stricter variant to them.
- **Threshold:** `DELTA3_MODE='exact'` uses $\delta_3=d_{dF}(A,B)$, as requested in the column definitions. `paper_ceil_delta3` records the paper's $\lceil d_{dF}(A,B)\rceil$ separately. Switching to `'ceil'` reproduces that threshold rule, but only guarantees distance $\leq\lceil d_{dF}(A,B)\rceil$, not $\leq d_{dF}(A,B)$.
- **Originals** always means the selected, structurally aligned, **unsimplified 30-point windows**. Their sizes, arcs and distances are not whole-protein measurements.

The requested equal-size and threshold conventions resolve ambiguities in the brief; both are visible and editable below. **CPS-2F is reserved for later.** All coordinates and distance/arc columns use Å; arc ratios are dimensionless. Run from this project root. The notebook includes executable algorithm definitions and executed results.

## 1. Load the data
""")
code("""
from pathlib import Path
import json, time, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, HTML
from numba import njit
from small_pair_data import load_prepared, prepare_small_pairs
from curve_algorithms import segment_frechet_decision, discrete_frechet
from small_pair_solvers import arc_length

ROOT=Path.cwd()
assert (ROOT/'data/prepared').is_dir(), 'Run from the project root; prepared residue maps are required.'
OUT=ROOT/'output/small_pairs'; FIGURES=OUT/'figures'
OUT.mkdir(parents=True,exist_ok=True);FIGURES.mkdir(parents=True,exist_ok=True)
chains=load_prepared(ROOT)
audit=pd.read_csv(ROOT/'output/identifier_audit.csv')
breaks=pd.read_csv(ROOT/'output/chain_breaks.csv')
N_TARGET=30
DELTAS=(2.,4.,8.)
DELTA3_MODE='exact'  # Alternative: 'ceil'; both exact and ceiling are always reported.
RECOMPUTE_ALIGNMENT=False
EPS=1e-9
plt.rcParams.update({'font.size':10,'figure.dpi':110})
""")
md('## 2. General properties\nThree short cells describe inventory and provenance. No arc lengths or inter-chain distances are computed before structural superposition.')
code("""
display(audit[['printed','verified','reported_length','prepared_ca_length','fragments']])
print('15 selected chains; 14 unique paper pairs (Tables 1 and 2 reuse the same seven pairs).')
print('Printed identifiers and historical length discrepancies remain unchanged from notebook 0.')
""")
code("""
display(audit[['verified','author_chain','label_chain','model','alternate_ca_residues','selected_noncanonical_altloc']])
print('Residue maps preserve model, author and label chain/sequence IDs, insertion codes, and alternate locations.')
print('The raw files are preserved; deposited PDB records supply mappings. Missing regions remain fragment boundaries.')
""")
code("""
inventory=pd.DataFrame([dict(chain=key,fragments=len(c['fragments']),
    longest_fragment=max(len(f['xyz']) for f in c['fragments']),target_window=N_TARGET)
    for key,c in chains.items()])
display(inventory)
print('Select the longest continuous fragment per chain; ties use the earliest fragment. No resampling or gap bridging.')
""")
md(r"""
## 3. Organize smaller curves and superimpose structurally

1. Use the longest continuous fragment of each chain. The reference is interpreted **1O7J A**, preserving printed `107j.a` in the results. Other printed labels are preserved in the audit and pipeline.
2. Run **US-align in rigid, sequential, sequence-independent structural alignment mode** on the two parent fragments. This searches for residue correspondence; it is not index-wise Kabsch.
3. Find the middle residue pair in that structural correspondence. Center one **contiguous 30-residue window** around each anchor (clamp at a fragment boundary). The reference A window is therefore **pair-specific**. Selection does not use simplification results or a target distance.
4. Structurally align the two selected windows with US-align again, then freeze the resulting rigid B-to-A transform for every method and threshold. Only now compute analysis distances and arc lengths. Retain every residue in each selected window, including structurally unaligned residues; deleting outliers would change the problem.

`pairs[pair_id]` contains aligned arrays `A`, `B`, residue maps `map_A`, `map_B`, parent fragment/window indices, and the recorded transform/correspondence. Each side has shape `(30,3)`. The same original arrays are used for independent and joint methods. US-align maximizes a structural similarity objective heuristically; it does not globally minimize discrete Fréchet distance, and poor fold agreement may remain.

These reduced curves cannot reproduce the paper's whole-protein distances or simplification sizes. We do not force dissimilar pairs toward 20–21 Å. Parent/window alignment logs, matrices, executable version/hash and selection metadata are cached for inspection and offline reruns.
""")
code("""
chains,pairs,pair_metadata=prepare_small_pairs(ROOT,n=N_TARGET,recompute=RECOMPUTE_ALIGNMENT)
alignment_rows=[]
for key,p in pairs.items():
    assert len(p['A'])==len(p['B'])==N_TARGET
    for side in ['A','B']:
        assert np.all(np.diff(p[f'map_{side}'].label_seq_id)==1)
    alignment_rows.append(dict(pair=key,A_fragment=p['A_fragment'],B_fragment=p['B_fragment'],
        A_author_start=p['map_A'].iloc[0].author_seq_id,A_author_end=p['map_A'].iloc[-1].author_seq_id,
        B_author_start=p['map_B'].iloc[0].author_seq_id,B_author_end=p['map_B'].iloc[-1].author_seq_id,
        window_TM=p['window_alignment']['tm_scores'][0],aligned_RMSD=p['window_alignment']['aligned_rmsd']))
alignment_audit=pd.DataFrame(alignment_rows)
display(alignment_audit.round(3))
alignment_audit.to_csv(OUT/'alignment_audit.csv',index=False)
usalign_source=json.loads((ROOT/'tools/usalign/source.json').read_text())
print(usalign_source['version'])
""")
code("""
first_pair=next(iter(pairs.values()))
display(first_pair['map_A'][['window_index','curve_index','fragment_id','label_seq_id',
    'author_seq_id','insertion_code','residue_name','model','altloc']].head(6))
print('map_A/map_B coordinates retain deposited values for provenance; A/B arrays contain the coordinates used for analysis.')
print('B_aligned = B_deposited @ rotation_B_to_A.T + translation_B_to_A')
""")
md(r"""
### The prepared 30-point pairs in 3D
Blue is A; orange is B. Both curves receive the **same display translation**, preserving their superposition. All pairs use the same two viewpoints, equal x/y/z scales, and common numerical axis limits. Circles mark starts and squares mark ends. No curve is individually recentered for these plots.
""")
code("""
from small_pair_plots import plot_aligned, pipeline_html, plot_outcomes, plot_case
plot_aligned(pairs,FIGURES)
""")
md(r"""
## 4. Algorithms and feasibility

### Shortest independent simplification: the requested LOCAL criterion
Build a DAG with all original vertices. Edge $i\to j$ is admissible precisely when
$$d_F(P[i..j],\overline{P_iP_j})\leq\delta.$$
A shortest path from vertex 0 to vertex $n-1$ minimizes the number of retained vertices. The segment decision computes, for each subcurve vertex, its allowed interval on the shortcut segment and checks for a monotone traversal; convexity covers the full original edges. Endpoints are retained. Ties use the first predecessor found in increasing index order. Counts are optimal for this **local** criterion, not claimed optimal under notebook 0's global criterion.
""")
source=(root/'curve_algorithms.py').read_text()
local=source[source.index('@njit(cache=True)\ndef local_shortcut_path'):source.index('@njit(cache=True)\ndef _gcs')]
local=local.replace('"""Restricted local baseline; used ONLY as a feasible link-count upper bound."""','"""Shortest endpoint-preserving path under the local shortcut criterion."""')
code(local.replace('cache=True','cache=False'))
md(r"""
### Exact synchronous discrete CPS and its feasibility floor

Minimize $k=|A'|=|B'|$ over increasing, endpoint-preserving subsequences such that
$$d_{dF}(A,A')\leq\delta_1,\quad d_{dF}(B,B')\leq\delta_2,\quad\max_t\|A'_t-B'_t\|\leq\delta_3.$$
The final bound supplies a valid diagonal discrete Fréchet coupling, so $d_{dF}(A',B')\leq\delta_3$. This is a **solver invariant**, not an empirical discovery.

The dynamic-programming state `(i,j,p,q)` gives the original positions `(i,j)` and retained positions `(p,q)`. Each original advances by zero or one; the retained indices either both stay or both increase strictly, which adds one paired vertex. Prefix minima over `(p,q)` avoid enumerating every jump. Its time and memory are $O(n_A^2n_B^2)$ for this synchronous variant; this is **not** the complexity or recurrence of Fan et al.'s general Algorithm 2.

Run the same DAG with a minimax objective to obtain the **exact feasibility floor**: the smallest achievable maximum paired distance satisfying the two discrete fidelity constraints. It depends on the pair, fidelity threshold, endpoint policy, and synchronous coupling rule. If `delta3 < floor` (within a `1e-9` Å comparison tolerance), record **infeasible**, leave joint output fields empty, and exclude that row from all aggregate comparisons. A floor can exceed original dDF because an optimal discrete Fréchet coupling need not be diagonal. The endpoint-distance lower bound alone is insufficient.

General Fan CPS with threshold equal to the original dDF always has the originals as a feasible candidate. The nontrivial infeasible rows here arise from the additional synchronous restriction; they are not failures of Fan CPS or failed solver runs.
""")
source=(root/'small_pair_solvers.py').read_text()
code(source[source.index('@njit(cache=True)'):].replace('cache=True','cache=False'))
md('### Validate before running proteins\nThe suite compares local shortest paths, joint cardinality optima, and minimax floors against exhaustive enumeration on 60 small random instances, including unequal lengths and explicit infeasibility.')
code("""
from verify_small_pairs import run_tests
test_report=run_tests()
# Also exercise the notebook-defined functions.
toy=np.array([[0.,0,0],[1,0,0],[2,0,0],[3,0,0]])
assert len(local_shortcut_path(toy,0))==2
assert synchronous_cps(toy,toy,0,0,0)['k']==4
assert synchronous_cps(toy,toy,1,1,0)['k']==2
""")
md(r"""
## 5. Run the full pipeline

One row per **unique pair × fidelity threshold**. Compute original dDF once after alignment; reuse its exact value (or ceiling if configured) for every threshold. For independent equal-size outputs, compute the index-wise maximum exactly. When sizes differ, write **`no coupling, |A'| != |B'|`** and leave the index-wise numeric/check fields empty; the discrete Fréchet distance is still reported in its own column, never substituted for that missing index-wise quantity.

For CPS, verify both fidelity bounds, index-wise bound and discrete Fréchet bound. The requested `cost_of_joining` is $k_{CPS}-\max(k_A,k_B)$. **It includes a change from continuous local fidelity to discrete fidelity, plus synchronization**, so it is not an isolated causal estimate of the cost of joint optimization.
""")
pipeline=(root/'small_pair_pipeline.py').read_text()
code(pipeline[pipeline.index('def run_pipeline'):pipeline.index("if __name__=='__main__':")])
code("""
pipeline,selected_paths,aggregates=run_pipeline(pairs,OUT,thresholds=DELTAS,delta3_mode=DELTA3_MODE)
assert len(pipeline)==len(pairs)*len(DELTAS)==42
assert not pipeline.duplicated(['pair','delta1']).any()
uncoupled=~pipeline.separate_coupled
assert pipeline.loc[uncoupled,'separate_index_max'].isna().all()
assert pipeline.loc[uncoupled,'separate_index_le_delta3'].isna().all()
feasible=pipeline[pipeline.included_in_aggregates]
assert feasible.CPS_solver_verified.all()
assert (pipeline.loc[~pipeline.included_in_aggregates,'delta3']+EPS < pipeline.loc[~pipeline.included_in_aggregates,'feasibility_floor']).all()
if DELTA3_MODE=='exact':
    assert np.allclose(pipeline.delta3,pipeline.original_dDF,rtol=0,atol=0)
    assert (feasible.CPS_dDF<=feasible.original_dDF+EPS).all()
print(f'{len(pipeline)} rows: {len(feasible)} feasible; {len(pipeline)-len(feasible)} below their exact feasibility floors.')
""")
md(r"""
## 6. The complete comparison table
**Independent simplification can move the simplified curves farther apart than their originals. Joint simplification cannot exceed the configured pair threshold by construction.** The table reports observed independent increases and separately verifies the joint invariant. With the default exact threshold, that invariant is also a bound by the originals' dDF.

Scroll horizontally for the complete pipeline; the header groups remain visible when scrolling vertically. Arc ratios mean **arc length**, not vertex retention. Empty index-wise fields are intentional. Gray rows are infeasible and excluded from aggregates; orange independent-distance cells exceed original dDF. Exact data (before display rounding) are exported to CSV.
""")
code("""
table_html=pipeline_html(pipeline)
display(HTML(table_html))
(OUT/'pipeline_table.html').write_text('<!doctype html><meta charset="utf-8"><title>Small protein pair pipeline</title><h1>Independent versus synchronous discrete CPS</h1>'+table_html,encoding='utf-8')
""")
md('### The three discrete distances side by side\nThis compact view makes the distance comparison readable without scrolling through all pipeline columns. Independent distance remains meaningful even when the two output counts differ; its index-wise coupling fields remain empty in the full table.')
code("""
distance_view=pipeline[['B_printed','delta1','original_dDF','separate_dDF','CPS_dDF','delta3',
    'separate_further_apart','status']]
display(distance_view.round(3))
display(aggregates.round(3))
print(f'Observed increases among feasible rows: {int(feasible.separate_further_apart.sum())}/{len(feasible)}.')
print('Aggregate denominators exclude every infeasible pair/threshold; all such rows remain visible in the full table and status plot.')
""")
code("""
plot_outcomes(pipeline,FIGURES)
""")
md('### An illustrative full pipeline\nThe following figure selects the largest observed independent distance increase **after** completing the unfiltered experiment. It is an illustration, not the rule for selecting input windows or aggregate rows. Original, independent and joint curves use the same frozen superposition and axes.')
code("""
example_row=plot_case(pairs,pipeline,selected_paths,FIGURES)
display(example_row[['pair','delta1','original_dDF','separate_dDF','CPS_dDF','separate_k_A','separate_k_B','CPS_k','cost_of_joining']].to_frame('value'))
""")
md(r"""
## 7. Planned comparisons and unresolved choices

| Method / extension | Present state | Planned comparison |
|---|---|---|
| Shortest independent, local continuous fidelity | Implemented and run on all selected windows | Current baseline at 2, 4, 8 Å; preserve its tie policy and local criterion |
| Synchronous discrete CPS | Exact solver, floor and reconstruction implemented; all feasible rows verified | Current joint comparison; equal size, diagonal pairing, fixed endpoints |
| General discrete CPS of Fan et al. | Different formulation; not substituted for the synchronous solver | Permit unequal cardinalities and non-diagonal coupling; report max size and leave undefined index-wise fields blank |
| Independent with discrete fidelity | Future metric-matched baseline | Separate the effect of changing fidelity metric from the effect of joining |
| CPS-2F | Explicit `NotImplementedError` interface | Define its two Fréchet constraints, metric versions, vertex domain, endpoints and objective when additional data arrive |

**Choices and limits to retain in future comparisons:**

- **Threshold convention:** exact original dDF is the primary threshold; the paper's ceiling is archived and available via `DELTA3_MODE`. A ceiling admits rounding slack and changes the meaning of the headline bound.
- **Meaning of CPS:** the requested equal-size/index-wise columns lead to a stricter synchronous variant. General Fan CPS and this variant must not share unlabeled optimum counts. Endpoints are fixed here for both methods.
- **Alignment:** US-align is a reproducible structural alignment heuristic, not a proof of optimal registration under dDF. Low-scoring dissimilar windows remain included and visible. We do not replace it with index-wise Kabsch, discard unmatched residues, or tune transforms separately for each method.
- **Window selection:** the rule chooses a single correspondence-centered continuous window per side, not a representative statistical sample of the whole fold. A varies across pairs. Future experiments should predeclare multiple windows/seeds and aggregate by parent protein pair, avoiding pseudo-replication.
- **Original data:** the verified interpretations and historical length discrepancies from notebook 0 still apply. No unresolved region is connected. A longest-fragment experiment excludes other fragments explicitly; whole-protein results require a declared fragment-pairing policy.
- **Interpretation of joining cost:** the requested difference combines joint constraints and different fidelity measures. Add the metric-matched baseline before attributing that difference solely to coupling.
- **Numerics:** double precision; local geometric primitive uses `1e-10` predicates, solver bound comparisons use `1e-9` Å tolerance. The floor is an exact finite-state minimax optimum subject to those numerical comparisons.

### Sources and saved artifacts

- Fan, Filtser, Katz, Wylie and Zhu, [_On the Chain Pair Simplification Problem_](https://arxiv.org/pdf/1409.2457), §7/Table 1: the paper uses the ceiling of original discrete Fréchet distance for the pair threshold, and 4 Å for the two fidelity bounds. The 2/4/8 Å sweep is this notebook's experiment. Their approximately 3.7–3.8 Å Cα spacing motivates 4 Å but does not make continuous and discrete fidelity equivalent.
- [Official US-align distribution](https://seq2fun.dcmb.med.umich.edu/US-align/) and [source repository](https://github.com/pylelab/USalign): structural alignment with inferred correspondence. Download provenance and binary SHA-256 are in `tools/usalign/source.json`.
- `data/prepared/*.residues.csv` and notebook 0's PDB audit retain original residue/model/chain identities. Current `output/small_pairs/alignments/` contains selected-window maps, aligned arrays, transforms and raw tool logs.

**Outputs:** `output/small_pairs/pipeline.csv`, `pipeline_table.html`, `aggregates_feasible_only.csv`, `selected_residues.csv`, `paths.json`, `alignment_audit.csv`, and `figures/`. Rerunning the notebook recomputes every solver row; structural fits are reused only when their input/source signature matches. Set `RECOMPUTE_ALIGNMENT=True` to refit. Notebook 0 and its results are preserved.
""")
nb=nbf.v4.new_notebook(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3 (protein analysis)',language='python',name='python3'),language_info=dict(name='python',version='3.14')))
nbf.validate(nb);nbf.write(nb,root/'1.analyze_independent.ipynb')
print(f'Created {len(cells)} notebook cells.')
