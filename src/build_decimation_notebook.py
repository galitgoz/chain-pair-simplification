from pathlib import Path
import nbformat as nbf
root=Path(__file__).resolve().parents[1];cells=[]
def md(s):cells.append(nbf.v4.new_markdown_cell(s.strip()))
def code(s):cells.append(nbf.v4.new_code_cell(s.strip()))

md(r"""
# Protein pairs: lockstep decimation → separate versus joint simplification

**Main comparison:** `sep_max_link` versus `cps_max_link` at the **same recomputed local continuous Fréchet fidelity threshold**. The joint method requires equal output counts and an index-wise link bound; independent outputs may lack such a coupling. We report every unique paper pair at every $(w,\alpha)$: **14 × 5 × 3 = 210 rows**, with full CSV and rendered table.

**Scope choices are explicit:** both methods use the previously specified local per-shortcut continuous fidelity criterion. “Joint CPS” below is a shortest path in the product of those two shortcut DAGs, constrained by retained-pair distances. It is not Fan et al.'s discrete-fidelity CPS-3F algorithm. Discrete Fréchet is used for the requested distance reporting and $\delta_3=\lceil d_{dF}(A_{dec},B_{dec})\rceil$.

**Data qualification is distinct from decimation.** Several full chains have unequal lengths or unresolved gaps. To obtain a continuous, equal-length input while preserving original row-index correspondence, this notebook explicitly selects the longest original-row interval continuous in both chains (earliest tie). It does not use shape, curvature, or simplification results to select it. All excluded residues and full versus qualified lengths are reported. **Full deposited-chain endpoints are not preserved when qualification excludes them.** The endpoints of the qualified input are retained in both chains at every stride. These results are therefore qualified-interval experiments, not undisclosed whole-chain comparisons.

All computations use **Å**. **Headline experiments: $w=1,2,4$ only. $w=8,16$ are RUNTIME ONLY**, always labelled; $w=1$ is the undecimated control in each comparison table and sweep. The ~20-residue PCA instance is display-only and supplies no table measurements.

## 1. Load files and provenance
""")
code("""
from pathlib import Path
import json, time, io, contextlib, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, HTML
from numba import njit
from protein_data import read_raw, paper_experiments
from decimation_data import prepare_decimation_pairs
from curve_algorithms import segment_frechet_decision, discrete_frechet

ROOT=Path.cwd(); assert (ROOT/'data/raw protein backbones').is_dir()
OUT=ROOT/'output/decimation';FIGURES=OUT/'figures'
OUT.mkdir(parents=True,exist_ok=True);FIGURES.mkdir(parents=True,exist_ok=True)
raw={p.stem:read_raw(p) for p in sorted((ROOT/'data/raw protein backbones').glob('*.raw'))}
sources=pd.DataFrame(json.loads((ROOT/'data/pdb_cache/sources.json').read_text()))
WS=(1,2,4,8,16);ALPHAS=(1,2,4)
RECOMPUTE_ALIGNMENT=False
plt.rcParams.update({'font.size':10,'figure.dpi':110})
""")
md('## 2. General properties — three brief cells\nCounts and source metadata only here. All geometry is measured after structural superposition.')
code("""
display(pd.DataFrame([dict(file=key,n=len(P),dimensions=P.shape[1],finite=bool(np.isfinite(P).all())) for key,P in raw.items()]))
print('The eight local raw files each contain 325 C-alpha coordinates; they have no embedded residue or gap metadata.')
""")
code("""
display(sources[['pdb_id','status','retrieved_utc']])
print('Current PDB/mmCIF records augment the local files; source URLs and SHA-256 hashes are cached.')
print('Literal 107J was unavailable (404); literal 4CEA exists and differs from 4ECA.')
""")
code("""
paper_rows=paper_experiments()
display(paper_rows.groupby('table').agg(rows=('B_printed','size')))
print('Tables 1 and 2 repeat seven pairs; Table 3 adds seven: 14 unique pairs for the sweep.')
print('Stride controls: 1 (undecimated), 2 and 4 (headline), 8 and 16 (RUNTIME ONLY). Alpha: 1, 2, 4.')
""")
md(r"""
## 3. Verify identifiers and organize usable pair inputs

Re-read the cached structural records and match every raw Cα coordinate to its deposited atom. Preserve original residue numbers, insertion codes, author versus label chains, model, alternate location and row index. Select the lowest deposited model. Raw alternate choices are retained; for chains without raw files choose the highest-occupancy Cα (ties: blank, A, lexical). Nonpolymer atoms are excluded. Unresolved residues, omitted raw residues and unverifiable/long peptide connections split fragments; they are never connected.

Printed `107j.a` is interpreted as **1O7J author A** by exact raw-coordinate matching and the unavailable literal ID. `4cea.b/d` are interpreted as **4ECA B/D**: literal 4CEA is HIV integrase, has only 150 Cα atoms in B and no protein D. These are supported interpretations, not confirmed author errata. `1toh` means author A; `1hfj.c` is author C but label B. Reported and verified lengths remain separate; no data are trimmed to match reported counts.

For each pair, form the maximal row-index runs that are continuous in **both** chains, choose the longest, and apply that **same interval** to A and B. This qualification is recorded before the user-specified decimation stage. The reference A interval can vary across pairs. There is no 30-point preselection, structural-alignment-based pruning, or resampling.

US-align performs sequence-independent, rigid structural alignment of B to A on the qualified inputs, without reordering or removing their vertices. The fitted transform is frozen before any distance/edge measurement and before all strides. It is not an index-wise Kabsch fit, and it is not refitted after decimation. Structural alignment is heuristic and does not guarantee the minimum discrete Fréchet distance.
""")
code("""
with contextlib.redirect_stdout(io.StringIO()) as preparation_log:
    pairs,audit,qualification=prepare_decimation_pairs(ROOT,recompute_alignment=RECOMPUTE_ALIGNMENT)
display(audit[['printed','verified','reported_length','polymer_sequence_length','deposited_ca_length','prepared_ca_length','author_chain','label_chain','model','fragments']])
assert audit.selected_zero_occupancy_ca.sum()==0
print('All raw coordinate matches:', audit.raw_max_coordinate_error_A.dropna().max(), 'Å maximum error.')
""")
code("""
display(qualification[['B_printed','B_verified','nA_full','nB_full','source_index_start','source_index_end',
    'nA_before','nB_before','excluded_A','excluded_B','full_chain_endpoints_preserved','regime']])
print('All qualification exclusions, including original residue IDs: output/decimation/qualification_excluded_residues.csv')
for p in pairs.values():
    assert len(p['A'])==len(p['B'])
    for side in ('A','B'):
        assert np.all(np.diff(p[f'map_{side}'].label_seq_id)==1)
        assert np.all(np.diff(p[f'map_{side}'].curve_index)==1)
first=next(iter(pairs.values()))
display(first['map_A'][['paired_input_index','curve_index','author_seq_id','insertion_code','residue_name','altloc']].head(5))
print('pairs[id] contains aligned A/B arrays, maps, qualified-source bounds, regime, and transform provenance.')
""")
md(r"""
## 4. Decimation — exactly the supplied rule, in lockstep
The function below is the requested implementation verbatim. Select the index set **once** by applying it to `np.arange(n)`, then index A and B with that same set. Assertions compare both arrays to the supplied function. `w=1` retains every vertex; both input endpoints survive every stride. No Douglas–Peucker, curvature, feature selection, or extra gap-boundary insertion is performed.
""")
alg=(root/'src/decimation_algorithms.py').read_text()
code(alg[alg.index('def decimate'):alg.index('@njit(cache=True)\ndef shortcut_dag')])
md(r"""
### Thresholds and edge statistics, recomputed on each decimated pair
$$e=\tfrac12\left(\mathrm{meanEdge}(A_{dec})+\mathrm{meanEdge}(B_{dec})\right),\quad
\delta_1=\delta_2=\alpha e,\quad\delta_3=\lceil d_{dF}(A_{dec},B_{dec})\rceil.$$
No threshold is copied from $w=1$. `edge_min_after` and `edge_max_after` are the min/max over **both** chains; separate mean/last-edge fields are also saved. `edge_growth_factor = edge_mean_after / edge_mean_before`. Length ratios are **arc length of the simplification divided by arc length of its decimated input**, not a vertex-count ratio. All `dF_*` report discrete Fréchet distance despite the short column name.

## 5. The simplification algorithms and prechecks

Both independent and joint methods use the same local continuous fidelity: edge $i\to j$ is allowed iff $d_F(P[i..j],\overline{P_iP_j})\leq\delta_1$. The independent optimum is a shortest shortcut-DAG path. Joint vertices are retained pairs $(i,j)$; a transition requires a valid shortcut in **each** chain and strict progress in **both** indices. Joint path length is the shared output count $k$, minimized subject to each retained-pair distance being at most $\delta_3$. This permits different source indices at corresponding **output** positions; only decimation uses a common source index set.

**The requested `floor`** is the bottleneck over all strictly increasing pairings from `(0,0)` to `(n-1,m-1)`, ignoring fidelity. A direct jump is legal in that definition, so it equals the maximum endpoint distance. It does **not** equal original dDF or a fixed index-wise maximum. We compute it with a prefix-minimum DP and verify the endpoint identity.

This floor is necessary, not sufficient once shortcut fidelity is imposed. The precheck additionally computes a clearly named `local_fidelity_floor` over valid shortcut-pair paths. A failed endpoint precheck is `infeasible_endpoint_floor`; a failed fidelity precheck is `infeasible_local_fidelity`. `feasible` retains your exact definition `delta3 >= floor`, while `solver_feasible` reports full feasibility. The cardinality solver is called only after both prechecks. If it then fails, the notebook raises a solver error rather than labelling the failure infeasible.

With the prescribed $\delta_3=\lceil d_{dF}(A_{dec},B_{dec})\rceil$, the endpoint-only precheck must pass, because discrete Fréchet includes both endpoints. It is retained to support lower user-specified bounds and tested on infeasible examples. Double-precision predicates use the existing segment tolerance `1e-10` and bound checks `1e-9` Å.
""")
code(alg[alg.index('@njit(cache=True)\ndef shortcut_dag'):].replace('cache=True','cache=False'))
md('### Validate decimation, optimality and infeasibility\nThe tests exercise exact stride sets and appended endpoints, lockstep application, exhaustive local-path and joint optima, and both infeasibility statuses.')
code("""
from verify_decimation import run_checks
checks=run_checks()
# Exercise the notebook-defined decimator as well.
test_A=np.arange(99).reshape(33,3);test_B=-test_A
for w in WS:
    a,b,idx=lockstep_decimate(test_A,test_B,w)
    assert np.array_equal(a,test_A[idx]) and np.array_equal(b,test_B[idx])
assert decimate(test_A,1) is test_A
# Compile the notebook-defined kernels before any reported runtime is measured.
warm=np.array([[0.,0,0],[1.,0,0],[2.,0,0]])
warm_edges=shortcut_dag(warm,1.0)
assert len(shortest_path(warm_edges))==2
assert joint_local(warm,warm,warm_edges,warm_edges,0.)['k']==2
""")
md(r"""
## 6. Run all 210 before/after pipelines
Independent and joint algorithms receive the **same decimated originals** and thresholds. If independent counts differ, link/check fields are empty and the coupling status states `no coupling, |A'| != |B'|`. The independently computed discrete Fréchet distance is still reported, without being substituted for an undefined diagonal coupling.

Joint output verifies the local shortcut edges, endpoints, strict order, equal sizes, maximum link and discrete Fréchet bound. The bound is a **solver verification**, not an empirical result. The joint method minimizes count under a hard link limit; it does not promise a smaller maximum link than every already-feasible independent result. `cost_of_joining = cps_k - max(sep_kA, sep_kB)` uses identical fidelity criteria here. Compilation is warmed up before timing; `graph_seconds` records shared shortcut construction separately, and `cps_seconds` includes its two prechecks plus the joint count solver.
""")
pipeline=(root/'src/decimation_pipeline.py').read_text()
code(pipeline[pipeline.index('def run_sweep'):pipeline.index("if __name__=='__main__':")])
code("""
started=time.perf_counter()
results,paths=run_sweep(pairs,OUT,ws=WS,alphas=ALPHAS,verbose=False)
assert len(results)==len(pairs)*len(WS)*len(ALPHAS)==210
assert not results.duplicated(['pair','w','alpha']).any()
for pairid,g in results.groupby('pair'):
    assert set(g.w)==set(WS) and (g.w==1).sum()==3
assert np.allclose(results.delta1,results.alpha*results.edge_mean_after)
assert np.array_equal(results.delta3,np.ceil(results.dF_after))
assert np.array_equal(results.feasible,results.delta3>=results.floor)
control=results[results.w==1]
assert (control.nA_before==control.nA_after).all() and (control.keep_frac==1).all()
assert np.allclose(control.dF_change,0) and np.allclose(control.edge_growth_factor,1)
assert results.loc[~results.sep_equal,['sep_max_link','sep_ok3']].isna().all().all()
successful=results[results.solver_feasible]
assert successful.cps_verified.all() and (successful.cps_max_link<=successful.delta3+1e-9).all()
print(f'{len(results)} rows; {len(successful)} joint solutions verified; {time.perf_counter()-started:.1f} seconds including graphs and prechecks.')
""")
md(r"""
## 7. Complete before/after report
The rendered table includes every pair, stride and alpha. It scrolls horizontally and vertically; the full unrounded CSV is `output/decimation/pipeline.csv`. Metadata and per-stage residue mappings make the input-qualification, raw-paired, decimated, separate and joint stages traceable. `w=1` control rows are highlighted. Every `w=8/16` row is labelled **RUNTIME ONLY**.
""")
code("""
from decimation_figures import full_table,F1,F2,F3,F4,headline_links
rendered=full_table(results)
display(HTML(rendered))
intro='<h1>Lockstep decimation: complete pair pipeline</h1><p>Headline: w≤4 only. w=8/16: RUNTIME ONLY. Before means qualified paired input; full-chain exclusions are in input_qualification.csv.</p>'
(OUT/'pipeline_table.html').write_text('<!doctype html><meta charset="utf-8"><title>Lockstep decimation pipeline</title>'+intro+rendered,encoding='utf-8')
""")
md('### Sampling diagnostics and headline coupling counts\nThese summary tables retain every stride, including its undecimated control, and label runtime-only rows. Alpha-independent sampling quantities are counted once per `(pair,w)`; coupling counts remain per alpha.')
code("""
sampling=results.drop_duplicates(['pair','w']).groupby(['w','stride_role']).agg(
    pairs=('pair','size'),mean_edge_Angstrom=('edge_mean_after','mean'),smallest_edge_Angstrom=('edge_min_after','min'),
    mean_growth_factor=('edge_growth_factor','mean'),stub_cases=('appended_stub','sum'),
    coarse_trace_cases=('coarse_trace','sum'),mean_dF_change_Angstrom=('dF_change','mean')).reset_index()
display(sampling.round(3))
sampling.to_csv(OUT/'sampling_summary_with_control.csv',index=False)
summary=results.assign(separate_violates=lambda t:t.sep_equal & (t.sep_ok3==False)).groupby(['w','stride_role','alpha']).agg(
    total_pairs=('pair','size'),separate_coupled=('sep_equal','sum'),separate_violates=('separate_violates','sum'),
    joint_feasible=('solver_feasible','sum'),mean_joining_cost=('cost_of_joining','mean'),mean_cps_seconds=('cps_seconds','mean')).reset_index()
display(summary.round(3))
summary.to_csv(OUT/'coupling_summary_with_control.csv',index=False)
""")
md(r"""
## 8. Figures

### F1 — vertices at each stage
One small panel per pair, with a line style for each stride. Separate and joint are distinguished by **colour and marker**. There is a complete figure for each alpha, avoiding an unreported choice of fidelity. “Raw*” means the qualified paired input; the full source counts and exclusions are in the qualification table. Y is logarithmic. RUNTIME ONLY strides remain labelled in every legend.
""")
code('F1(results,FIGURES)')
md('### F2 — when decimation changes the distance being measured\nBoth axes are logarithmic. Colour encodes an operational regime from the **fixed pre-decimation structural alignment**: TM-score ≥0.5 similar, 0.3–0.5 intermediate, <0.3 dissimilar. These are declared plotting bins, not a new biological classification. The shaded region is RUNTIME ONLY.\n\nA flat distance can be dominated by retained endpoints; it does not prove that helical detail survived. Read F2 together with edge growth and the coarse-trace flags.')
code('F2(results,FIGURES)')
md('### F3 — usable separate coupling at every `(w, alpha)`\nColour **and** symbol distinguish no coupling, a violating coupling, and a coupling meeting δ₃. There is no proxy for unequal output sizes.')
code('F3(results,FIGURES)')
md('### Headline comparison — maximal index-wise links at the same δ₁\nOnly w≤4 appears here, including the w=1 control. Unequal separate sizes have no comparable maximum link, so they remain explicit in F3 and the full table rather than being plotted as zero.')
code('headline_links(results,FIGURES)')
md('### F4 — a separate display-only instance\nTake 20 consecutive residues from the middle of the first dissimilar pair (or the first pair if none is dissimilar), at w=1 and α=2. Reuse the frozen structural transform. Compute PCA jointly for both 3D windows and draw A alone, B alone, and the joint result with index-wise links. The caption reports retained variance. Solving stays in 3D; **no measurement from this display window is added to any results table**.')
code('F4(pairs,FIGURES)')
md(r"""
## 9. Caveats and planned extensions

- **Stub edges:** when w does not divide n−1, the appended endpoint creates a final step shorter than w in source-index space. Its geometric chord may be short and is not averaged away: inspect `edge_min_after`, `last_stride`, `edge_last_A_after`, `edge_last_B_after`, and `appended_stub`.
- **Folded chords:** mean edge grows sublinearly with stride because a chord is shorter than the intervening arc in a folded chain. `edge_growth_factor` and `growth_per_stride` quantify this. Once mean edge exceeds roughly **10 Å**, interpret the input as a coarse trace that has lost helical resolution, rather than a resolved protein backbone. `coarse_trace` flags that warning; it is not a sharp physical theorem.
- **Headline scope:** only **w≤4**, including w=1 control, supports headline comparisons. **w=8/16 are RUNTIME ONLY** in tables, legends, figures and summary outputs. Coarse-trace warnings still apply to any w≤4 row exceeding 10 Å. Runtime summaries include w=1 as a reference.
- **Threshold rescaling is part of the experiment:** alpha is fixed but delta1 changes with measured edge length at every stride. These are bond-scaled experiments, not experiments at one constant Å error. Delta3 also changes after decimation; neither threshold is carried over from the control.
- **Distance versus coupling:** ceil(dDF) allows rounding slack, so even a verified joint dF may slightly exceed unrounded original dF. The actual guarantee is `cps_max_link <= delta3`, and consequently `cps_dF <= delta3`. A smaller dF after decimation can reflect lost disagreement, not improved accuracy. A feasible independent coupling need not be worsened/improved monotonically by the count-optimal joint solution.
- **Input qualification limits the scope:** all 14 named pairs appear, but unequal lengths and gaps require the explicitly reported common continuous intervals. Full-chain endpoints excluded at this stage are not silently claimed to survive; only the qualified-input endpoints are retained at all strides. No claim reproduces the paper's whole-chain preprocessing or historical counts. Another preparation policy requires an explicit new dataset version.
- **Alignment and ties:** US-align provides a fixed structural superposition, not the exact dF-optimal transform. Index correspondence remains original row order; structural correspondence is used only to obtain the rigid transform. The DAG algorithms resolve ties by the first increasing predecessor, and every retained index is exported.
- **Feasibility versus software failure:** endpoint and local-fidelity prechecks are separate from the optimal-count solver. Infeasible rows retain their reason and have blank joint measurements; they are excluded from joint aggregates. A contradiction after prechecks raises an error. Bound satisfaction is reported as solver verification.

**Future comparisons:** retain this local-continuous baseline; add general discrete-fidelity CPS-3F as a separately named method (different problem); define CPS-2F's constraints and objective before implementing it; evaluate additional predetermined intervals and larger datasets without changing the exact decimator. Do not compare unlike fidelity rules as a pure cost of joining.

### Reproducibility and sources

The notebook embeds the decimator, local/product-DAG solvers and complete pipeline. Companion modules are independently tested. Native US-align logs, transforms, and hashes are cached under `output/decimation/inputs/` and `tools/usalign/source.json`. PDB identities are verified directly against `data/pdb_cache/`; original files remain intact.

Fan et al., [_On the Chain Pair Simplification Problem_](https://arxiv.org/pdf/1409.2457), §7 supplies the protein pairs and ceiling threshold convention; this notebook's local continuous joint solver and stride/alpha sweeps are separately specified experiments. Structural alignment uses the [official US-align tool](https://seq2fun.dcmb.med.umich.edu/US-align/) ([source](https://github.com/pylelab/USalign)). Notebook 0 remains unchanged.

Outputs: `pipeline.csv`, `pipeline_table.html`, `input_qualification.csv`, `qualification_excluded_residues.csv`, `stage_residue_mapping.csv`, `paths.json`, control-containing summaries, and all F1–F4 PNGs under `output/decimation/`. F4 is display-only. Run from the project root with the existing Python environment. Every notebook run recomputes solver results; only unchanged structural alignments are cached.
""")
nb=nbf.v4.new_notebook(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3 (protein analysis)',language='python',name='python3'),language_info=dict(name='python',version='3.14')))
nbf.validate(nb);nbf.write(nb,root/'notebooks/legacy/1.analyze_independent.ipynb');print(f'Built {len(cells)} cells')
