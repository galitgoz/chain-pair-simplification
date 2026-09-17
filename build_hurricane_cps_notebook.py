from pathlib import Path
import nbformat as nbf
root=Path(__file__).resolve().parent
cells=[]
def md(s):cells.append(nbf.v4.new_markdown_cell(s))
def code(s):cells.append(nbf.v4.new_code_cell(s))
md('''# Hurricane tracks: Independent GCS, CPS-2F and CPS-3F

Ordered 2D storm tracks from `data/hurricanes/`. Distances are **Euclidean kilometres in the supplied basin-specific projection**, not great-circle distances. Track order is retained; timestamps are not present, so coupling is geometric and is not synchronization in physical time.

All tracks are loaded. The executed default selects four pairs per basin across the supplied ball-packedness rankings (eight pairs total), at full supplied resolution, with α=1 and 2. This is a selected comparison, not a census of all 2,222 candidate pairs. Set `PAIRS_PER_BASIN=None` to run every supplied pair. Read the [dataset provenance and preprocessing](data/hurricanes/README.md) for historical HURDAT2 versions, projection, filtering, and measured packedness definitions.''')
md('## 1. Load the data')
code('''from pathlib import Path
import json, hashlib, inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, HTML, Markdown
from hurricane_cps import load_hurricanes, organize, run, comparison_plots, track_plots
ROOT=Path.cwd()
DATA=ROOT/'data/hurricanes'
OUT=ROOT/'output/hurricane_cps_expanded'
OUT.mkdir(parents=True,exist_ok=True)
PAIRS_PER_BASIN=16  # None runs all 2,222 supplied candidate pairs
ALPHAS=(1,2)
BUDGETS=dict(max_states=300_000,max_transitions=20_000_000,seconds=120.,max_chain_edges=2_000_000)
tracks, metadata, candidate_pairs=load_hurricanes(DATA)''')
md('## 2. General properties — three cells')
code('''display(pd.DataFrame({'property':['Tracks','Candidate pairs','Dimensions','Coordinate units','Supplied continuous > discrete distances'],
    'value':[len(tracks),len(candidate_pairs),2,'projected km',int(candidate_pairs.supplied_distance_order_violation.sum())]}))
display(metadata.groupby('basin').agg(tracks=('id','size'),vertices_min=('n','min'),vertices_median=('n','median'),vertices_max=('n','max')))''')
code('''display(metadata[['n','length_km','extent_km','c_ball','c_cube']].describe().round(3))
display(candidate_pairs[['curve_A','curve_B','basin','frechet_km','discrete_frechet_km']].head())''')
code('''fig,axes=plt.subplots(1,2,figsize=(10,3))
for basin,g in metadata.groupby('basin'):
    axes[0].hist(g.n,bins=20,alpha=.5,label=basin)
    axes[1].scatter(g.n,g.c_ball,s=8,alpha=.4,label=basin)
axes[0].set(xlabel='Vertices per track',ylabel='Track count')
axes[1].set(xlabel='Vertices per track',ylabel='Supplied measured ball packedness')
axes[0].legend();axes[1].legend();fig.tight_layout();plt.show();plt.close(fig)''')
md('''## 3. Organize tracks into reproducible paired inputs

`tracks[id]` holds an n×2 array; `data['instances'][pair+'@1']` holds A and B and their basin metadata. The suffix `@1` means stride 1: no further decimation. Selection takes evenly spaced ranks of supplied `c_ball_pair` within each basin, breaking ties by track IDs. It is independent of solver outcomes. Supplied packedness values are measurements, not exact constants; they are used for cohort coverage, never substituted into theoretical graph bounds.

The supplied pair metadata contains distance-order inconsistencies. Selected-pair distances are recomputed from the coordinate files. δ₃ is rounded upward to the next kilometre from the recomputed original discrete Fréchet distance, ensuring the original pair is feasible. δ₁=δ₂=α times the average of the two tracks' mean edge lengths. The pair table records every threshold.''')
code('''data,selected=organize(tracks,candidate_pairs,PAIRS_PER_BASIN)
selected.to_csv(OUT/'selected_pairs.csv',index=False)
display(selected[['curve_A','name_A','year_A','curve_B','name_B','year_B','basin','n_A','n_B','c_ball_pair']])
parameters=[]
from curve_algorithms import discrete_frechet
for item in data['instances'].values():
    A,B=item['A'],item['B_xyz']
    mean_edge=.5*(np.linalg.norm(np.diff(A,axis=0),axis=1).mean()+np.linalg.norm(np.diff(B,axis=0),axis=1).mean())
    for alpha in ALPHAS:
        parameters.append(dict(pair=item['pair'],basin=item['basin'],alpha=alpha,nA=len(A),nB=len(B),
            delta1=alpha*mean_edge,delta2=alpha*mean_edge,delta3=float(np.ceil(discrete_frechet(A,B))),stride=1))
parameters=pd.DataFrame(parameters)
parameters.to_csv(OUT/'parameters.csv',index=False)
display(parameters.round(3))''')
md('''## 4. Algorithms and validation

All methods retain original vertices and anchor both endpoints; the objective is k=max(|A′|,|B′|).

| Method | A→A′, B→B′ fidelity | A′↔B′ coupling |
|---|---|---|
| Independent GCS | Global continuous Fréchet ≤ δ₁/δ₂ | Measured, unconstrained |
| CPS-2F | Continuous Fréchet ≤ δ₁/δ₂ | Discrete Fréchet ≤ δ₃ |
| CPS-3F | Discrete Fréchet ≤ δ₁/δ₂ | Discrete Fréchet ≤ δ₃ |

The same implementation as the protein notebook is reused: [van de Kerkhof et al., Algorithm 1](papers/global%20curve%20simplification.pdf), [Fan et al., Algorithm 1](papers/chain_pair_simplification.pdf), and [Galit Gozoltzani thesis, Chapters 3–5](papers/M_Sc__Thesis___Galit_Gozoltzani.pdf). CPS-2F uses auxiliary sphere/edge intersections as owner positions; selected simplified vertices remain original vertices. The local thesis is an annotated draft; its discretization claim is a source assumption, not established by these experiments.

Exact lower-bound certificates can establish joint optimality without constructing the product graph. Otherwise the sparse configuration DP runs with explicit resource limits. A resource limit is reported as missing, never infeasible or zero. The time limit is soft and excludes any guarantee of interruption during preprocessing or certificate search.''')
code('''from curve_algorithms import shortest_independent, continuous_decision
from cps_paper_algorithms import solve_cps
from verify_cps_papers import candidates
import itertools
rng=np.random.default_rng(20260915)
checks=0
for trial in range(4):
    A=np.cumsum(rng.normal(size=(4,2)),axis=0);B=A+rng.normal(scale=.2,size=(4,2))
    d1=d2=.8;d3=max(np.linalg.norm(A[0]-B[0]),np.linalg.norm(A[-1]-B[-1]))+.8
    for method in ['CPS-2F','CPS-3F']:
        ca=list(candidates(A,d1,method));cb=list(candidates(B,d2,method))
        optimum=min((max(len(a),len(b)) for a in ca for b in cb if discrete_frechet(A[a],B[b])<=d3+1e-9),default=None)
        actual=solve_cps(A,B,d1,d2,d3,method,certificate=False,seconds=30)
        assert actual.get('k')==optimum
        if method=='CPS-2F':assert len(shortest_independent(A,d1))==min(map(len,ca))
        checks+=1
print(f'{checks} exhaustive 2D joint-optimization checks passed.')''')
md('## 5. Run the three methods')
code('''results,paths,distance_audit=run(data,selected,ALPHAS,BUDGETS,OUT)
display(results.groupby(['method','status']).size().rename('cases').reset_index())
display(distance_audit.round(4))''')
md('''## 6. Results and geometric accuracy

The fidelity columns measure continuous Fréchet error in km for both methods, providing a common geometric comparison. They do not measure topology, landfall preservation, forecast error, or time alignment. Compression and coupling are reported separately; smaller k does not imply smaller achieved coupling. In these open storm tracks, no protein contact metric is applied.''')
code('''columns=['pair','alpha','method','nA','nB','kA','kB','k','delta1','delta2','delta3',
    'removed_pct','budget_ratio','fidelity_A_upper','fidelity_B_upper','vertices_over_GCS','seconds','status','reason']
display(HTML(results[columns].to_html(index=False,float_format=lambda x:f'{x:.3f}')))
solved=results[results.status.str.startswith('optimal')]
assert (solved[solved.method.ne('GCS')].budget_ratio<=1+1e-8).all()
assert (solved.fidelity_ratio<=1+1e-5).all()
matched=results.pivot(index='case',columns='method',values='k').dropna()
assert (matched.GCS<=matched['CPS-2F']).all() and (matched['CPS-2F']<=matched['CPS-3F']).all()
matched['CPS2F_minus_CPS3F']=matched['CPS-2F']-matched['CPS-3F']
matched.to_csv(OUT/'matched_vertex_comparison.csv')
display(matched)''')
md('''## 7. Coupling budget and compression

Grouped bars identify each pair; parameters.csv lists the exact thresholds. NA marks unsolved cases. Boxplots use only pair–threshold cases solved by all three methods, with the same cohort for each method. Repeated thresholds and storms are not independent statistical observations. The exported cohort identifies inclusions.''')
code("comparison_plots(results,OUT/'figures')")
md('''## 8. Graph size and auxiliary utilization

Graph size here is **valid vertices before reachability pruning**, not edges or only discovered states. CPS-3F has no auxiliary points and one bound m²n². CPS-2F has measured-extension bound m*m n*n (where m*=|A*| and n*=|B*|), and the explicit worst-case bound [m+2m(m−1)]m[n+2n(n−1)]n, of order O(m³n³). Legends show exact formulas. Auxiliary utilization is 100 × added points / [2n(n−1)] separately for each curve. All original vertices are excluded from the auxiliary numerator and denominator. Work counters remain separate from full valid graph counts.''')
code('''from cps_graph_sizes import graph_size_table,graph_size_plots,auxiliary_utilization_plots
graph_sizes=graph_size_table(data,results,OUT)
(OUT/'graph_sizes.html').write_text((OUT/'graph_sizes.html').read_text(encoding='utf-8').replace('All distances are in angstroms.','All distances are in projected kilometres.'),encoding='utf-8')
graph_size_plots(graph_sizes,OUT/'figures')
utilization=auxiliary_utilization_plots(graph_sizes,OUT/'figures')
display(HTML(utilization.to_html(index=False,float_format=lambda x:f'{x:.3f}')))
display(Markdown('[Graph table](output/hurricane_cps/graph_sizes.csv) · [Utilization table](output/hurricane_cps/auxiliary_utilization.csv)'))''')
md('''## 9. Track geometry and coupling for every selected pair

At α=1, each pair has three method panels sharing the same coordinates and scale. Faint lines show inputs, bold lines show simplified tracks, and dotted lines show the discrete coupling, including repeated indices. Stars mark starts. Basin-projected axes are kilometres; no map background or geographical interpretation is inferred. Change the final argument to inspect α=2.''')
code("track_plots(data,results,paths,OUT/'figures',alpha=ALPHAS[0])")
md('''## 10. Interpretation and reproducibility

- Compare Independent with CPS-2F to measure the vertex cost of adding the coupling constraint under continuous fidelity. CPS-2F versus CPS-3F changes the fidelity metric; both enforce the same discrete coupling budget.
- Bounds are upper limits, not predictions. Auxiliary utilization and valid graph size depend on geometry and thresholds. These experiments do not establish a c-packed runtime theorem.
- Input curves already reflect the dataset's filtering, projection and removal of repeated consecutive coordinates. No missing observations or physical times are reconstructed.
- This default cohort is selected for coverage across supplied packedness, not a population estimate. Resource-limited results remain visible; boxplot exclusions can change the represented cohort.

CSV tables, paths, and publication-ready PDF/SVG/PNG figures are saved under `output/hurricane_cps/`. The manifest records inputs, software sources, selection and limits.''')
code('''source_files=[ROOT/name for name in ['hurricane_cps.py','curve_algorithms.py','cps_paper_algorithms.py','cps_graph_sizes.py']]
input_files=sorted(DATA.glob('*.txt'))+[DATA/'pairs.csv',DATA/'packedness.csv',DATA/'README.md']
manifest=dict(pairs_per_basin=PAIRS_PER_BASIN,alphas=list(ALPHAS),stride=1,units='projected km',
    loaded_tracks=len(tracks),candidate_pairs=len(candidate_pairs),selected_pairs=len(selected),
    budgets=BUDGETS,algorithm_checks=checks,rows=len(results),status_counts=results.status.value_counts().to_dict(),
    sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files+input_files})
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
display({k:v for k,v in manifest.items() if k!='sha256'})''')
for cell in cells:
    cell.source=cell.source.replace('output/hurricane_cps/','output/hurricane_cps_expanded/')
    cell.source=cell.source.replace('four pairs per basin','16 pairs per basin').replace('(eight pairs total)','(32 pairs total)')
nb=nbf.v4.new_notebook(cells=cells)
nb.metadata.update(kernelspec=dict(display_name='Protein analysis',language='python',name='protein-analysis'),language_info=dict(name='python'))
nbf.write(nb,root/'1.analyze_CPS.ipynb')
print(f'Created {len(cells)} cells')
