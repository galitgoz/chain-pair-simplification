"""Execute only new graph reporting; preserve existing solver results and cells."""

# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path[:0] = [str(_Path(__file__).resolve().parents[2] / 'src'), str(_Path(__file__).resolve().parents[2])]

from pathlib import Path
import copy, json, sys
import nbformat as nbf
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

ROOT = Path(__file__).resolve().parents[2]
heading = r'''### 8.3 Configuration graph G: measured size versus theory

Here graph size means **the number of vertices (valid configurations), $|V(G)|$**, before anchored-endpoint reachability pruning. All pairs are counted, including certificate and resource-limited cases. No joint optimization is rerun.

Write $m=|A|$, $n=|B|$, $m^*=|A^*|$, and $n^*=|B^*|$. CPS-3F has **two bars** (valid vertices and one bound); CPS-2F has **three bars**. Each panel has its own legend with the calculation:

- **Exact valid vertices:** $\sum_{p,q:\|a_p-b_q\|\leq\delta_3} c_A(p)c_B(q)$, where $c_A(p)$ counts owner positions of $A^*$ within $\delta_1$ of $a_p$, and similarly for $B$. This factorization counts the graph without materializing its Cartesian product. Counts use the solver's numerical tolerance.
- **Cartesian upper bound using measured extended sizes (CPS-2F only):** $m^*m n^*n$.
- **Worst-case theoretical upper bound:** CPS-2F uses $[m+2m(m-1)]m[n+2n(n-1)]n=O(m^3n^3)$; CPS-3F uses $m^2n^2$. CPS-3F has no auxiliary points and uses original vertices as owner positions, so its coincident bounds are shown by a single bar. These are bounds, not predictions of equality.

Source: [Galit Gozoltzani thesis, §5.1](papers/M_Sc__Thesis___Galit_Gozoltzani.pdf), together with the four-index discrete configuration graph in [Fan et al., Algorithm 1](papers/chain_pair_simplification.pdf). The log scale makes the large gap to the worst-case bound visible. Independent GCS uses a different graph and is omitted from this joint-graph comparison.

The table also records **discovered states** and **examined transitions** from the existing solver run. These measure solver work, not full graph size: certificate runs build no joint graph, and missing counters were not recorded for interrupted runs. Full $|E(G)|$ is not measured; the thesis gives the asymptotic edge bound $O(|V|^2)$. Auxiliary counts exclude original vertices. Every row specifies the pair, stride, and all three distance budgets in ångströms.
'''
source = '''from cps_graph_sizes import graph_size_table, graph_size_plots, check_graph_counts, auxiliary_utilization_plots
graph_count_checks = check_graph_counts()
graph_sizes = graph_size_table(data, results, OUT)
print(f'Graph-count verification: {graph_count_checks} brute-force checks passed.')
graph_size_plots(graph_sizes, OUT/'figures')
display(Markdown('#### Auxiliary-point utilization percentage\\n\\n'
    'For each CPS-2F pair, the bars show **100 × added auxiliary points / theoretical auxiliary bound**, '
    'separately for A and B. The bounds are 2|A|(|A|−1) and 2|B|(|B|−1); original vertices are excluded. '
    'CPS-3F has no auxiliary points and is omitted. Zero denominators are reported as N/A.'))
auxiliary_utilization = auxiliary_utilization_plots(graph_sizes, OUT/'figures')
display(HTML(auxiliary_utilization.to_html(index=False, float_format=lambda x:f'{x:.3f}')))
display(Markdown(f'[Auxiliary utilization and exact budgets]({OUT.relative_to(ROOT).as_posix()}/auxiliary_utilization.csv)'))
display(HTML(graph_sizes.to_html(index=False, float_format=lambda x:f'{x:.5g}')))
display(Markdown(f'[Graph counts and exact budgets]({OUT.relative_to(ROOT).as_posix()}/graph_sizes.csv)'))
graph_manifest = dict(ws=list(WS), alphas=list(ALPHAS), checks=graph_count_checks,
    rows=len(graph_sizes), definition='full valid vertices before reachability pruning',
    sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in ['src/cps_graph_sizes.py','src/cps_paper_algorithms.py','src/curve_algorithms.py']},
    comparison_sha256=hashlib.sha256((OUT/'comparison.csv').read_bytes()).hexdigest())
(OUT/'graph_size_manifest.json').write_text(json.dumps(graph_manifest,indent=2))
'''

if __name__ == '__main__':
    path = ROOT/'notebooks/legacy/0.analyze_CPS.ipynb'
    nb = nbf.read(path, as_version=4)
    backup = ROOT/'output/cps_papers_w8/0.analyze_CPS_before_graph_sizes.ipynb'
    if not backup.exists(): backup.write_bytes(path.read_bytes())
    loading = next(c for c in nb.cells if c.cell_type=='code' and 'data = load_inputs(' in c.source)
    mini = nbf.v4.new_notebook(cells=[copy.deepcopy(loading),
        nbf.v4.new_code_cell("results = pd.read_csv(OUT/'comparison.csv')"),
        nbf.v4.new_code_cell(source)])
    kernel_root = ROOT/'.analysis_kernels'
    km = KernelManager(kernel_name='protein-analysis', kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(kernel_root)]))
    NotebookClient(mini, km=km, timeout=600, resources={'metadata':{'path':str(ROOT)}}).execute()
    newcell = mini.cells[-1]
    newcell.execution_count = max((c.execution_count or 0 for c in nb.cells if c.cell_type=='code'), default=0)+1
    newcell.metadata['reporting_execution'] = 'Executed in fresh kernel using current inputs and saved comparison.csv'
    positions = [i for i,c in enumerate(nb.cells) if c.source.startswith('### 8.3 Configuration graph')]
    if positions:
        nb.cells[positions[0]:positions[0]+2] = [nbf.v4.new_markdown_cell(heading), newcell]
    else:
        i = next(i for i,c in enumerate(nb.cells) if c.source.startswith('## 9.'))
        nb.cells[i:i] = [nbf.v4.new_markdown_cell(heading), newcell]
    nbf.validate(nb); nbf.write(nb,path)
    print('Added executed graph-size plots and table; preserved all prior outputs.')
