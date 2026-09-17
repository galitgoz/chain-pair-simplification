"""Configuration-graph vertex counts without materializing the product graph."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cps_paper_algorithms import extended_curve, EPS


def valid_vertex_count(A, B, XA, XB, d1, d2, d3):
    """Count all valid (owner A, retained A, owner B, retained B) states.

    This counts the full constraint-valid graph before endpoint reachability
    pruning, not just the states discovered by the anchored solver.
    """
    ca = (np.linalg.norm(XA[:, None] - A[None, :], axis=2) <= d1 + EPS).sum(axis=0)
    cb = (np.linalg.norm(XB[:, None] - B[None, :], axis=2) <= d2 + EPS).sum(axis=0)
    allowed = np.linalg.norm(A[:, None] - B[None, :], axis=2) <= d3 + EPS
    return sum(int(ca[p]) * int(cb[q]) for p, q in zip(*np.nonzero(allowed)))


def graph_size_table(data, results, out):
    out = Path(out)
    rows = []
    for r in results.loc[results.method.isin(['CPS-2F', 'CPS-3F'])].itertuples():
        item = data['instances'][f'{r.pair}@{r.w}']
        A, B = item['A'], item['B_xyz']
        m, n = len(A), len(B)
        if r.method == 'CPS-2F':
            XA, _ = extended_curve(A, r.delta1)
            XB, _ = extended_curve(B, r.delta2)
            ma, nb = m + 2*m*(m-1), n + 2*n*(n-1)
        else:
            XA, XB, ma, nb = A, B, m, n
        actual = valid_vertex_count(A, B, XA, XB, r.delta1, r.delta2, r.delta3)
        candidate = len(XA)*m*len(XB)*n
        upper = ma*m*nb*n
        assert 0 <= actual <= candidate <= upper
        if pd.notna(r.states):
            assert r.states <= actual
        rows.append(dict(pair=r.pair, w=r.w, alpha=r.alpha, method=r.method,
            nA=m, nB=n, delta1=r.delta1, delta2=r.delta2, delta3=r.delta3,
            auxiliary_A=len(XA)-m, auxiliary_B=len(XB)-n,
            auxiliary_upper_A=ma-m, auxiliary_upper_B=nb-n,
            extended_A=len(XA), extended_B=len(XB),
            valid_vertices=actual, cartesian_bound=candidate, theoretical_vertex_bound=upper,
            fraction_of_bound=actual/upper, discovered_states=r.states,
            examined_transitions=r.transitions, status=r.status))
    table = pd.DataFrame(rows)
    table.to_csv(out/'graph_sizes.csv', index=False)
    document = '<!doctype html><meta charset="utf-8"><h1>Configuration graph sizes</h1>'
    document += ('<p>Exact valid vertices before reachability pruning; theory is an upper bound. '
                 'Discovered states and examined transitions are solver work counters, not full graph sizes. '
                 'Certificate runs have zero work counters; missing counters were not recorded. '
                 'Full graph edges are not counted. All distances are in angstroms.</p>')
    document += table.to_html(index=False, float_format=lambda x: f'{x:.5g}')
    (out/'graph_sizes.html').write_text(document, encoding='utf-8')
    return table


def graph_size_plots(table, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({'font.family':'serif', 'mathtext.fontset':'stix',
                         'pdf.fonttype':42, 'svg.fonttype':'none', 'font.size':10}):
        for (w, alpha), group in table.groupby(['w', 'alpha']):
            fig, axes = plt.subplots(1, 2, figsize=(16, 9), sharey=True)
            for ax, method in zip(axes, ['CPS-3F', 'CPS-2F']):
                g = group[group.method.eq(method)]
                x = np.arange(len(g))
                exact_label = (r'Valid vertices: $\sum_{(p,q)\in C_3}c_A(p)c_B(q)$'
                               '\nAll three distance constraints; before reachability pruning')
                if method == 'CPS-3F':
                    width = .34
                    bars = [(-.5, 'valid_vertices', exact_label, '#4477AA', ''),
                            (.5, 'theoretical_vertex_bound',
                             r'Theoretical bound: $m^2n^2$'+'\nOriginal owner positions; no auxiliary points', '#EE7733', '..')]
                else:
                    width = .25
                    bars = [(-1, 'valid_vertices', exact_label, '#4477AA', ''),
                            (0, 'cartesian_bound', r'Cartesian bound: $m^*m\,n^*n$'
                             +'\nMeasured extended sizes, including auxiliary points', '#228833', '//'),
                            (1, 'theoretical_vertex_bound',
                             r'Worst-case bound: $[m+2m(m-1)]m[n+2n(n-1)]n$'
                             +'\n'+r'At most $2r(r-1)$ auxiliary points per $r$-vertex curve; $O(m^3n^3)$', '#EE7733', '..')]
                for shift, field, label, color, hatch in bars:
                    ax.bar(x+shift*width, g[field], width, label=label, color=color, hatch=hatch)
                ax.set_yscale('log'); ax.set_title(method)
                ax.set_xticks(x, g.pair.str.replace('__',' / ',regex=False), rotation=90)
                ax.grid(axis='y', which='major', alpha=.2); ax.set_axisbelow(True)
                ax.spines[['top','right']].set_visible(False)
                ax.legend(loc='lower left', bbox_to_anchor=(0,1.10), frameon=False,
                          fontsize=9, borderaxespad=0, labelspacing=.8)
            axes[0].set_ylabel(r'Configuration graph vertices $|V(G)|$ (log scale)')
            fig.suptitle(f'Graph size versus theoretical upper bounds — w={w}, α={alpha:g}')
            fig.text(.5,.02,
                     r'$m=|A|,\ n=|B|,\ m^*=|A^*|,\ n^*=|B^*|;\quad C_3=\{(p,q):\|a_p-b_q\|\leq\delta_3\}$.'
                     +'\n'+r'$c_A(p)$ counts owner positions within $\delta_1$ of $a_p$; $c_B(q)$ uses $\delta_2$.'
                     +'\nOwners are original vertices for CPS-3F and extended-curve vertices for CPS-2F.'
                     +'\nExact pair-specific δ₁, δ₂, δ₃ and numerical counts: graph_sizes.csv. Bounds are upper limits, not predictions.',
                     ha='center', fontsize=9)
            fig.subplots_adjust(left=.07,right=.99,top=.68,bottom=.32,wspace=.12)
            stem = out/f'graph_sizes_w{w}_alpha{alpha:g}'
            for ext in ['pdf', 'svg', 'png']:
                fig.savefig(stem.with_suffix('.'+ext), dpi=300, bbox_inches='tight')
            plt.show()
            plt.close(fig)


def check_graph_counts():
    """Independent brute-force quadruple count on small 3D curves."""
    from itertools import product
    A = np.array([[0.,0,0],[1,1,0],[2,0,0]])
    B = np.array([[0.,.2,0],[1,.8,0],[2,.2,0]])
    checks = 0
    for d in [.3, 1.1]:
        for continuous in [False, True]:
            XA = extended_curve(A,d)[0] if continuous else A
            XB = extended_curve(B,d)[0] if continuous else B
            brute = sum(np.linalg.norm(XA[i]-A[p]) <= d+EPS and
                        np.linalg.norm(XB[j]-B[q]) <= d+EPS and
                        np.linalg.norm(A[p]-B[q]) <= .9+EPS
                        for i,p,j,q in product(range(len(XA)),range(len(A)),range(len(XB)),range(len(B))))
            assert brute == valid_vertex_count(A,B,XA,XB,d,d,.9)
            checks += 1
    return checks


def auxiliary_utilization_plots(table, out):
    """Compare actual added owner positions with the per-curve auxiliary bound."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    columns = ['pair','w','alpha','nA','nB','delta1','delta2','delta3',
               'auxiliary_A','auxiliary_upper_A','auxiliary_B','auxiliary_upper_B']
    utilization = table.loc[table.method.eq('CPS-2F'), columns].copy()
    for curve in ['A', 'B']:
        bound = utilization[f'auxiliary_upper_{curve}']
        utilization[f'utilization_pct_{curve}'] = 100 * utilization[f'auxiliary_{curve}'] / bound.where(bound > 0)
        values = utilization[f'utilization_pct_{curve}'].dropna()
        assert values.between(0,100).all()
    utilization.to_csv(out.parent/'auxiliary_utilization.csv', index=False)
    with plt.rc_context({'font.family':'serif','mathtext.fontset':'stix',
                         'pdf.fonttype':42,'svg.fonttype':'none','font.size':10}):
        for (w, alpha), g in utilization.groupby(['w','alpha']):
            fig, ax = plt.subplots(figsize=(12,7))
            x = np.arange(len(g)); width = .36
            for offset, curve, color, hatch in [(-.5,'A','#4477AA',''),(.5,'B','#228833','//')]:
                values = g[f'utilization_pct_{curve}']
                bars = ax.bar(x+offset*width, values, width, label=f'Curve {curve}', color=color, hatch=hatch)
                ax.bar_label(bars, labels=[f'{v:.1f}%' if pd.notna(v) else '' for v in values],
                             fontsize=8, padding=3, rotation=90)
            ax.set_xticks(x,g.pair.str.replace('__',' / ',regex=False),rotation=90)
            ax.set_ylabel(r'Utilization percentage: $100\times\mathrm{auxiliary\ points}/\mathrm{theoretical\ bound}$')
            top = g[['utilization_pct_A','utilization_pct_B']].max().max()
            ax.set_ylim(0, max(5, float(top)*1.3) if pd.notna(top) else 5)
            ax.set_title(f'CPS-2F auxiliary-point utilization — w={w}, α={alpha:g}')
            ax.legend(frameon=False); ax.grid(axis='y',alpha=.2); ax.set_axisbelow(True)
            ax.spines[['top','right']].set_visible(False)
            fig.text(.5,.02,
                r'Curve A: $100a_A/[2m(m-1)]$; curve B: $100a_B/[2n(n-1)]$.'
                +'\n'+r'$a_A,a_B$: added auxiliary points only; $m=|A|$, $n=|B|$. A zero bound gives N/A.'
                +'\nCPS-3F adds no auxiliary points and is omitted. Exact δ₁, δ₂, δ₃: auxiliary_utilization.csv.',
                ha='center',fontsize=9)
            fig.subplots_adjust(left=.09,right=.99,top=.91,bottom=.36)
            stem = out/f'auxiliary_utilization_w{w}_alpha{alpha:g}'
            for ext in ['pdf','svg','png']:
                fig.savefig(stem.with_suffix('.'+ext),dpi=300,bbox_inches='tight')
            plt.show()
            plt.close(fig)
    return utilization
