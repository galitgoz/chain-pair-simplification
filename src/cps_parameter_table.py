"""A single, prominent table of input sizes, budgets and CPS-2F auxiliary counts."""
from pathlib import Path
import pandas as pd
import numpy as np


def configured_parameter_table(data,alphas,out):
    """Input-only diagnostics: available even if the joint solvers hit a limit."""
    from cps_paper_algorithms import extended_curve
    from curve_algorithms import discrete_frechet
    rows=[]
    for item in data['instances'].values():
        A=item['A'];B=item['B_xyz']
        mean_edge=.5*(np.linalg.norm(np.diff(A,axis=0),axis=1).mean()+np.linalg.norm(np.diff(B,axis=0),axis=1).mean())
        d3=float(np.ceil(discrete_frechet(A,B)))
        for alpha in alphas:
            delta=alpha*mean_edge
            na=len(A);nb=len(B)
            a=len(extended_curve(A,delta)[0])-na;b=len(extended_curve(B,delta)[0])-nb
            rows.append(dict(pair=item['pair'],w=item['w'],alpha=alpha,nA=na,nB=nb,
                delta1=delta,delta2=delta,delta3=d3,auxiliary_A=a,auxiliary_B=b,
                auxiliary_total=a+b,extended_A=na+a,extended_B=nb+b,method='CPS-2F'))
    return parameter_table(pd.DataFrame(rows),out)


def parameter_table(results,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    columns=['pair','w','alpha','nA','nB','delta1','delta2','delta3',
             'auxiliary_A','auxiliary_B','auxiliary_total','extended_A','extended_B']
    table=results.loc[results.method.eq('CPS-2F'),columns].copy()
    assert not table.duplicated(['pair','w','alpha']).any()
    assert table[['auxiliary_A','auxiliary_B','auxiliary_total']].notna().all().all()
    assert (table.auxiliary_total==table.auxiliary_A+table.auxiliary_B).all()
    # Thesis Observation 4: at most two sphere/edge intersections per pair.
    table['auxiliary_upper_A']=2*table.nA*(table.nA-1)
    table['auxiliary_upper_B']=2*table.nB*(table.nB-1)
    table['auxiliary_upper_total']=table.auxiliary_upper_A+table.auxiliary_upper_B
    assert (table.auxiliary_A<=table.auxiliary_upper_A).all()
    assert (table.auxiliary_B<=table.auxiliary_upper_B).all()
    table=table[['pair','w','alpha','nA','nB','delta1','delta2','delta3',
                 'auxiliary_A','auxiliary_upper_A','auxiliary_B','auxiliary_upper_B',
                 'auxiliary_total','auxiliary_upper_total','extended_A','extended_B']]
    names={'pair':'Protein pair (A / B)','w':'Stride w','alpha':'α',
           'nA':'|A|','nB':'|B|','delta1':'δ1 (Å)','delta2':'δ2 (Å)','delta3':'δ3 (Å)',
           'auxiliary_A':'Auxiliary points A','auxiliary_B':'Auxiliary points B',
           'auxiliary_total':'Auxiliary points total','extended_A':'|A*|','extended_B':'|B*|',
           'auxiliary_upper_A':'Theoretical aux upper bound A',
           'auxiliary_upper_B':'Theoretical aux upper bound B',
           'auxiliary_upper_total':'Theoretical aux upper bound total'}
    integer_columns=['w','nA','nB','auxiliary_A','auxiliary_B','auxiliary_total','extended_A','extended_B',
                     'auxiliary_upper_A','auxiliary_upper_B','auxiliary_upper_total']
    table[integer_columns]=table[integer_columns].astype(int)
    table['pair']=table['pair'].str.replace('__',' / ',regex=False)
    table.rename(columns=names).to_csv(out/'pair_parameters.csv',index=False)
    pieces=['''<!doctype html><meta charset="utf-8"><title>Protein pair parameters and auxiliary counts</title>
<style>body{font:15px Georgia,serif;color:#192a35;margin:24px}h1{font-size:24px}h2{font-size:18px;margin-top:30px}
table{border-collapse:collapse;font-size:13px}th,td{padding:7px 9px;border-bottom:1px solid #d6dde2;text-align:right;white-space:nowrap}
th{background:#edf2f5}td:first-child,th:first-child{text-align:left}.scroll{overflow-x:auto}tr:nth-child(even){background:#f8fafb}</style>
<h1>Per-pair input sizes, distance budgets and CPS-2F auxiliary points</h1>
<p>All distances are in Å. Auxiliary points A and B are additional owner positions in the CPS-2F extended curves;
they exclude the original vertices. |A*| = |A| + auxiliary A; |B*| = |B| + auxiliary B.
Counts are evaluated for every case, including cases solved by an independent-bound certificate.
They are not selected simplification vertices. CPS-3F adds no auxiliary points.</p>
<p><b>Theoretical counts are worst-case upper bounds, not predicted or exact counts.</b>
For a curve with n input vertices and n−1 edges, at most 2n(n−1) auxiliary points arise:
each vertex-centered sphere intersects each nondegenerate edge at most twice.
Thus U<sub>A</sub>=2|A|(|A|−1) and U<sub>B</sub>=2|B|(|B|−1).
These bounds exclude original vertices; extended sizes satisfy |A*|≤|A|+U<sub>A</sub>
and |B*|≤|B|+U<sub>B</sub>. The bounds depend on input size; measured counts also depend
on δ1/δ2 and curve geometry. Assumption: finite sphere–edge intersections (no edge contained
in a sphere boundary). Source: Galit Gozoltzani thesis, Chapter 3, Observation 4.</p>''']
    for (w,alpha),g in table.groupby(['w','alpha'],sort=True):
        pieces.append(f'<h2>Stride w={w}; α={alpha:g} — {len(g)} pairs</h2><div class="scroll">')
        pieces.append(g.drop(columns=['w','alpha']).rename(columns=names).to_html(index=False,
            float_format=lambda x:f'{x:.3f}',border=0))
        pieces.append('</div>')
    document='\n'.join(pieces)
    (out/'pair_parameters.html').write_text(document,encoding='utf-8')
    return table,document
