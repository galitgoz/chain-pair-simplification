"""Aligned 3D windows, complete pipeline tables, and outcome diagnostics."""
from pathlib import Path
import html
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

VIEWS=[(25,35),(25,125)]
COLORS=['#17689b','#d26526']

def axes3d(ax,half,view):
    ax.set(xlim=(-half,half),ylim=(-half,half),zlim=(-half,half),xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)')
    ax.set_box_aspect((1,1,1));ax.view_init(*view);ax.tick_params(labelsize=7,pad=0)

def plot_aligned(pairs,directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    items=list(pairs.items())
    centers={key:np.vstack([p['A'],p['B']]).mean(0) for key,p in items}
    half=max(np.abs(np.vstack([p['A'],p['B']])-centers[key]).max() for key,p in items)*1.08
    for page in range(0,len(items),4):
        subset=items[page:page+4];fig=plt.figure(figsize=(12,4.3*len(subset)))
        fig.subplots_adjust(left=.04,right=.96,bottom=.05,top=.92,wspace=.15,hspace=.40)
        for r,(key,p) in enumerate(subset):
            for c,view in enumerate(VIEWS):
                ax=fig.add_subplot(len(subset),2,r*2+c+1,projection='3d')
                for side,color in zip(['A','B'],COLORS):
                    P=p[side]-centers[key]
                    ax.plot(*P.T,color=color,lw=1.5,marker='.',markersize=3,label=p[f'{side}_key'])
                    ax.scatter(*P[0],color=color,s=32,marker='o',edgecolor='black',linewidth=.5)
                    ax.scatter(*P[-1],color=color,s=34,marker='s',edgecolor='black',linewidth=.5)
                axes3d(ax,half,view)
                ax.set_title(f'{p["A_key"]} / {p["B_key"]} · n=30/30 · az={view[1]}°',fontsize=10)
                ax.legend(loc='upper left',fontsize=8)
        fig.suptitle('Structurally superimposed original windows\nShared center per pair; common viewpoints and one common Å scale',fontsize=12)
        fig.savefig(directory/f'aligned_windows_{page//4+1}.png',dpi=130)
        plt.show()

def pipeline_html(frame):
    groups=[('Pair / threshold',[('B_printed','B printed'),('pair','Verified pair'),('delta1','δ₁=δ₂'),('status','Status'),('feasibility_floor','Floor')]),
        ('Originals',[('original_n_A','n A'),('original_n_B','n B'),('original_arc_A','Arc A'),('original_arc_B','Arc B'),('original_dDF','dDF(A,B)'),('delta3','δ₃ used'),('paper_ceil_delta3','ceil dDF')]),
        ('Separate',[('separate_k_A','k A'),('separate_k_B','k B'),('separate_arc_A','Arc A′'),('separate_arc_B','Arc B′'),('separate_arc_ratio_A','A′/A'),('separate_arc_ratio_B','B′/B'),('separate_dDF','dDF(A′,B′)')]),
        ('Separate coupling',[('separate_coupled','Equal k?'),('separate_coupling_status','Coupling'),('separate_index_max','max index d'),('separate_index_le_delta3','≤ δ₃?')]),
        ('Synchronous discrete CPS',[('CPS_k','k'),('CPS_arc_A','Arc A′'),('CPS_arc_B','Arc B′'),('CPS_arc_ratio_A','A′/A'),('CPS_arc_ratio_B','B′/B'),('CPS_index_max','max index d'),('CPS_dDF','dDF(A′,B′)'),('CPS_solver_verified','Solver verified')]),
        ('Cost of joining',[('cost_of_joining','k CPS − max k separate')])]
    out=['<div style="overflow:auto;max-height:760px;border:1px solid #ccd5dd"><table style="border-collapse:collapse;font:12px system-ui;white-space:nowrap">',
         '<thead style="position:sticky;top:0;background:#e7edf4;z-index:2"><tr>']
    for group,columns in groups: out.append(f'<th colspan="{len(columns)}" style="padding:9px;border:1px solid #bcc9d5">{html.escape(group)}</th>')
    out.append('</tr><tr>')
    cols=[item for _,columns in groups for item in columns]
    for col,label in cols: out.append(f'<th style="padding:8px;border:1px solid #bcc9d5">{html.escape(label)}</th>')
    out.append('</tr></thead><tbody>')
    for _,row in frame.iterrows():
        background='#fff' if row.included_in_aggregates else '#eeeeee'
        out.append(f'<tr style="background:{background}">')
        for col,label in cols:
            value=row[col]
            if pd.isna(value): rendered=''
            elif isinstance(value,(bool,np.bool_)): rendered='yes' if value else 'no'
            elif isinstance(value,(float,np.floating)): rendered=f'{value:.3f}' if np.isfinite(value) else '∞'
            else: rendered=str(value)
            emphasis=';background:#ffe0c8;font-weight:600' if col=='separate_dDF' and row.separate_further_apart else ''
            out.append(f'<td style="padding:7px;border:1px solid #d8dfe5{emphasis}">{html.escape(rendered)}</td>')
        out.append('</tr>')
    out.append('</tbody></table></div><p>Distances and arc lengths: Å. Ratios: dimensionless. Blank index-wise fields mean no equal-size coupling. Gray rows are infeasible and excluded from aggregates. Orange cells exceed the originals’ dDF.</p>')
    return ''.join(out)

def plot_outcomes(frame,directory):
    directory=Path(directory);feasible=frame[frame.included_in_aggregates]
    fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
    for delta,color in zip([2,4,8],['#17689b','#d26526','#3e8d57']):
        g=feasible[feasible.delta1==delta]
        axes[0].scatter(g.original_dDF,g.separate_dDF,label=f'δ={delta} Å',color=color,s=36,alpha=.8)
        axes[1].scatter(g.original_dDF,g.CPS_dDF,color=color,s=36,alpha=.8)
        axes[2].scatter(g.delta1+np.linspace(-.13,.13,len(g)),g.cost_of_joining,color=color,s=34,alpha=.8)
    limit=max(frame.original_dDF.max(),frame.separate_dDF.max())*1.06
    for ax in axes[:2]:
        ax.plot([0,limit],[0,limit],'--',color='#6b7279',lw=1)
        ax.set(xlim=(0,limit),ylim=(0,limit),xlabel='Original dDF(A,B) (Å)',ylabel='Simplified dDF(A′,B′) (Å)')
    axes[0].set_title('Independent can increase separation');axes[0].legend(fontsize=8)
    axes[1].set_title('Joint bound: solver verification')
    axes[2].set(xlabel='δ₁ = δ₂ (Å)',ylabel='k CPS − max(k A, k B)',title='Requested cost of joining',xticks=[2,4,8])
    axes[2].axhline(0,color='gray',lw=1)
    fig.suptitle(f'Feasible rows only ({len(feasible)}/{len(frame)}); infeasible rows excluded',fontsize=12)
    fig.savefig(directory/'comparison_outcomes.png',dpi=150);plt.show()
    # Display every pair/threshold, including infeasible cells, to prevent survivorship concealment.
    pivot=frame.pivot(index='B_printed',columns='delta1',values='included_in_aggregates')
    fig,ax=plt.subplots(figsize=(7,7),layout='constrained')
    ax.imshow(pivot.to_numpy(dtype=int),cmap=plt.matplotlib.colors.ListedColormap(['#dce1e5','#d2eadc']),vmin=0,vmax=1,aspect='auto')
    for r,key in enumerate(pivot.index):
        for c,delta in enumerate(pivot.columns):
            row=frame[(frame.B_printed==key)&(frame.delta1==delta)].iloc[0]
            label=f'k={int(row.CPS_k)}' if row.included_in_aggregates else f'infeasible\nfloor={row.feasibility_floor:.2f}'
            ax.text(c,r,label,ha='center',va='center',fontsize=8)
    ax.set(xticks=range(3),xticklabels=['2 Å','4 Å','8 Å'],yticks=range(len(pivot)),yticklabels=pivot.index,
        xlabel='δ₁ = δ₂',title='All pair/threshold statuses (printed B identifiers)')
    fig.savefig(directory/'feasibility_status.png',dpi=150);plt.show()

def plot_case(pairs,frame,paths,directory):
    eligible=frame[frame.included_in_aggregates]
    # Selection AFTER analysis for illustration only; the complete unfiltered table remains primary.
    row=eligible.sort_values('separate_dDF_minus_original',ascending=False).iloc[0]
    p=pairs[row.pair];path=paths[f'{row.pair}@{row.delta1:g}'];center=np.vstack([p['A'],p['B']]).mean(0)
    half=np.abs(np.vstack([p['A'],p['B']])-center).max()*1.1
    fig=plt.figure(figsize=(14,9));fig.subplots_adjust(left=.03,right=.97,bottom=.07,top=.88,wspace=.15,hspace=.35)
    for c,method in enumerate(['original','separate','CPS']):
        for r,view in enumerate(VIEWS):
            ax=fig.add_subplot(2,3,r*3+c+1,projection='3d')
            for side,color in zip(['A','B'],COLORS):
                P=p[side]-center
                ids=np.arange(len(P)) if method=='original' else path[f'{method}_{side}']
                if method!='original':ax.plot(*P.T,color=color,lw=.7,alpha=.2)
                ax.plot(*P[ids].T,color=color,lw=1.8,marker='o',markersize=3,label=f'{side}: {len(ids)} vertices')
            axes3d(ax,half,view);ax.legend(fontsize=7,loc='upper left')
            distance=row.original_dDF if method=='original' else row[f'{method}_dDF']
            ax.set_title(f'{method} · dDF={distance:.3f} Å',fontsize=10)
    fig.suptitle(f'{p["A_key"]} / {p["B_key"]}, δ={row.delta1:g} Å\nLargest observed independent distance increase among feasible rows (illustration, not sample selection)',fontsize=12)
    fig.savefig(Path(directory)/'illustrative_pipeline.png',dpi=140);plt.show()
    return row
