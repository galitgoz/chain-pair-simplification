from pathlib import Path
import html
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from decimation_algorithms import edge_stats,shortcut_dag,shortest_path,joint_local,discrete_frechet

STRIDES=[1,2,4,8,16]
ROLE={1:'w=1 control',2:'w=2 headline',4:'w=4 headline',8:'w=8 RUNTIME ONLY',16:'w=16 RUNTIME ONLY'}
STYLES={1:'-',2:'--',4:':',8:'-.',16:(0,(5,2,1,2,1,2))}
COLORS={'separate':'#17689b','CPS':'#cd6229'}

def full_table(frame):
    groups=[('Instance',[('B_printed','B printed'),('w','w'),('stride_role','Role'),('alpha','α')]),
        ('Before stride decimation (qualified input)',[(c,c) for c in ['nA_before','nB_before','edge_mean_before','dF_before']]),
        ('After decimation',[(c,c) for c in ['nA_after','nB_after','edge_mean_after','edge_min_after','edge_max_after','keep_frac','dF_after','dF_change','edge_growth_factor','appended_stub','last_stride','coarse_trace']]),
        ('Thresholds / prechecks',[(c,c) for c in ['delta1','delta2','delta3','floor','feasible','local_fidelity_floor','solver_feasible']]),
        ('Separate',[(c,c) for c in ['sep_kA','sep_kB','sep_equal','sep_lenratio_A','sep_dF','sep_max_link','sep_ok3','sep_coupling']]),
        ('Joint CPS',[(c,c) for c in ['cps_k','cps_lenratio_A','cps_dF','cps_max_link','cost_of_joining','cps_status','cps_verified']])]
    parts=['<div style="overflow:auto;max-height:760px;border:1px solid #b9c6d1"><table style="border-collapse:collapse;font:12px system-ui;white-space:nowrap">',
        '<thead style="position:sticky;top:0;background:#e6eef6;z-index:2"><tr>']
    for name,cols in groups:parts.append(f'<th colspan="{len(cols)}" style="padding:9px;border:1px solid #bdc9d4">{html.escape(name)}</th>')
    cols=[v for _,cc in groups for v in cc];parts.append('</tr><tr>')
    for c,label in cols:parts.append(f'<th style="padding:7px;border:1px solid #bdc9d4">{html.escape(label)}</th>')
    parts.append('</tr></thead><tbody>')
    for row in frame.itertuples():
        bg='#f0f0f0' if row.w>=8 else ('#fff' if row.solver_feasible else '#fce5df')
        parts.append(f'<tr style="background:{bg}">')
        for col,label in cols:
            v=getattr(row,col)
            if pd.isna(v):s=''
            elif isinstance(v,(bool,np.bool_)):s='yes' if v else 'no'
            elif isinstance(v,(float,np.floating)):s=f'{v:.3f}' if np.isfinite(v) else '∞'
            else:s=str(v)
            style=''
            if col=='sep_max_link' and row.sep_equal and not row.sep_ok3:style=';background:#ffdcb8;font-weight:700'
            if col=='w' and row.w==1:style=';font-weight:700;background:#daecfb'
            parts.append(f'<td style="padding:7px;border:1px solid #d7dfe5{style}">{html.escape(s)}</td>')
        parts.append('</tr>')
    parts.append('</tbody></table></div><p>Å throughout; length ratios are arc-length ratios to the decimated input. Blank separate link/check fields mean unequal output sizes, not a replacement by dF. Gray w=8/16 rows are RUNTIME ONLY; w=1 is the undecimated control.</p>')
    return ''.join(parts)

def F1(frame,directory):
    directory=Path(directory);keys=list(frame.pair.drop_duplicates())
    for alpha in [1,2,4]:
        fig,axes=plt.subplots(4,4,figsize=(16,13),sharey=True)
        fig.subplots_adjust(left=.06,right=.99,bottom=.08,top=.85,wspace=.2,hspace=.5)
        for ax,key in zip(axes.flat,keys):
            data=frame[(frame.pair==key)&(frame.alpha==alpha)]
            for row in data.itertuples():
                x=np.array([0,1,2],float)
                ax.plot(x,[row.nA_before,row.nA_after,row.sep_kA],color=COLORS['separate'],marker='o',ms=3.5,lw=1.3,ls=STYLES[row.w],alpha=.85)
                if row.solver_feasible:
                    ax.plot(x,[row.nA_before,row.nA_after,row.cps_k],color=COLORS['CPS'],marker='s',ms=3.2,lw=1.1,ls=STYLES[row.w],alpha=.8)
            ax.set_yscale('log');ax.set_xticks([0,1,2],['raw*','decimated','simplified'],rotation=18)
            ax.set_title(data.B_verified.iloc[0],fontsize=10);ax.grid(alpha=.2,which='both')
        for ax in list(axes.flat)[len(keys):]:ax.set_visible(False)
        for ax in axes[:,0]:ax.set_ylabel('Vertices of A (log scale)')
        handles=[Line2D([0],[0],color=COLORS['separate'],marker='o',label='Separate'),Line2D([0],[0],color=COLORS['CPS'],marker='s',label='Joint CPS')]
        handles += [Line2D([0],[0],color='#555',ls=STYLES[w],label=ROLE[w]) for w in STRIDES]
        fig.legend(handles=handles,ncol=4,loc='upper center',bbox_to_anchor=(.5,.945),fontsize=9)
        fig.suptitle(f'F1 · vertices through the pipeline · α={alpha}\nOne style per stride; colour AND marker distinguish methods',fontsize=14,y=.995)
        fig.text(.06,.018,'* raw = qualified, aligned paired input; qualification exclusions are reported separately. Headline: w≤4. w=8/16: RUNTIME ONLY.',fontsize=10)
        fig.savefig(directory/f'F1_vertices_alpha{alpha}.png',dpi=150);plt.show()

def F2(frame,directory):
    data=frame.drop_duplicates(['pair','w']);palette={'similar':'#247b77','intermediate':'#d18026','dissimilar':'#8d5189'}
    markers=['o','s','^','D','v','P','X','<','>','h','*','p','8','d']
    fig,ax=plt.subplots(figsize=(13,7));fig.subplots_adjust(left=.09,right=.69,bottom=.17,top=.89)
    ax.axvspan(5.65,18,color='#ececec',zorder=0)
    for marker,(key,g) in zip(markers,data.groupby('pair',sort=False)):
        g=g.sort_values('w')
        ax.plot(g.w,g.dF_after,color=palette[g.regime.iloc[0]],marker=marker,ms=6,lw=1.4,
            label=f'{g.B_verified.iloc[0]} ({g.regime.iloc[0]})')
    ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xlim(.9,18)
    ax.set_xticks(STRIDES,['1\ncontrol','2','4','8\nruntime only','16\nruntime only'])
    ax.set(xlabel='Uniform stride w (log scale)',ylabel='dDF after decimation (Å, log scale)',
        title='F2 · how stride changes the measured separation')
    ax.grid(alpha=.2,which='both');ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),fontsize=9)
    fig.text(.09,.025,'Regimes use fixed pre-decimation US-align TM score: similar ≥0.5; intermediate 0.3–0.5; dissimilar <0.3.\nA distance decrease is not automatically better fidelity: decimation may remove the residues exposing disagreement.',fontsize=10)
    fig.savefig(Path(directory)/'F2_distance_vs_stride.png',dpi=150);plt.show()

def F3(frame,directory):
    keys=list(frame.pair.drop_duplicates());columns=[(w,a) for w in STRIDES for a in [1,2,4]]
    values=np.zeros((len(keys),len(columns)),int)
    for i,key in enumerate(keys):
        for j,(w,a) in enumerate(columns):
            r=frame[(frame.pair==key)&(frame.w==w)&(frame.alpha==a)].iloc[0]
            values[i,j]=0 if not r.sep_equal else (2 if r.sep_ok3 else 1)
    colors=['#dce1e5','#f3c087','#9bd3c6'];symbols=['×','!','✓']
    fig,ax=plt.subplots(figsize=(16,8));fig.subplots_adjust(left=.10,right=.99,top=.87,bottom=.20)
    ax.imshow(values,cmap=ListedColormap(colors),vmin=0,vmax=2,aspect='auto')
    for i in range(len(keys)):
        for j in range(len(columns)):ax.text(j,i,symbols[values[i,j]],ha='center',va='center',fontsize=13)
    ax.set_yticks(range(len(keys)),[frame[frame.pair==key].B_verified.iloc[0] for key in keys])
    ax.set_xticks(range(len(columns)),[f'{w}:{a}'+('\nRUNTIME' if w>=8 else ('\ncontrol' if w==1 else '')) for w,a in columns],fontsize=9)
    for x in [2.5,5.5,8.5,11.5]:ax.axvline(x,color='white',lw=3)
    ax.set(title='F3 · did separate simplification produce a usable index-wise coupling?',xlabel='Column = w : α; w=1 is control; w=8/16 are RUNTIME ONLY')
    handles=[Line2D([0],[0],color=c,marker='s',ms=14,lw=0,label=f'{s} {label}') for c,s,label in zip(colors,symbols,['no coupling (unequal k)','coupled, violates δ₃','coupled, meets δ₃'])]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.025),ncol=3,fontsize=10)
    fig.savefig(Path(directory)/'F3_separate_coupling.png',dpi=150);plt.show()

def F4(pairs,directory):
    """Display-only 20-residue computation; never append values to result tables."""
    pair=next((p for p in pairs.values() if p['regime']=='dissimilar'),next(iter(pairs.values())))
    n=min(20,len(pair['A']));start=(len(pair['A'])-n)//2
    A=pair['A'][start:start+n];B=pair['B'][start:start+n]
    e=edge_stats(A,B)['mean'];d3=float(np.ceil(discrete_frechet(A,B)))
    alpha=2;EA=shortcut_dag(A,alpha*e);EB=shortcut_dag(B,alpha*e)
    ia=shortest_path(EA);ib=shortest_path(EB);joint=joint_local(A,B,EA,EB,d3)
    if joint['ia'] is None:raise ValueError('Predeclared F4 display window is infeasible; do not insert a fabricated joint drawing')
    combined=np.vstack([A,B]);center=combined.mean(0);_,s,Vt=np.linalg.svd(combined-center,full_matrices=False)
    variance=float(np.sum(s[:2]**2)/np.sum(s**2));pa=(A-center)@Vt[:2].T;pb=(B-center)@Vt[:2].T
    lower=np.vstack([pa,pb]).min(0);upper=np.vstack([pa,pb]).max(0);pad=(upper-lower).max()*.08
    fig,axes=plt.subplots(1,3,figsize=(15,5.5));fig.subplots_adjust(left=.06,right=.98,bottom=.23,top=.80,wspace=.25)
    for ax,P,ids,color,name in [(axes[0],pa,ia,COLORS['separate'],'A alone'),(axes[1],pb,ib,COLORS['CPS'],'B alone')]:
        ax.plot(*P.T,color=color,lw=.9,alpha=.45,label='original 20-residue window')
        ax.plot(*P[ids].T,color=color,lw=2,marker='o',ms=4,label='separate simplification')
        ax.set_title(name);ax.legend(fontsize=8)
    ja=joint['ia'];jb=joint['ib'];ax=axes[2]
    for P,ids,color,label in [(pa,ja,COLORS['separate'],'A joint'),(pb,jb,COLORS['CPS'],'B joint')]:
        ax.plot(*P.T,color=color,lw=.7,alpha=.2);ax.plot(*P[ids].T,color=color,lw=1.8,marker='s',ms=4,label=label)
    for a,b in zip(pa[ja],pb[jb]):ax.plot([a[0],b[0]],[a[1],b[1]],'--',color='#555',lw=.85)
    ax.set_title('Joint result with index-wise links');ax.legend(fontsize=8)
    for ax in axes:
        ax.set_aspect('equal',adjustable='box');ax.set(xlim=(lower[0]-pad,upper[0]+pad),ylim=(lower[1]-pad,upper[1]+pad),xlabel='PC1 (Å)',ylabel='PC2 (Å)');ax.grid(alpha=.15)
    fig.suptitle(f'F4 · display-only, w=1 · {pair["A_verified"]} / {pair["B_verified"]}\n{n} consecutive qualified residues, indices {start}–{start+n-1}',fontsize=13)
    fig.text(.06,.055,f'Common PCA of both 3D windows retains {variance:.2%} of coordinate variance. α=2, δ₁=δ₂={alpha*e:.2f} Å, δ₃={d3:g} Å.\nAll solving is in 3D. This projection can shorten apparent links. No F4 measurement enters a results table.',fontsize=10)
    fig.savefig(Path(directory)/'F4_display_window.png',dpi=160);plt.show()

def headline_links(frame,directory):
    data=frame[(frame.w<=4)&frame.sep_equal&frame.solver_feasible]
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for ax,w in zip(axes,[1,2,4]):
        g=data[data.w==w]
        for alpha,color in zip([1,2,4],['#17689b','#d26526','#44855c']):
            h=g[g.alpha==alpha];ax.scatter(h.sep_max_link,h.cps_max_link,label=f'α={alpha}',color=color,s=30)
        lim=max([1.]+g.sep_max_link.tolist()+g.cps_max_link.tolist())*1.05
        ax.plot([0,lim],[0,lim],'--',color='gray',lw=1)
        ax.set(xlim=(0,lim),ylim=(0,lim),xlabel='Separate max link (Å)',ylabel='CPS max link (Å)',title=ROLE[w]);ax.legend(fontsize=8)
    fig.suptitle('Headline: the index-wise link constraint, compared at the same δ₁\nOnly equal-size separate outputs and feasible joint rows are plottable; all other statuses remain in F3',fontsize=11)
    fig.savefig(Path(directory)/'headline_max_links.png',dpi=150);plt.show()
