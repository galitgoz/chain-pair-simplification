"""Shared plotting helpers: axes are in angstroms and never bridge fragments."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

VIEWS=[(25,45),(25,135)]

def equal_axes(ax,half,view):
    ax.set(xlim=(-half,half),ylim=(-half,half),zlim=(-half,half),xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)')
    ax.set_box_aspect((1,1,1)); ax.view_init(*view)
    ax.tick_params(labelsize=7,pad=0)

def plot_fragments(chains,directory):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    keys=list(chains)
    centers={key:np.mean(c['xyz'],axis=0) for key,c in chains.items()}
    half=max(np.max(np.abs(c['xyz']-centers[key])) for key,c in chains.items())*1.07
    colors=plt.get_cmap('tab10')
    for page in range(0,len(keys),5):
        fig=plt.figure(figsize=(12,22))
        fig.subplots_adjust(left=.03,right=.97,bottom=.03,top=.94,wspace=.12,hspace=.38)
        for row,key in enumerate(keys[page:page+5]):
            c=chains[key]
            for col,view in enumerate(VIEWS):
                ax=fig.add_subplot(5,2,row*2+col+1,projection='3d')
                for f,fragment in enumerate(c['fragments']):
                    P=fragment['xyz']-centers[key]
                    ax.plot(*P.T,color=colors(f),lw=1.15,label=f'f{f}: {len(P)} Cα')
                    ax.scatter(*P[0],color=colors(f),s=25,marker='o',edgecolors='black',linewidths=.5)
                    ax.scatter(*P[-1],color=colors(f),s=30,marker='s',edgecolors='black',linewidths=.5)
                equal_axes(ax,half,view)
                ax.set_title(f'{key.upper()} · {len(c["xyz"])} Cα · azimuth {view[1]}°',fontsize=10)
                if col==0: ax.legend(loc='upper left',fontsize=7,framealpha=.85)
        fig.suptitle('Prepared continuous fragments — common views and common Å scale\nCentroid translation for display only; circle = start, square = end; colors separate fragments',fontsize=12)
        fig.savefig(directory/f'prepared_fragments_{page//5+1}.png',dpi=130)
        plt.show()

def plot_simplified(chains,paths,delta,directory):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    selected=['1o7j.a','1d9q.d','4eca.c','3jq0.a']
    half=max(np.max(np.abs(chains[key]['xyz']-chains[key]['xyz'].mean(0))) for key in selected)*1.07
    fig=plt.figure(figsize=(12,18))
    fig.subplots_adjust(left=.03,right=.97,bottom=.04,top=.94,wspace=.12,hspace=.38)
    for row,key in enumerate(selected):
        c=chains[key]; center=c['xyz'].mean(0)
        count=0
        for col,view in enumerate(VIEWS):
            ax=fig.add_subplot(4,2,row*2+col+1,projection='3d')
            count=0
            for f in c['fragments']:
                ids=paths[f'{f["id"]}@{delta:g}']['fragment_indices']
                P=f['xyz']-center; Q=P[ids]; count+=len(ids)
                ax.plot(*P.T,color='#8a949d',lw=.9,alpha=.7,label='Original' if f is c['fragments'][0] else None)
                ax.plot(*Q.T,color='#d05a2a',lw=1.7,marker='o',markersize=2.5,label='Global simplification' if f is c['fragments'][0] else None)
            equal_axes(ax,half,view)
            ax.set_title(f'{key.upper()} · {len(c["xyz"])} → {count} vertices across {len(c["fragments"])} fragment(s)',fontsize=10)
            if col==0: ax.legend(fontsize=7,loc='upper left')
    fig.suptitle(f'Independent global continuous Fréchet simplification · δ = {delta:g} Å\nNo line is drawn across a fragment boundary',fontsize=12)
    fig.savefig(directory/'independent_overlays.png',dpi=140)
    plt.show()

def plot_metrics(results,directory):
    directory=Path(directory)
    summary=results.groupby(['chain','delta_A'],sort=False)[['n','k']].sum().reset_index()
    pivot=summary.assign(retained=lambda t:t.k/t.n).pivot(index='chain',columns='delta_A',values='retained')
    fig,axes=plt.subplots(1,3,figsize=(17,6),layout='constrained')
    im=axes[0].imshow(pivot.to_numpy(),vmin=0,vmax=1,cmap='viridis',aspect='auto')
    axes[0].set(yticks=range(len(pivot)),yticklabels=pivot.index,xticks=range(len(pivot.columns)),xticklabels=[f'{x:g} Å' for x in pivot.columns],title='Retained fraction (sum over fragments)',xlabel='Continuous fidelity threshold')
    for i in range(len(pivot)):
        for j in range(len(pivot.columns)):
            x=pivot.iloc[i,j]; axes[0].text(j,i,f'{x:.0%}',ha='center',va='center',fontsize=8,color='white' if x<.55 else 'black')
    fig.colorbar(im,ax=axes[0],shrink=.7)
    for delta,g in results.groupby('delta_A'):
        axes[1].scatter(g.continuous_F_upper_A,g.discrete_F_A,label=f'δ={delta:g} Å',s=25,alpha=.8)
        axes[2].scatter(g.n,g.runtime_s,label=f'δ={delta:g} Å',s=25,alpha=.8)
    lim=max(results.discrete_F_A.max(),results.continuous_F_upper_A.max())*1.03
    axes[1].plot([0,lim],[0,lim],'--',color='gray',lw=1)
    axes[1].set(xlabel='Continuous Fréchet upper bracket (Å)',ylabel='Discrete Fréchet (Å)',title='Different fidelity measures',xlim=(0,lim),ylim=(0,lim))
    axes[1].legend(fontsize=8)
    axes[2].set(xlabel='Fragment vertices',ylabel='Solver wall time (seconds)',title='Observed runtime (first call may include JIT)')
    axes[2].grid(alpha=.2)
    fig.savefig(directory/'independent_metrics.png',dpi=150)
    plt.show()
    return summary
