from pathlib import Path
import html
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

METHODS=['GCS','CPS-3F','CPS-2F']
COLORS=['#2878B5','#D76A25','#24866C']
MARKERS=['o','s','^']

def render_table(frame,out):
    cols=['pair','tables','w','alpha','delta1','delta2','delta3','method','nA','nB','kA','kB','k',
          'auxiliary_A','auxiliary_B','auxiliary_total','extended_A','extended_B',
          'removed_pct_A','removed_pct_B','compression_factor_A','compression_factor_B',
          'contact_F1_A','contact_F1_B','contact_F1_pooled',
          'lenratio_A','lenratio_B','dDF_original','dDF_simplified','coupling_max_link','pair_bound_verified',
          'continuous_fidelity_A_upper','discrete_fidelity_A','vertices_over_GCS','seconds','status']
    names={'nA':'|A|','nB':'|B|','kA':"|A′|",'kB':"|B′|",'extended_A':'|A*|','extended_B':'|B*|',
           'auxiliary_A':'Aux added A','auxiliary_B':'Aux added B','auxiliary_total':'Aux added total',
           'delta1':'δ1 (Å)','delta2':'δ2 (Å)','delta3':'δ3 (Å)'}
    table=frame[cols].rename(columns=names).to_html(index=False,na_rep='—',float_format=lambda x:f'{x:.3f}',escape=True)
    document='''<!doctype html><meta charset="utf-8"><title>GCS / CPS comparison</title>
<style>body{font:14px system-ui;margin:24px;color:#193041}h1{font-size:24px}table{border-collapse:collapse;font-size:12px}th,td{padding:7px;border-bottom:1px solid #dae2e6;text-align:right;white-space:nowrap}th{position:sticky;top:0;background:#eaf1f5}tr:nth-child(even){background:#f6f8fa}.scroll{overflow:auto;max-height:680px}</style>
<h1>Global independent GCS · CPS-3F · thesis CPS-2F</h1>
<p>Å throughout. Pairing is discrete Fréchet, allowing repeated coupling indices and unequal output sizes. All methods retain the prepared input endpoints.</p>
<p>w ≥ 8 is a coarse-trace/runtime experiment. Blank solver metrics indicate a recorded resource limit or infeasibility, never a fabricated solution.</p><div class="scroll">'''+table+'</div>'
    Path(out).write_text(document,encoding='utf-8');return document

def summary_plots(frame,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for w,group in frame.groupby('w'):
        pairs=list(group.B.unique());alphas=list(group.alpha.unique())
        fig,axes=plt.subplots(len(alphas),2,figsize=(16,4*len(alphas)),squeeze=False)
        x=np.arange(len(pairs))
        for r,alpha in enumerate(alphas):
            part=group[group.alpha==alpha]
            for method,c,marker,offset in zip(METHODS,COLORS,MARKERS,[-.23,0,.23]):
                s=part[part.method==method].set_index('B').reindex(pairs)
                axes[r,0].bar(x+offset,s.k,width=.22,color=c,label=method)
                axes[r,1].plot(x,s.dDF_simplified/s.delta3,color=c,marker=marker,label=method,lw=1.3)
            axes[r,1].axhline(1,color='black',ls='--',lw=1,label='pair bound')
            axes[r,0].set_ylabel(f'α={alpha}: max(kA,kB)');axes[r,1].set_ylabel('d_dF(simplified) / δ₃')
            for ax in axes[r]:
                ax.set_xticks(x,pairs,rotation=40,ha='right');ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
        axes[0,0].legend(ncol=3);axes[0,1].legend(ncol=4)
        label='COARSE TRACE / RUNTIME ONLY' if w>=8 else ('UNDECIMATED CONTROL' if w==1 else 'HEADLINE')
        fig.suptitle(f'All pairs · w={w} · {label}\nVertex cost and the actual discrete coupling constraint',fontsize=16)
        fig.tight_layout(rect=(0,0,1,.94));fig.savefig(out/f'comparison_w{w}.png',dpi=150,bbox_inches='tight');plt.show()
        fig,ax=plt.subplots(figsize=(12,5))
        for method,c,marker in zip(METHODS,COLORS,MARKERS):
            s=group[group.method==method]
            ax.scatter(s.delta1,s.seconds,color=c,marker=marker,label=method,alpha=.7)
        ax.set(xlabel='Fidelity threshold δ₁ = δ₂ (Å)',ylabel='Elapsed seconds (log scale)',yscale='log',title=f'Solver time · w={w} · {label}')
        ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/f'runtime_w{w}.png',dpi=150,bbox_inches='tight');plt.show()

def spatial_plots(data,paths,out,alpha=1,w=16):
    """One row per pair; identical view and spatial limits across its four panels."""
    out=Path(out);items=[x for x in data['instances'].values() if x['w']==w]
    for page in range(0,len(items),4):
        subset=items[page:page+4];fig=plt.figure(figsize=(17,4.3*len(subset)))
        for r,item in enumerate(subset):
            A=item['A'];B=item['B_xyz'];both=np.vstack([A,B]);center=(both.max(0)+both.min(0))/2
            from curve_algorithms import discrete_frechet
            d=alpha*.5*(np.linalg.norm(np.diff(A,axis=0),axis=1).mean()+np.linalg.norm(np.diff(B,axis=0),axis=1).mean())
            d3=np.ceil(discrete_frechet(A,B))
            radius=max(np.ptp(both,axis=0).max()/2,1)*1.08
            for col,method in enumerate(['Input']+METHODS):
                ax=fig.add_subplot(len(subset),4,4*r+col+1,projection='3d')
                a=np.arange(len(A));b=np.arange(len(B));coupling=None
                if method!='Input':
                    result=paths.get(f"{item['instance']}@{alpha}@{method}")
                    if result:
                        a=np.array(result['indices_A']);b=np.array(result['indices_B']);coupling=result['coupling']
                    else:
                        ax.text2D(.1,.5,'No solver result\nSee status table',transform=ax.transAxes)
                        a=np.array([],int);b=np.array([],int)
                ax.plot(*A.T,color=COLORS[0],alpha=.15,lw=.8);ax.plot(*B.T,color=COLORS[1],alpha=.15,lw=.8)
                if len(a):ax.plot(*A[a].T,color=COLORS[0],marker='o',ms=3,lw=1.7)
                if len(b):ax.plot(*B[b].T,color=COLORS[1],marker='s',ms=3,lw=1.7)
                if coupling is not None:
                    for u,v in coupling:
                        ax.plot(*np.vstack([A[a[u]],B[b[v]]]).T,color='#52646d',alpha=.55,ls=':',lw=.7)
                for setter,k in [(ax.set_xlim,0),(ax.set_ylim,1),(ax.set_zlim,2)]:setter(center[k]-radius,center[k]+radius)
                ax.set_box_aspect((1,1,1));ax.view_init(elev=22,azim=-58)
                ax.set_xlabel('x (Å)',labelpad=-1);ax.set_ylabel('y (Å)',labelpad=-1);ax.set_zlabel('z (Å)',labelpad=-1)
                ax.tick_params(labelsize=7,pad=0)
                ax.set_title(f"{item['pair'].replace('__',' / ')} · {method}\n"
                             f"δ1={d:.3f}, δ2={d:.3f}, δ3={d3:.3f} Å\n"
                             f"|A|={len(A)}, |B|={len(B)}; output={len(a)} / {len(b)}",fontsize=8)
        fig.suptitle(f'Prepared input and all three simplifications · α={alpha}, w={w}\nBlue circles: A (1o7j.a); orange squares: B; dotted lines: discrete coupling\n'+('COARSE TRACE / RUNTIME ONLY' if w>=8 else 'HEADLINE / CONTROL'),fontsize=14)
        fig.subplots_adjust(top=.90,bottom=.03,left=.02,right=.98,hspace=.27,wspace=.10)
        fig.savefig(out/f'geometry_w{w}_alpha{alpha}_page{page//4+1}.png',dpi=140,bbox_inches='tight');plt.show()
