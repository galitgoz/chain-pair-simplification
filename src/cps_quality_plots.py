"""Every comparison has a pair label and all three tolerance values."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def academic_dual_bars(results,out):
    """One publication figure per (stride, alpha); no averaging across budgets."""
    from matplotlib.ticker import MaxNLocator
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    methods=[('GCS','Independent','#4477AA',''),('CPS-2F','CPS-2F','#228833','//'),
             ('CPS-3F','CPS-3F','#EE7733','..')]
    tables=[]
    style={'font.family':'serif','font.serif':['STIXGeneral','DejaVu Serif'],
           'mathtext.fontset':'stix','font.size':10,'axes.titlesize':12,
           'axes.labelsize':12,'xtick.labelsize':8.5,'ytick.labelsize':9,
           'axes.linewidth':.7,'pdf.fonttype':42,'ps.fonttype':42,
           'hatch.linewidth':.45,'savefig.facecolor':'white'}
    with plt.rc_context(style):
        for (w,alpha),group in results.groupby(['w','alpha'],sort=True):
            pairs=list(group.pair.unique());x=np.arange(len(pairs));width=.24
            params=group.drop_duplicates('pair')[['pair','w','alpha','delta1','delta2','delta3']].copy()
            name=f'dual_bars_w{w}_alpha{alpha:g}'
            params['figure']=name;tables.append(params)
            fig,axes=plt.subplots(1,2,figsize=(11.5,5.7),sharex=True)
            for offset,(method,label,color,hatch) in zip([-width,0,width],methods):
                rows=group[group.method==method].set_index('pair').reindex(pairs)
                solved=rows.status.str.startswith('optimal')
                coupling=(rows.dDF_simplified/rows.delta3).where(solved & rows.delta3.gt(0))
                k=rows.k.where(solved)
                for ax,values in zip(axes,[coupling,k]):
                    ax.bar(x+offset,values,width=width*.94,label=label,color=color,
                           edgecolor='white',linewidth=.4,hatch=hatch,zorder=3)
                    missing=~np.isfinite(values.to_numpy())
                    for pos in x[missing]+offset:
                        ax.text(pos,.018,'NA',transform=ax.get_xaxis_transform(),
                                ha='center',va='bottom',rotation=90,fontsize=6.5,color=color)
            axes[0].axhline(1.,color='#222222',linestyle=(0,(5,3)),linewidth=1.7,zorder=5)
            axes[0].set_title('Joint solution stays inside the budget',pad=13)
            axes[0].set_ylabel(r"$d_{\mathrm{dF}}(A',B')\,/\,\delta_3$ budget")
            axes[1].set_title('at little or no cost in vertices',pad=13)
            axes[1].set_ylabel(r"$k = \max(|A'|,\,|B'|)$")
            axes[1].legend(loc='upper center',bbox_to_anchor=(.5,1.015),ncol=3,
                           frameon=False,fontsize=9,handlelength=1.4,columnspacing=1.2)
            axes[1].yaxis.set_major_locator(MaxNLocator(integer=True))
            for ax in axes:
                ax.set_xticks(x,[p.replace('__',' / ') for p in pairs],rotation=90,ha='center')
                ax.set_xlim(-.6,len(pairs)-.4)
                ax.set_ylim(0,ax.get_ylim()[1]*1.18)
                ax.spines[['top','right']].set_visible(False)
                ax.grid(axis='y',color='#dddddd',linewidth=.55,zorder=0)
                ax.tick_params(axis='x',length=0,pad=5)
            fig.text(.5,.015,rf'$w={w},\ \alpha={alpha:g}$'
                     +' · Thresholds vary by pair; see accompanying parameter table. NA = unsolved.',
                     ha='center',fontsize=9)
            fig.subplots_adjust(left=.07,right=.99,bottom=.32,top=.87,wspace=.26)
            for extension in ('pdf','svg','png'):
                fig.savefig(out/f'{name}.{extension}',dpi=300,bbox_inches='tight',pad_inches=.04)
            plt.show();plt.close(fig)
    table=pd.concat(tables,ignore_index=True)
    table.to_csv(out.parent/'dual_bar_parameters.csv',index=False)
    return table


def budget_vertex_boxplots(results,out):
    """Matched-case boxplots, plus explicit paired overhead over independent GCS."""
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    methods=['GCS','CPS-3F','CPS-2F'];colors=['#2878B5','#D76A25','#24866C']
    cohorts=[]
    for w,group in results.groupby('w'):
        valid=group[group.status.str.startswith('optimal') & group.delta3.gt(0)].copy()
        counts=valid.groupby('case').method.nunique()
        keep=counts[counts==len(methods)].index
        cohort=group.drop_duplicates('case')[['case','pair','w','alpha','delta1','delta2','delta3']].copy()
        cohort['included']=cohort.case.isin(keep)
        cohort['exclusion_reason']=np.where(cohort.included,'','At least one method unsolved or δ3 not positive')
        cohorts.append(cohort)
        sample=valid[valid.case.isin(keep)].copy()
        if sample.empty:continue
        baseline=sample[sample.method=='GCS'].set_index('case').k
        sample['k_independent']=sample.case.map(baseline)
        sample['coupling_budget_fraction']=sample.dDF_simplified/sample.delta3
        sample['extra_vertices']=sample.k-sample.k_independent
        sample['vertex_cost_ratio']=sample.k/sample.k_independent
        sample.to_csv(out.parent/f'boxplot_values_w{w}.csv',index=False)
        assert sample.groupby('method').size().nunique()==1
        assert sample[sample.method!='GCS'].coupling_budget_fraction.le(1+1e-8).all()
        assert sample.extra_vertices.ge(0).all()

        def draw(ax,metric,ylabel,title,reference=None):
            values=[sample[sample.method==m][metric].to_numpy() for m in methods]
            boxes=ax.boxplot(values,tick_labels=['Independent\nGCS','CPS-3F','CPS-2F'],
                patch_artist=True,widths=.52,showfliers=False,
                medianprops={'color':'#17252a','linewidth':2})
            for patch,c in zip(boxes['boxes'],colors):patch.set(facecolor=c,alpha=.28)
            # Fixed jitter; each dot is one matched (pair, stride, alpha) case.
            rng=np.random.default_rng(2457)
            for i,(values,c) in enumerate(zip(values,colors),1):
                ax.scatter(i+rng.uniform(-.15,.15,len(values)),values,s=22,color=c,alpha=.65,zorder=3)
            if reference is not None:ax.axhline(reference,color='#444444',ls='--',lw=1)
            ax.set_ylabel(ylabel);ax.set_title(title,loc='left',fontweight='bold',pad=12)
            ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
        note=(f'w={w} · {len(keep)} matched cases from {sample.pair.nunique()} protein pairs · '
              f'{len(cohort)-len(keep)} cases excluded from every method\n'
              'Boxes: median and IQR; whiskers: 1.5 IQR; dots: individual cases. '
              'Exact pair and δ1/δ2/δ3: accompanying cohort table.')
        fig,axes=plt.subplots(1,2,figsize=(12.5,5.8))
        draw(axes[0],'coupling_budget_fraction','d_dF(A′, B′) / δ3',
             '(a) Achieved coupling / pairing budget',1)
        axes[0].text(.98,1.015,'Budget = 1',transform=axes[0].get_yaxis_transform(),ha='right',fontsize=9)
        draw(axes[1],'k','k = max(|A′|, |B′|)', '(b) Vertex cost by method')
        axes[1].set_ylim(bottom=0)
        fig.suptitle('Coupling and vertex cost against independent GCS',fontsize=15)
        fig.text(.5,.02,note,ha='center',fontsize=9)
        fig.tight_layout(rect=(0,.13,1,.94))
        fig.savefig(out/f'coupling_vertex_boxplots_w{w}.png',dpi=180,bbox_inches='tight')
        fig.savefig(out/f'coupling_vertex_boxplots_w{w}.pdf',bbox_inches='tight')
        plt.show();plt.close(fig)
        fig,ax=plt.subplots(figsize=(8,5))
        draw(ax,'extra_vertices','k_method − k_independent','Paired vertex overhead relative to independent GCS',0)
        fig.text(.5,.015,f'Same {len(keep)} matched cases · w={w} · 0 means no additional vertices',ha='center',fontsize=9)
        fig.tight_layout(rect=(0,.06,1,1))
        fig.savefig(out/f'vertex_overhead_boxplots_w{w}.png',dpi=180,bbox_inches='tight')
        plt.show();plt.close(fig)
    cohort=pd.concat(cohorts,ignore_index=True)
    cohort.to_csv(out.parent/'boxplot_cohort.csv',index=False)
    return cohort


def quality_plots(results,out):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    metrics=[('removed_pct_total','Vertices removed (%) ↑'),
             ('contact_F1_pooled','Contact preservation F1 ↑'),
             ('max_continuous_error','Max continuous Fréchet error (Å) ↓')]
    for w,group in results.groupby('w'):
        pairs=list(group.pair.unique())
        for page in range(0,len(pairs),3):
            chunk=pairs[page:page+3]
            fig,axes=plt.subplots(len(chunk),3,figsize=(18,4.5*len(chunk)),squeeze=False)
            for r,pair in enumerate(chunk):
                rows=group[group.pair==pair];cases=rows.drop_duplicates('case').sort_values('alpha')
                x=np.arange(len(cases))
                labels=[f"α={v.alpha:g}\nδ1={v.delta1:.3f}\nδ2={v.delta2:.3f}\nδ3={v.delta3:.3f}" for v in cases.itertuples()]
                for c,(metric,label) in enumerate(metrics):
                    ax=axes[r,c]
                    for method,color,marker in [('CPS-2F','#24866C','^'),('CPS-3F','#D76A25','s')]:
                        m=rows[rows.method==method].set_index('case').reindex(cases.case)
                        ax.plot(x,m[metric],marker=marker,color=color,label=method,lw=1.5,ms=7)
                        for pos,(_,v) in enumerate(m.iterrows()):
                            if not str(v.status).startswith('optimal'):
                                ax.text(pos,.05,f'{method}\nlimit',transform=ax.get_xaxis_transform(),ha='center',fontsize=8,color=color)
                            elif not np.isfinite(v[metric]):
                                ax.text(pos,.05,'F1 undefined',transform=ax.get_xaxis_transform(),ha='center',fontsize=8,color=color)
                    ax.set_title(pair.replace('__',' / ')+f' · w={w}',fontsize=10)
                    ax.set_xticks(x,labels,fontsize=8);ax.set_ylabel(label);ax.grid(alpha=.2)
                    if c==0:ax.set_ylim(-3,103)
                    if c==1:ax.set_ylim(-.04,1.04)
                    ax.legend(fontsize=8)
            fig.suptitle('CPS-2F versus CPS-3F · identical pair and tolerances in each comparison\n'
                         'Contact F1 is a structural proxy, not a knot invariant. All δ values are in Å.',fontsize=13)
            fig.tight_layout(rect=(0,0,1,.95))
            fig.savefig(out/f'quality_w{w}_page{page//3+1}.png',dpi=140,bbox_inches='tight')
            plt.show();plt.close(fig)
