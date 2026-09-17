"""Hurricane data, exact CPS comparisons, and domain-specific reporting."""
from pathlib import Path
import json, time, hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from curve_algorithms import shortest_independent, discrete_frechet, continuous_distance, continuous_decision
from cps_paper_algorithms import solve_cps, ResourceLimit, discrete_coupling

METHODS = ['GCS','CPS-2F','CPS-3F']
COLORS = ['#4477AA','#228833','#EE7733']

def load_hurricanes(folder):
    folder = Path(folder)
    tracks = {p.stem:np.ascontiguousarray(np.loadtxt(p),dtype=float) for p in sorted(folder.glob('*.txt'))}
    metadata = pd.read_csv(folder/'packedness.csv')
    pairs = pd.read_csv(folder/'pairs.csv')
    assert set(metadata.id) == set(tracks)
    for r in metadata.itertuples():
        P = tracks[r.id]
        assert P.shape == (r.n,2) and np.isfinite(P).all()
        assert (np.linalg.norm(np.diff(P,axis=0),axis=1)>0).all()
    assert set(pairs.curve_A)|set(pairs.curve_B) <= set(tracks)
    info = metadata.set_index('id')
    assert all(info.loc[r.curve_A,'basin']==info.loc[r.curve_B,'basin']==r.basin for r in pairs.itertuples())
    pairs['supplied_distance_order_violation'] = pairs.frechet_km > pairs.discrete_frechet_km + .02
    return tracks, metadata, pairs


def organize(tracks, pairs, pairs_per_basin=4):
    """Deterministic ranks across supplied ball-packedness, not solver outcomes."""
    chosen=[]
    for basin,g in pairs.groupby('basin',sort=True):
        g=g.sort_values(['c_ball_pair','curve_A','curve_B'])
        ids=np.arange(len(g)) if pairs_per_basin is None else np.unique(np.linspace(0,len(g)-1,min(pairs_per_basin,len(g))).round().astype(int))
        chosen.append(g.iloc[ids])
    selected=pd.concat(chosen,ignore_index=True)
    instances={}
    for r in selected.itertuples():
        pair=f'{r.curve_A}__{r.curve_B}'
        instances[f'{pair}@1']=dict(pair=pair,w=1,A=tracks[r.curve_A],B_xyz=tracks[r.curve_B],
            basin=r.basin,name_A=r.name_A,name_B=r.name_B)
    return dict(instances=instances),selected


def run(data, selected, alphas, budgets, out, resume=True):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rows=[]; paths={}; audit=[]
    cached={};cached_paths={}
    if resume and all((out/name).exists() for name in ['comparison.csv','paths.json','manifest.json']):
        manifest=json.loads((out/'manifest.json').read_text())
        root=Path(__file__).resolve().parent
        # Reuse only exact optimal results from unchanged algorithms and inputs.
        relevant={name:digest for name,digest in manifest['sha256'].items()
                  if name.startswith('data/hurricanes/') or name in ['curve_algorithms.py','cps_paper_algorithms.py']}
        assert relevant and all(hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
                                for name,digest in relevant.items()),'Cached inputs/algorithms changed; use resume=False.'
        previous=pd.read_csv(out/'comparison.csv')
        cached={(r['case'],r['method']):r for r in previous.to_dict('records') if r['status'].startswith('optimal')}
        cached_paths=json.loads((out/'paths.json').read_text())
    for key,item in data['instances'].items():
        A,B=item['A'],item['B_xyz']
        original=float(discrete_frechet(A,B)); clo,chi=continuous_distance(A,B)
        assert clo<=original+1e-5
        source=selected[(selected.curve_A+'__'+selected.curve_B).eq(item['pair'])].iloc[0]
        audit.append(dict(pair=item['pair'],supplied_continuous=source.frechet_km,
            supplied_discrete=source.discrete_frechet_km,recomputed_continuous_lower=clo,
            recomputed_continuous_upper=chi,recomputed_discrete=original,
            supplied_order_violation=source.supplied_distance_order_violation))
        edge=.5*(np.linalg.norm(np.diff(A,axis=0),axis=1).mean()+np.linalg.norm(np.diff(B,axis=0),axis=1).mean())
        for alpha in alphas:
            d=float(alpha*edge);d3=float(np.ceil(original));case=f'{key}@{alpha}'
            start=time.perf_counter();ia=shortest_independent(A,d);ib=shortest_independent(B,d)
            independent=dict(indices_A=ia,indices_B=ib,status='optimal_global_GCS',seconds=time.perf_counter()-start)
            for method in METHODS:
                prior=cached.get((case,method))
                if prior is not None and f'{case}@{method}' in cached_paths:
                    assert np.allclose([prior['delta1'],prior['delta2'],prior['delta3']],[d,d,d3],rtol=0,atol=1e-8)
                    assert (prior['nA'],prior['nB'])==(len(A),len(B))
                    prior=prior.copy();prior['reused_verified_result']=True
                    rows.append(prior);paths[f'{case}@{method}']=cached_paths[f'{case}@{method}']
                    continue
                start=time.perf_counter()
                if method=='GCS': result=independent
                else:
                    try: result=solve_cps(A,B,d,d,d3,method,**budgets)
                    except ResourceLimit as e: result=dict(status='resource_limit',reason=str(e),seconds=time.perf_counter()-start)
                row=dict(pair=item['pair'],case=case,w=1,alpha=alpha,method=method,basin=item['basin'],
                    nA=len(A),nB=len(B),delta1=d,delta2=d,delta3=d3,dDF_original=original,
                    c_ball_pair=source.c_ball_pair,c_cube_pair=source.c_cube_pair,
                    status=result['status'],seconds=result.get('seconds',time.perf_counter()-start),
                    reason=result.get('reason',''),states=result.get('states'),transitions=result.get('transitions'),
                    reused_verified_result=False,soft_time_limit_seconds=budgets['seconds'],
                    transition_limit=budgets['max_transitions'])
                if 'indices_A' in result:
                    a,b=result['indices_A'],result['indices_B'];SA,SB=A[a],B[b]
                    coupling,dist=discrete_coupling(SA,SB)
                    ca=continuous_distance(A,SA);cb=continuous_distance(B,SB)
                    assert continuous_decision(A,SA,d) and continuous_decision(B,SB,d)
                    if method!='GCS': assert dist<=d3+1e-9
                    if method=='CPS-3F': assert discrete_frechet(A,SA)<=d+1e-9 and discrete_frechet(B,SB)<=d+1e-9
                    row.update(kA=len(a),kB=len(b),k=max(len(a),len(b)),
                        budget_ratio=dist/d3,dDF_simplified=dist,
                        removed_pct=100*(1-(len(a)+len(b))/(len(A)+len(B))),
                        fidelity_A_upper=ca[1],fidelity_B_upper=cb[1],
                        fidelity_ratio=max(ca[1]/d,cb[1]/d),
                        vertices_over_GCS=max(len(a),len(b))-max(len(ia),len(ib)))
                    paths[f'{case}@{method}']=dict(indices_A=a.tolist(),indices_B=b.tolist(),coupling=coupling.tolist())
                rows.append(row)
            print(f"{item['pair']} alpha={alpha}: "+', '.join(f"{r['method']}={r.get('k',r['status'])}" for r in rows[-3:]),flush=True)
            pd.DataFrame(rows).to_csv(out/'comparison.csv',index=False)
            (out/'paths.json').write_text(json.dumps(paths,indent=2))
    pd.DataFrame(audit).to_csv(out/'distance_audit.csv',index=False)
    return pd.DataFrame(rows),paths,pd.DataFrame(audit)


def save_figure(fig,out,name):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for ext in ['pdf','svg','png']:fig.savefig(out/f'{name}.{ext}',dpi=220,bbox_inches='tight')
    plt.show();plt.close(fig)


def comparison_plots(results,out,page_size=8,make_boxplot=True):
    paged=page_size is not None and results.pair.nunique()>page_size
    if paged:
        pairs=results.pair.drop_duplicates().tolist()
        for start in range(0,len(pairs),page_size):
            comparison_plots(results[results.pair.isin(pairs[start:start+page_size])],
                             Path(out)/f'page_{start//page_size+1:02d}',page_size=None,make_boxplot=False)
    with plt.rc_context({'font.family':'serif','mathtext.fontset':'stix','pdf.fonttype':42}):
        for alpha,g in ([] if paged else results.groupby('alpha')):
            pairs=g.pair.drop_duplicates().tolist(); x=np.arange(len(pairs)); width=.25
            fig,axes=plt.subplots(1,2,figsize=(14,7))
            for mi,method in enumerate(METHODS):
                r=g[g.method.eq(method)].set_index('pair').reindex(pairs)
                for ax,field in zip(axes,['budget_ratio','k']):
                    ax.bar(x+(mi-1)*width,r[field],width,color=COLORS[mi],label='Independent' if method=='GCS' else method)
                    for xx,val in zip(x+(mi-1)*width,r[field]):
                        if pd.isna(val):ax.text(xx,.03,'NA',rotation=90,transform=ax.get_xaxis_transform(),fontsize=7)
            axes[0].axhline(1,color='black',ls='--',lw=1.5)
            axes[0].set_title('Joint solution stays inside the budget')
            axes[0].set_ylabel(r'$d_{dF}(A^\prime,B^\prime)/\delta_3$ budget')
            axes[1].set_title('Vertex cost of the joint constraint')
            axes[1].set_ylabel(r'$k=\max(|A^\prime|,|B^\prime|)$');axes[1].legend(frameon=False)
            for ax in axes:
                ax.set_xticks(x,[p.replace('__',' / ') for p in pairs],rotation=90)
                ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False)
            fig.suptitle(f'Hurricane tracks, full supplied resolution; α={alpha:g}')
            fig.text(.5,.015,'Pair IDs identify storms and years. Exact δ₁/δ₂/δ₃ in km: parameters.csv. NA = solver limit.',ha='center')
            fig.subplots_adjust(bottom=.37,wspace=.25)
            save_figure(fig,out,f'dual_bars_alpha{alpha:g}')
        if not make_boxplot:return
        pivot=results.pivot(index='case',columns='method',values='k')
        complete=pivot.dropna().index
        cohort=results[results.case.isin(complete)]
        fig,axes=plt.subplots(1,3,figsize=(13,4))
        for ax,field,label in zip(axes,['budget_ratio','k','vertices_over_GCS'],
            [r'$d_{dF}(A^\prime,B^\prime)/\delta_3$',r'$k=\max(|A^\prime|,|B^\prime|)$',r'$k-k_{Independent}$']):
            values=[cohort[cohort.method.eq(m)][field].dropna().to_numpy() for m in METHODS]
            ax.boxplot(values,tick_labels=['Independent','CPS-2F','CPS-3F'])
            for i,v in enumerate(values):ax.scatter(np.full(len(v),i+1),v,s=15,color=COLORS[i],alpha=.7)
            ax.set_ylabel(label);ax.grid(axis='y',alpha=.2)
        axes[0].axhline(1,color='black',ls='--');axes[2].axhline(0,color='black',ls='--')
        fig.suptitle(f'Matched solved cases only: {len(complete)} pair–threshold cases')
        fig.tight_layout();save_figure(fig,out,'matched_boxplots')
        cohort[['case','pair','alpha','delta1','delta2','delta3']].drop_duplicates().to_csv(Path(out).parent/'boxplot_cohort.csv',index=False)


def track_plots(data,results,paths,out,alpha=1):
    from matplotlib.ticker import MaxNLocator
    for key,item in data['instances'].items():
        A,B=item['A'],item['B_xyz'];case=f'{key}@{alpha}'
        g=results[results.case.eq(case)].set_index('method');r=g.iloc[0]
        fig,axes=plt.subplots(1,3,figsize=(15,6),sharex=True,sharey=True)
        combined=np.vstack([A,B]);center=.5*(combined.min(axis=0)+combined.max(axis=0))
        half=.55*np.ptp(combined,axis=0).max()
        for ax,method in zip(axes,METHODS):
            ax.plot(A[:,0],A[:,1],color='#4477AA',alpha=.25,lw=2)
            ax.plot(B[:,0],B[:,1],color='#EE7733',alpha=.25,lw=2)
            sol=paths.get(f'{case}@{method}')
            if sol:
                SA=A[sol['indices_A']];SB=B[sol['indices_B']]
                ax.plot(SA[:,0],SA[:,1],'-o',color='#4477AA',ms=4,label='A simplified')
                ax.plot(SB[:,0],SB[:,1],'-s',color='#EE7733',ms=4,label='B simplified')
                for i,j in sol['coupling']:ax.plot([SA[i,0],SB[j,0]],[SA[i,1],SB[j,1]],':',color='gray',alpha=.5,lw=.7)
                detail=f"|A′|={len(SA)}, |B′|={len(SB)}; coupling/budget={g.loc[method,'budget_ratio']:.2f}"
            else:detail='Resource limit; no solution reported'
            ax.scatter(*A[0],marker='*',s=75,color='#4477AA');ax.scatter(*B[0],marker='*',s=75,color='#EE7733')
            ax.set_title(('Independent GCS' if method=='GCS' else method)+'\n'+detail,fontsize=10)
            ax.set_xlim(center[0]-half,center[0]+half);ax.set_ylim(center[1]-half,center[1]+half)
            ax.set_aspect('equal',adjustable='box');ax.set_xlabel('Projected easting (km)');ax.grid(alpha=.2)
            ax.tick_params(labelsize=9)
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(6))
        axes[0].set_ylabel('Projected northing (km)')
        handles,labels=axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.90),ncol=2,frameon=False,fontsize=9)
        fig.suptitle(f"{item['pair'].replace('__',' / ')} — {item['name_A']} / {item['name_B']} ({item['basin']})\n"
            f"|A|={len(A)}, |B|={len(B)}; δ₁={r.delta1:.2f}, δ₂={r.delta2:.2f}, δ₃={r.delta3:.2f} km; α={alpha}")
        fig.text(.5,-.055,'Faint lines: input tracks. Dotted links: discrete coupling. Stars: start points. Same coordinates and scales across methods.',ha='center',fontsize=9)
        fig.tight_layout(rect=[0,.06,1,.83]);save_figure(fig,out,f'tracks_{item["pair"]}_alpha{alpha}')
