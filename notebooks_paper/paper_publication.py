"""Publication reporting from saved evidence only; never runs an optimizer.

Run after execution is stopped/completed. Historical artifacts are read-only.
"""
from pathlib import Path
import json, hashlib, shutil, re
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'output/paper_v1/expanded_comparison'
TIMING = SOURCE/'timing_study'
OUT = ROOT/'output/paper_v1/experimental_chapter'
METHODS = ['Independent continuous','Independent discrete','CPS-2F','CPS-3F']
SHORT = dict(zip(METHODS,['IC','ID','CPS-2F','CPS-3F']))
COLORS = dict(zip(METHODS,['#377eb8','#984ea3','#e66101','#1b9e77']))
CAPTIONS = {}

def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
    from repository_paths import resolve_recorded_path
    return hashlib.sha256(resolve_recorded_path(p).read_bytes()).hexdigest()
def dump(p,x): Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def tex(x):
    s=str(x)
    return ''.join({'_':r'\_','%':r'\%','&':r'\&','#':r'\#'}.get(c,c) for c in s)
def num(x,d=3): return '--' if pd.isna(x) else f'{float(x):.{d}g}'
def savefig(fig,name,caption):
    CAPTIONS[name]=caption
    fig.savefig(OUT/'figures'/f'{name}.pdf',bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/'figures'/f'{name}.png',dpi=180,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def method_legend(ax):
    ax.legend(handles=[Line2D([],[],color=COLORS[m],marker='o',ls='',label=m) for m in METHODS],fontsize=8)
def latex_table(frame,path,caption,label):
    # Generate directly from data without adding Jinja2 to the frozen environment.
    header=' & '.join(tex(c) for c in frame.columns)+r' \\'
    lines=[r'\begin{longtable}{'+'l'*len(frame.columns)+'}',r'\caption{'+tex(caption)+r'}\label{'+label+r'}\\',r'\toprule',header,r'\midrule',r'\endfirsthead',r'\toprule',header,r'\midrule',r'\endhead']
    for row in frame.itertuples(index=False,name=None):lines.append(' & '.join('--' if pd.isna(x) else tex(x) for x in row)+r' \\')
    lines.extend([r'\bottomrule',r'\end{longtable}']);text='\n'.join(lines)+'\n'
    path.write_text('\\begingroup\n\\small\n'+text+'\\endgroup\n',encoding='utf-8')

def summaries():
    original=pd.read_csv(SOURCE/'results.csv'); timing=pd.read_csv(TIMING/'results.csv')
    manifest=pd.read_csv(SOURCE/'input_manifest.csv')
    aliases={}; order=[]
    for domain,prefix in [('protein','P'),('hurricane','H')]:
        for i,pair in enumerate(manifest[manifest.domain.eq(domain)].pair,1):
            aliases[pair]=f'{prefix}{i}';order.append(pair)
    for frame in [original,timing,manifest]:frame['pair_alias']=frame.pair.map(aliases)
    original.to_csv(OUT/'tables/original_cohort.csv',index=False)
    timing.to_csv(OUT/'tables/timing_observations.csv',index=False)
    rows=[]
    for cid,g in timing.groupby('base_configuration_id',sort=False):
        g=g.sort_values('repetition');v=g[g.validated.eq(True)];resources=g[g.resource_limited.eq(True)]
        row=dict(configuration_id=cid,pair=g.iloc[0].pair,pair_alias=g.iloc[0].pair_alias,method=g.iloc[0].method,
                 planned=3,recorded=int(g.origin.eq('new').sum()),successful=len(v),
                 successful_solver_median=v.solver_seconds.median(),successful_solver_min=v.solver_seconds.min(),
                 successful_solver_max=v.solver_seconds.max(),resource_exit_count=len(resources),
                 resource_exit_seconds=json.dumps(resources.solver_seconds.dropna().tolist()),
                 objectives=json.dumps(v.k.tolist()),objective_consistent=len(set(v.k))<=1,
                 peak_rss_mib_max=g.peak_rss_mib.max())
        for _,r in g.iterrows():
            for key in ['status','solver_seconds','warmup_seconds','validation_seconds','total_seconds','peak_rss_mib','k','attempt_path']:
                row[f'r{int(r.repetition)}_{key}']=r[key]
        selections=[]
        for _,r in v.iterrows():
            result=read(ROOT/r.attempt_path/'raw_solver.json')['raw_result']
            selections.append(json.dumps([result['indices_A'],result['indices_B']]))
        row['distinct_validated_output_pairs']=len(set(selections))
        rows.append(row)
    summary=pd.DataFrame(rows);summary.to_csv(OUT/'tables/timing_summary.csv',index=False)
    paired=[]
    for pair in order:
        group=original[original.pair.eq(pair)].set_index('method')
        for fidelity,base,joint in [('continuous',METHODS[0],METHODS[2]),('discrete',METHODS[1],METHODS[3])]:
            b,c=group.loc[base],group.loc[joint];valid=bool(b.validated and c.validated)
            paired.append(dict(pair=pair,pair_alias=aliases[pair],fidelity=fidelity,baseline=base,joint=joint,
                baseline_validated=b.validated,joint_validated=c.validated,baseline_k=b.k,joint_k=c.k,
                delta_k=c.k-b.k if valid else np.nan,baseline_coupling_ratio=b.coupling_distance/b.delta3,
                joint_coupling_ratio=c.coupling_distance/c.delta3 if c.validated else np.nan,
                equal_k_coupling_repair=bool(valid and not b.coupling_feasible and c.coupling_feasible and b.k==c.k),
                baseline_source=b.attempt_path,joint_source=c.attempt_path))
    paired=pd.DataFrame(paired);paired.to_csv(OUT/'tables/paired_quality.csv',index=False)
    thresholds=original.drop_duplicates('pair')[['pair','delta1','delta2','delta3']]
    inputs=manifest.merge(thresholds,on='pair');inputs.to_csv(OUT/'tables/inputs.csv',index=False)
    bounds=pd.read_csv(SOURCE/'graph_bounds.csv');bounds['pair_alias']=bounds.pair.map(aliases)
    bounds.to_csv(OUT/'tables/graph_bounds.csv',index=False)
    return original,timing,summary,paired,inputs,bounds,aliases,order

def coupling_figure(paired,aliases,order):
    fig,axes=plt.subplots(2,2,figsize=(10.6,6.4),sharex='col',layout='constrained')
    for col,fidelity in enumerate(['continuous','discrete']):
        d=paired[paired.fidelity.eq(fidelity)].set_index('pair').loc[order];x=np.arange(len(d))
        base=METHODS[col];joint=METHODS[col+2]
        ax=axes[0,col];ax.axhline(1,color='black',ls='--',lw=1.3)
        for offset,key,m in [(-.14,'baseline_coupling_ratio',base),(.14,'joint_coupling_ratio',joint)]:
            ax.scatter(x+offset,d[key],color=COLORS[m],label=m,s=28,zorder=3)
        ax.set_title(f'{fidelity.capitalize()} fidelity');ax.set_ylabel(r'$d_{dF}(A^\prime,B^\prime)/\delta_3$');ax.legend(fontsize=8)
        ax=axes[1,col];ax.axhline(0,color='.6',lw=.7)
        ax.scatter(x,d.delta_k,color=COLORS[joint],s=30)
        repair=d.equal_k_coupling_repair.to_numpy()
        ax.scatter(x[repair],d.delta_k.to_numpy()[repair],marker='*',s=125,facecolor='gold',edgecolor='black',label='Equal k; coupling repaired',zorder=4)
        for pos,missing in enumerate(d.delta_k.isna()):
            if missing:ax.text(pos,.15,'no output',ha='center',va='bottom',rotation=90,fontsize=7,color='darkred')
        ax.set_ylabel(r'$k_{\rm CPS}-k_{\rm independent}$');ax.set_xticks(x,[aliases[p] for p in order]);ax.set_ylim(-.3,max(1.4,float(d.delta_k.max())+.4 if d.delta_k.notna().any() else 1.4))
        ax.set_xlabel('Pair aliases in the saved input manifest')
        if repair.any():ax.legend(fontsize=8)
    savefig(fig,'coupling_vertex_cost','Does enforcing coupling require more vertices? Output discrete coupling divided by the shared threshold (top), and matched-fidelity objective difference (bottom). Dashed line is the coupling bound. Independent violations are permitted. Stars identify equal-objective coupling repairs. Missing comparisons are not zero cost; their outcome or pending status is in the result table. Pair aliases and exact thresholds are in the input manifest.')

def bound_figure(bounds,aliases,order):
    fig,axes=plt.subplots(2,2,figsize=(16 if len(order)>30 else 11,7),layout='constrained')
    for col,m in enumerate(METHODS[2:]):
        b=bounds[bounds.method.eq(m)];g=b[b.quantity.eq('configuration_states_input_bound')].set_index('pair').loc[order]
        a=b[b.quantity.eq('configuration_states_active_bound')].set_index('pair').loc[order]
        x=np.arange(len(g));ax=axes[0,col]
        ax.plot(x,g.finite_upper_bound,'s',color='.5',label=r'Input-only bound $B_{\rm input}$')
        ax.plot(x,a.finite_upper_bound,'^',color='#386cb0',label=r'Active product $H_AH_B$')
        done=g.search_completed.eq(True).to_numpy()
        ax.scatter(x[done],g.observed.to_numpy()[done],color=COLORS[m],label='Discovered; completed search')
        ax.scatter(x[~done],g.observed.to_numpy()[~done],color='black',marker='x',s=65,label='Observed before termination')
        ax.set_yscale('log');ax.set_title(m);ax.set_ylabel('Configuration states / finite bound');ax.legend(fontsize=7)
        ax=axes[1,col]
        for off,frame,label,color in [(-.12,g,r'$100D/B_{\rm input}$','.5'),(.12,a,r'$100D/(H_AH_B)$','#386cb0')]:
            ax.scatter(x+off,frame.percent_of_bound,color=color,label=label,s=28)
            ax.scatter((x+off)[~done],frame.percent_of_bound.to_numpy()[~done],marker='x',color='black',s=70)
        if (g.percent_of_bound.dropna()>0).any() or (a.percent_of_bound.dropna()>0).any():ax.set_yscale('log')
        else:ax.text(.5,.5,'No observed counts available',transform=ax.transAxes,ha='center')
        ax.set_ylabel('Discovered states / named bound (%)');ax.legend(fontsize=8)
        for row in range(2):axes[row,col].set_xticks(x,[aliases[p] for p in order],rotation=90 if len(order)>12 else 0)
    savefig(fig,'graph_size_levels','How large is the explored state set relative to finite construction bounds? D counts discovered configurations, including the source. The input-only bound follows sphere–edge intersection counting; H_A H_B is the product of active single-curve vertices after anchored pruning. Missing counters remain missing. Crosses mark counts observed before termination. Neither completed search nor bound occupancy implies enumeration of all valid states; percentages do not measure progress. Formulas and assumptions are in graph_bounds.csv. Input-only bounds can be shown before execution; they are not observed graphs.')

def timing_figures(original,timing,summary,aliases,order):
    fig,axes=plt.subplots(1,2,figsize=(10.5,3.8),layout='constrained')
    for ax,f,title,denom in [(axes[0],original,'Original configurations',10),(axes[1],timing,'Fresh scheduled repetitions',30)]:
        bottom=np.zeros(4)
        groups=[('Validated',f.validated.eq(True),'#4daf4a'),('Resource exit',f.resource_limited.eq(True)&~f.validated.eq(True),'#e69f00'),('Pending',f.origin.eq('unexecuted'),'.85')]
        covered=groups[0][1]|groups[1][1]|groups[2][1]
        groups.append(('Other / error',~covered,'#e41a1c'))
        for label,mask,color in groups:
            vals=np.array([int((mask&f.method.eq(m)).sum()) for m in METHODS]);ax.bar(range(4),vals,bottom=bottom,label=label,color=color);bottom+=vals
        ax.set_xticks(range(4),[SHORT[m] for m in METHODS]);ax.set_ylim(0,denom*1.15);ax.set_title(title);ax.set_ylabel(f'Observations (denominator {denom} per method)')
        for i,m in enumerate(METHODS):ax.text(i,denom+.3,f'{int((f.validated.eq(True)&f.method.eq(m)).sum())}/{denom}',ha='center',fontsize=9)
    handles,labels=axes[1].get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=9,loc='outside lower center',ncol=4,frameon=False)
    savefig(fig,'completion','What completed under the fixed resource policy? Original configurations and fresh scheduled repetitions have separate denominators. Green means a saved output passed external feasibility validation; orange denotes resource termination, never mathematical infeasibility. Pending repetitions have no measured outcome. Exact termination reasons are retained in timing_observations.csv and the completion table.')
    fig,axes=plt.subplots(2,2,figsize=(11,6.5),sharex=True,layout='constrained')
    for ax,m in zip(axes.flat,METHODS):
        for i,pair in enumerate(order):
            g=timing[timing.pair.eq(pair)&timing.method.eq(m)]
            for _,r in g.iterrows():
                if r.origin=='unexecuted':continue
                x=i+(r.repetition-2)*.15
                ax.scatter(x,r.solver_seconds,marker='o' if r.validated else 'x',color=COLORS[m] if r.validated else 'black',s=30)
            v=g[g.validated.eq(True)].solver_seconds
            if len(v):ax.errorbar(i,v.median(),yerr=[[v.median()-v.min()],[v.max()-v.median()]],fmt='_',color='black',capsize=5,lw=1.2)
        ax.set_title(m);ax.set_yscale('log');ax.set_ylabel('Solver time (s)');ax.set_xticks(range(len(order)),[aliases[p] for p in order])
    axes[0,0].legend(handles=[Line2D([],[],marker='o',ls='',color='.3',label='Individual validated solution'),Line2D([],[],marker='x',ls='',color='black',label='Time to resource termination'),Line2D([],[],marker='_',color='black',label='Success median and range')],fontsize=7)
    savefig(fig,'timing_stability','How stable are repeated outcomes and solver times? Each point is one scheduled fresh repetition, offset by round. Bars show median and min–max only among validated successes, with counts and all outcomes in timing_summary.csv. Crosses are time to resource termination and are excluded from success statistics. Blank locations are unexecuted. Solver time retains internal checks; warm-up and external validation are separate. Repetitions are not independent input instances.')
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.1),layout='constrained',gridspec_kw={'width_ratios':[2,1]})
    for m in METHODS:
        for domain,marker in [('protein','o'),('hurricane','^')]:
            g=timing[timing.method.eq(m)&timing.domain.eq(domain)&timing.validated.eq(True)]
            axes[0].scatter(g.solver_seconds,g.retained_total,color=COLORS[m],marker=marker,s=30,alpha=.7)
        g=timing[timing.method.eq(m)&timing.resource_limited.eq(True)]
        axes[1].scatter(g.solver_seconds,[SHORT[m]]*len(g),marker='x',color=COLORS[m],s=50)
    axes[0].set_xscale('log');axes[0].set_xlabel('Successful solver time (s)');axes[0].set_ylabel(r'Retained vertices $(k_A+k_B)/(m+n)$');method_legend(axes[0])
    axes[0].set_title('Validated outputs: circles proteins; triangles hurricanes',fontsize=9)
    axes[1].set_xlabel('Time to resource termination (s)');axes[1].set_title('Resource exits: quality unknown',fontsize=9)
    if not timing.resource_limited.eq(True).any():axes[1].text(.5,.5,'No recorded resource exits yet',ha='center',transform=axes[1].transAxes,fontsize=8)
    savefig(fig,'compression_cost','What compression is obtained at the observed computational cost? Only validated fresh repetitions appear in the quality panel. Protein and hurricane markers differ; colors identify methods. Continuous fidelity and discrete fidelity define different feasible sets, so shorter CPS-2F outputs do not establish superiority for an identical optimization problem. Resource exits appear separately because no solution quality is available. No historical runtime is pooled into this timing plot.')
    fig,axes=plt.subplots(1,3,figsize=(11,3.7),layout='constrained')
    measured=timing[timing.origin.eq('new')].copy();measured['input_size']=measured.nA+measured.nB
    for ax,x,label in zip(axes,['input_size','states','transitions'],['Input vertices m+n','Discovered configurations','Candidate transition checks']):
        for m in METHODS:
            for valid,marker in [(True,'o'),(False,'x')]:
                g=measured[measured.method.eq(m)&measured.validated.eq(valid)]
                ax.scatter(g[x],g.solver_seconds,color=COLORS[m],marker=marker,s=25,alpha=.7)
        ax.set_xlabel(label);ax.set_yscale('log');ax.set_ylabel('Solver time (s)')
        if x!='input_size':ax.set_xscale('log')
    method_legend(axes[0])
    savefig(fig,'work_indicators','Which recorded work indicators accompany solver time? Fresh repetitions only; circles are validated successes and crosses are recorded unsuccessful terminations. Historical missing counters are not filled. Discovered states are not all valid states; transition checks are not unique graph edges or all machine operations. These are descriptive associations with no complexity-law fit. All protein inputs have the same sizes and share one structure.')
    fig,axes=plt.subplots(1,2,figsize=(10.5,3.6),layout='constrained')
    for i,m in enumerate(METHODS):
        g=measured[measured.method.eq(m)];axes[0].scatter([i]*len(g),g.peak_rss_mib,color=COLORS[m],alpha=.6,s=23)
        for j,key in enumerate(['warmup_seconds','solver_seconds','validation_seconds','total_seconds']):
            if g[key].notna().any():axes[1].scatter([j+(i-1.5)*.12]*len(g),g[key],color=COLORS[m],alpha=.5,s=16)
    axes[0].axhline(512,color='black',ls='--');axes[0].set_xticks(range(4),[SHORT[m] for m in METHODS]);axes[0].set_ylabel('Sampled peak process-tree RSS (MiB)')
    axes[1].set_xticks(range(4),['Warm-up','Solver','Validation','Total process'],rotation=20);axes[1].set_yscale('log');axes[1].set_ylabel('Phase / process time (s)')
    savefig(fig,'supplement_memory_phases','Supplement: sampled peak process-tree RSS and separate recorded times for fresh repetitions. Memory includes interpreter, compilation, allocator and monitoring-visible descendants, and is an approximate sampled peak, not exact algorithm memory. Phase values and unsuccessful outcomes remain individually available in the observation table; total process time is not solver time.')

def geometry(original):
    examples=[('AL031854__AL011991',METHODS[0]),('AL031854__AL011991',METHODS[2]),('EP111954__EP021984',METHODS[2]),('EP111954__EP021984',METHODS[3])]
    fig,axes=plt.subplots(2,2,figsize=(10.7,8),layout='constrained')
    for ax,(pair,m) in zip(axes.flat,examples):
        r=original[original.pair.eq(pair)&original.method.eq(m)].iloc[0];assert r.validated
        raw=read(ROOT/r.attempt_path/'raw_solver.json');result=raw['raw_result']
        inp=np.load(ROOT/r.input_file);A,B=inp['A'],inp['B']
        ia=result['indices_A'];ib=result['indices_B']
        for curve,indices,color,label in [(A,ia,'#2166ac','A'),(B,ib,'#b2182b','B')]:
            ax.plot(curve[:,0],curve[:,1],color=color,lw=.9,alpha=.4,label=f'{label} input')
            simp=curve[indices];ax.plot(simp[:,0],simp[:,1],color=color,marker='o',ms=3,lw=1.6,label=f"{label}' output")
        validation=read(ROOT/r.attempt_path/'validation.json');path=np.asarray(validation.get('coupling',[]),dtype=int)
        if len(path):
            left=A[ia][path[:,0]];right=B[ib][path[:,1]];worst=int(np.argmax(np.linalg.norm(left-right,axis=1)))
            line=np.array([left[worst],right[worst]])
            ax.plot(line[:,0],line[:,1],color='black',ls='--',lw=1,label='Bottleneck of saved coupling')
        ax.set_aspect('equal',adjustable='datalim');ax.set_xlabel('Basin-frame x (km)');ax.set_ylabel('Basin-frame y (km)')
        ax.set_title(f'{pair}\n{m}; k={int(r.k)}; coupling={r.coupling_distance:.2f} km\n'+r'$\delta_3$'+f'={r.delta3:.2f} km; ratio={r.coupling_distance/r.delta3:.3f}',fontsize=9)
    axes[0,0].legend(fontsize=7,ncol=2)
    savefig(fig,'geometric_examples','Illustrative selections, not representative averages. Top: Atlantic AL031854/AL011991 has equal k=5 for independent continuous and CPS-2F, but only the coupled output satisfies the coupling bound. Bottom: EP111954/EP021984 has CPS-2F k=16 and CPS-3F k=79, illustrating the change in feasible set between continuous and discrete fidelity. Pale curves are supplied inputs; strong curves and markers are selected original vertices. The dashed segment marks a bottleneck pair in the saved discrete coupling. Distances are rounded for display only; full-precision saved thresholds are unchanged. No independent storm alignment or temporal synchronization is inferred.')

def write_chapter(original,timing,summary,paired,inputs,bounds):
    schedule=read(TIMING/'schedule.json');finished=int(timing.origin.eq('new').sum());valid=int(timing.validated.eq(True).sum())
    resource=int(timing.resource_limited.eq(True).sum());pending=len(timing)-finished
    cont=paired[paired.fidelity.eq('continuous')];disc=paired[paired.fidelity.eq('discrete')]
    claims=[];macros={}
    def claim(name,value,source,selection,operation):
        macros[name]=str(value);claims.append(dict(claim=name,value=value,source=source,selection=selection,operation=operation))
    claim('PrimaryCount',len(original),'tables/original_cohort.csv','all','row count')
    claim('PrimaryValid',int(original.validated.eq(True).sum()),'tables/original_cohort.csv','validated=True','row count')
    claim('TimingPlanned',len(timing),'tables/timing_observations.csv','all','row count')
    claim('TimingFinished',finished,'tables/timing_observations.csv','origin=new','row count')
    claim('TimingValid',valid,'tables/timing_observations.csv','validated=True','row count')
    claim('TimingResource',resource,'tables/timing_observations.csv','resource_limited=True','row count')
    claim('TimingPending',pending,'tables/timing_observations.csv','origin=unexecuted','row count')
    claim('ContinuousMatched',int(cont.delta_k.notna().sum()),'tables/paired_quality.csv','fidelity=continuous; delta_k available','row count')
    claim('ContinuousEqual',int(cont.delta_k.eq(0).sum()),'tables/paired_quality.csv','fidelity=continuous; delta_k=0','row count')
    claim('ContinuousExtra',int(cont.delta_k.gt(0).sum()),'tables/paired_quality.csv','fidelity=continuous; delta_k>0','row count')
    claim('DiscreteEqual',int(disc.delta_k.eq(0).sum()),'tables/paired_quality.csv','fidelity=discrete; delta_k=0','row count')
    claim('EqualRepair',int(cont.equal_k_coupling_repair.sum()),'tables/paired_quality.csv','fidelity=continuous; equal_k_coupling_repair=True','row count')
    claim('BaselineViolations',int((cont.baseline_coupling_ratio>1+1e-9).sum()),'tables/paired_quality.csv','fidelity=continuous; coupling ratio>1+1e-9','row count')
    claim('ObjectiveDifferences',int(timing.objective_consistency.eq('unexpected_objective_difference').sum()),'tables/timing_observations.csv','objective_consistency=unexpected_objective_difference','row count')
    complete_groups=[g for _,g in timing.groupby('base_configuration_id') if g.origin.eq('new').all()]
    claim('CompleteTriples',len(complete_groups),'tables/timing_observations.csv','group by base_configuration_id; all three origin=new','group count')
    claim('StableTriples',sum(g.status.nunique()==1 for g in complete_groups),'tables/timing_observations.csv','all three recorded; distinct status count=1','group count')
    claim('MixedTriples',sum(g.status.nunique()>1 for g in complete_groups),'tables/timing_observations.csv','all three recorded; distinct status count>1','group count')
    active=bounds[bounds.quantity.eq('configuration_states_active_bound')]
    claim('GraphComparisons',len(active),'tables/graph_bounds.csv','quantity=configuration_states_active_bound','row count')
    claim('ActiveAvailable',int(active.finite_upper_bound.notna().sum()),'tables/graph_bounds.csv','quantity=configuration_states_active_bound; finite_upper_bound available','row count')
    a2=active[active.method.eq('CPS-2F')].percent_of_bound.dropna()
    claim('ActiveTwoMin',f'{a2.min():.2f}','tables/graph_bounds.csv','quantity=configuration_states_active_bound; method=CPS-2F','minimum percentage; rounded to 2 decimals')
    claim('ActiveTwoMax',f'{a2.max():.2f}','tables/graph_bounds.csv','quantity=configuration_states_active_bound; method=CPS-2F','maximum percentage; rounded to 2 decimals')
    for m,key in zip(METHODS,['IC','ID','TwoF','ThreeF']):
        claim('Valid'+key,int((original.validated.eq(True)&original.method.eq(m)).sum()),'tables/original_cohort.csv',f'method={m}; validated=True','row count')
    (OUT/'tables/numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros.items())+'\n',encoding='utf-8')
    pd.DataFrame(claims).to_csv(OUT/'tables/claim_traceability.csv',index=False)
    input_table=inputs[['pair_alias','pair','nA','nB','delta1','delta2','delta3']].copy()
    for key in ['delta1','delta2','delta3']:input_table[key]=input_table[key].map(lambda x:f'{x:.4f}')
    input_table.columns=['ID','Pair','m','n','delta1','delta2','delta3']
    latex_table(input_table,OUT/'tables/inputs.tex','Frozen inputs and thresholds. Protein units are angstroms; hurricane units are projected kilometres. Values here are rounded to four decimals for display; the CSV and frozen configurations retain full precision.','tab:inputs')
    counts=[]
    for cohort,frame in [('Original',original),('Timing',timing)]:
        for m in METHODS:
            g=frame[frame.method.eq(m)]
            counts.append(dict(Cohort=cohort,Method=SHORT[m],Planned=len(g),Validated=int(g.validated.eq(True).sum()),Resource=int(g.resource_limited.eq(True).sum()),Pending=int(g.origin.eq('unexecuted').sum())))
    counts=pd.DataFrame(counts);counts.to_csv(OUT/'tables/completion.csv',index=False)
    domain_counts=[]
    for cohort,frame in [('Original',original),('Timing',timing)]:
        for (domain,method),g in frame.groupby(['domain','method']):
            domain_counts.append(dict(cohort=cohort,domain=domain,method=method,planned=len(g),validated=int(g.validated.sum()),
                resource_limited=int(g.resource_limited.sum()),pending=int(g.origin.eq('unexecuted').sum())))
    pd.DataFrame(domain_counts).to_csv(OUT/'tables/completion_by_domain.csv',index=False)
    latex_table(counts,OUT/'tables/completion.tex','Validated coverage and resource exits, keeping the original cohort separate from fresh repetitions. Pending is not a failure.','tab:completion')
    tt=summary[['pair_alias','method','successful']].copy();tt['method']=tt.method.map(SHORT)
    for rep in [1,2,3]:tt[f'R{rep}']=summary[f'r{rep}_status'].replace({'validated':'valid','unexecuted':'pending','transition_budget_limit':'T-limit'})
    tt['Success median [min,max] s']=summary.apply(lambda r:'--' if not r.successful else f'{r.successful_solver_median:.4g} [{r.successful_solver_min:.4g},{r.successful_solver_max:.4g}]',axis=1)
    tt.columns=['Pair','Method','Valid n','R1','R2','R3','Success median [min,max] s']
    latex_table(tt,OUT/'tables/timing.tex','Fresh timing outcomes, with successful solver times only. Individual resource-exit times and all phase timings are available in timing_observations.csv. Three repetitions were scheduled per configuration; blank statistics do not mean zero.','tab:timing')
    quality=original[['pair_alias','method','kA','kB','k','retained_total','coupling_distance','delta3','status']].copy()
    quality['coupling_ratio']=quality.coupling_distance/quality.delta3
    quality.to_csv(OUT/'tables/quality_details.csv',index=False)
    qt=quality[['pair_alias','method','kA','kB','k','retained_total','coupling_ratio']].copy();qt['method']=qt.method.map(SHORT)
    for k in ['kA','kB','k']:qt[k]=qt[k].map(lambda x:num(x,4))
    for k in ['retained_total','coupling_ratio']:qt[k]=qt[k].map(lambda x:num(x,4))
    latex_table(qt,OUT/'tables/quality.tex','Original-cohort solution details. Independent coupling ratios above one are descriptive violations of a constraint the baseline does not enforce. Missing values correspond to no validated output.','tab:quality')
    # All per-repetition outcomes, including pending slots, have one record and explicit provenance.
    reasons=timing.groupby(['method','repetition','termination_reason'],dropna=False).size().rename('count').reset_index()
    reasons.to_csv(OUT/'tables/termination_reasons.csv',index=False)
    chapters=r'''\chapter{Experimental evaluation}
\label{chap:experiments}
\input{tables/numbers.tex}
\section{Experimental questions}
We evaluate the compression and computational cost of enforcing discrete output coupling under continuous or discrete fidelity. The primary questions are whether coupling changes the required number of vertices, how the two fidelity definitions change compression, and how many configuration states are discovered relative to finite construction bounds. Computational observations are bounded by a predefined resource policy, rather than a guarantee of exhaustive completion within that policy.

The completed primary cohort and the fresh timing study are distinct evidence sets. The primary cohort has \PrimaryValid{} externally validated solutions among \PrimaryCount{} configurations. At this draft's data snapshot, \TimingFinished{} of \TimingPlanned{} scheduled fresh repetitions are recorded: \TimingValid{} validated outputs, \TimingResource{} resource exits, and \TimingPending{} pending repetitions. Pending work is not counted as a solver failure. The timing results must not be described as a completed three-repetition study while any repetitions remain pending.

\section{Datasets and preprocessing}
We use ten previously frozen pairs, comprising five protein pairs and five hurricane pairs. Table~\ref{tab:inputs} identifies every pair and its numerical thresholds. The four original pilot pairs are P1, P2, H1 and H2; the other pairs were selected before their optimization results were observed, using input-size and geometric diversity among the supplied pairs, with individual reuse discouraged where possible. The recorded selection rule used standardized size, size-imbalance, chord-to-path and relative separation features with farthest-first selection. It could not avoid the common protein anchor.

Every protein curve has 22 vertices after the existing $w=16$ subsampling of the supplied prepared curves. All five pairs share the structure/chain \texttt{1o7j.a}; previously prepared alignment to that anchor is retained. These are neither full-resolution backbone experiments nor five independent biological replicates. They cannot establish protein input-size scaling. Hurricane inputs retain all supplied track vertices, without additional subsampling, cropping or independent alignment. Coordinates remain in the common projected basin frame. Spatial coordinates alone do not establish temporal synchronization. The manifest records preparation, input hashes, provenance and reuse.
\input{tables/inputs.tex}

\section{Compared methods and their constraints}
Let $A=(a_1,\ldots,a_m)$ and $B=(b_1,\ldots,b_n)$ denote the saved polygonal input curves. Outputs $A'$ and $B'$ are strictly ordered subsequences of original vertices, with both endpoints preserved. Define $k_A=|A'|$, $k_B=|B'|$, and the objective $k=\max(k_A,k_B)$. Unequal output lengths are allowed. A discrete coupling may repeat an index. Let $d_F$ denote continuous Fr\'echet distance and $d_{dF}$ discrete Fr\'echet distance, using the manuscript's notation~\cite{manuscript}.
\begin{itemize}
\item Independent continuous (IC) optimizes the two curves independently under $d_F(A,A')\leq\delta_1$ and $d_F(B,B')\leq\delta_2$.
\item Independent discrete (ID) uses $d_{dF}(A,A')\leq\delta_1$ and $d_{dF}(B,B')\leq\delta_2$ instead.
\item CPS-2F imposes the continuous fidelity constraints and $d_{dF}(A',B')\leq\delta_3$.
\item CPS-3F imposes the discrete fidelity constraints and the same discrete output-coupling constraint.
\end{itemize}
The independently optimal vertex counts minimize their maximum in the uncoupled problem. Their outputs need not satisfy $\delta_3$, and a coupling violation is not a baseline failure. Comparing CPS-2F with IC isolates the additional coupling requirement under continuous fidelity; comparing CPS-3F with ID gives the discrete counterpart. Comparing the two CPS methods changes the feasible set. Shorter CPS-2F outputs reflect its weaker continuous fidelity constraint, not superiority on an identical optimization problem. No local shortcut predicate replaces global continuous fidelity.

\section{Thresholds, resources, implementation and timing}
For each saved curve, $s_A$ and $s_B$ are its respective mean edge lengths. We use $\alpha=0.5$, $\delta_1=\alpha s_A$, $\delta_2=\alpha s_B$, and $\delta_3=d_{dF}(A,B)$ computed on the saved inputs without physical-unit rounding. Alpha was selected using earlier pilot evidence, not an independent held-out selection. The original curves are a feasible joint solution at these thresholds; a resource exit does not establish infeasibility.

All methods share a 60-second external solver allowance and a 512-MiB sampled process-tree RSS cap. CPS uses a 55-second internal soft limit, 100,000 configuration states, 20,000,000 candidate transition checks and 800,000 main single-curve edges. The edge cap applies during edge enumeration before anchored pruning; the reported completed-graph edge count is the post-pruning outgoing-entry count including waits. Independent discrete and the independent preprocessing within CPS-3F retain the existing 2,000,000-edge default per curve. Independent continuous has no analogous configuration-state or transition limit. Thus shared external limits coexist with method-specific internal work limits. The soft CPS time check is not an inner-transition-loop interrupt; the external supervisor is the backstop. These practical budgets are not theoretical sufficiency guarantees.

Separate budgets are retained for startup (60 s), input loading (15 s), warm-up (120 s), checkpointing (5 s), external validation (30 s), optional measurements (30 s), and reporting (10 s). The solver clock excludes explicit warm-up and external validation but includes any checks performed internally before return; no estimated overhead is subtracted. The continuous predicate's tolerance and numerical epsilon remain $10^{-10}$ and $10^{-9}$ respectively, with all other implementation tolerances unchanged and code hashes recorded.

The timing schedule fixes three rounds of every configuration and rotates method order across rounds while preserving pair order. Exactly one worker runs at a time. Each fresh process uses the same interpreter, thread settings (one OpenBLAS and OMP thread), isolated output-local cache policy and relevant-path warm-up. Existing process-aware readiness sampling precedes each four-method pair/round block and every resume; active workers are checked before every attempt. Readiness is outside measured phases. Pauses and timestamps remain in the progress history. No original observation is counted as a new repetition, and resource exits are not replaced by additional attempts.

Both certificate shortcuts are disabled with \texttt{certificate=False}. The implementation still computes its independent preprocessing inside the CPS solver, and that cost remains in solver time; it cannot bypass graph execution. No greedy preprocessing, alternative-subsequence search, vertex replacement or reference-guided pruning is used. Successful CPS records require explicit configuration-graph route evidence and raw status \texttt{optimal\_configuration\_graph}. Raw selected indices and counters are atomically saved before external validation. Validation checks original-vertex membership, bounds, strict order, endpoints, vertex counts, objective consistency, the applicable fidelity predicate and CPS output coupling. Optional achieved-distance searches occur afterwards. Feasibility validation alone does not independently prove optimality: optimality remains the algorithm's claim, supported separately by existing exhaustive small-instance correctness checks, not by the measured output's feasibility alone.

\section{Completion and computational performance}
Table~\ref{tab:completion} and Figure~\ref{fig:completion} retain explicit denominators. In the original cohort, IC, ID, CPS-2F and CPS-3F validated \ValidIC{}, \ValidID{}, \ValidTwoF{} and \ValidThreeF{} configurations respectively, each out of ten. The one original missing CPS-2F output is H4 (\texttt{AL072006\_\_AL011900}), which reached the transition budget after completing its single-curve graphs. Its reported work counts are observations before termination, not a completed product-graph size. Resource-limited solver times are time to termination, not time to a successful solution.
\input{tables/completion.tex}
@@FIG:completion@@
Figure~\ref{fig:timing_stability} displays individual fresh repetitions, with medians and ranges computed only from validated successes. Table~\ref{tab:timing} exposes all scheduled outcomes and the successful count for each statistic. Among \CompleteTriples{} configurations with all repetitions recorded, \StableTriples{} have the same outcome in every round and \MixedTriples{} have mixed outcomes. Across the recorded repetitions, \ObjectiveDifferences{} unexpected objective differences were detected relative to the saved primary results or other validated repetitions. This statement applies only to recorded outputs; missing repetitions cannot support a stability conclusion. Different selected subsequences can be legitimate alternative optima.
@@FIG:timing_stability@@
Figure~\ref{fig:compression_cost} separates the quality--time observations from resource exits with unknown output quality. Figure~\ref{fig:work_indicators} relates time to recorded input size, discovered states and transition checks without fitting a scaling law. The input families and computational quantities are too restricted to infer general asymptotic behavior. Sampled memory and phase times appear in the supplement; sampled process-tree RSS includes runtime and compilation overhead and can miss unsampled peaks.
@@FIG:compression_cost@@
@@FIG:work_indicators@@

\section{Simplification quality and the effect of coupling}
The original cohort supplies \ContinuousMatched{} matched validated continuous-fidelity comparisons. CPS-2F and IC have equal $k$ on \ContinuousEqual{} of these, while \ContinuousExtra{} requires an additional vertex in the maximum output length. All \DiscreteEqual{} discrete-fidelity comparisons have equal CPS-3F and ID objectives. This null vertex-cost result is retained, rather than selecting only cases with a difference. The continuous independent outputs violate coupling in \BaselineViolations{} cases, and \EqualRepair{} matched cases restore coupling with unchanged $k$. Figure~\ref{fig:coupling_vertex_cost} shows both the normalized coupling and the matched-fidelity vertex cost. These counts describe this selected cohort, not a general guarantee that coupling is free.
@@FIG:coupling_vertex_cost@@
Figure~\ref{fig:geometric_examples} uses explicitly illustrative selections. The Atlantic example demonstrates that equal $k$ does not imply equal coupling feasibility. The second hurricane example illustrates the compression difference induced by continuous versus discrete fidelity. Neither is a representative average. The reported retained-vertex fraction is $(k_A+k_B)/(m+n)$; separate output sizes remain in Table~\ref{tab:quality}.
@@FIG:geometric_examples@@

\section{Observed graph size relative to finite bounds}
For a nondegenerate curve with $m$ original vertices, each sphere around an original vertex intersects each straight edge in at most two locations. The construction assumes no edge is contained in a sphere boundary; the implementation skips degenerate edges and retains original endpoints. Thus the number of added matching locations satisfies $a_A\leq 2m(m-1)$, and similarly $a_B\leq 2n(n-1)$~\cite[Observation 4, printed p.~18 (PDF p.~26)]{thesis}. Counts refer to distinct interior curve parameters after tolerance-based merging, not distinct spatial coordinates. Auxiliary matching locations are not selectable output vertices.

Single-curve dense index spaces have sizes $S_A=m(m+a_A)$ and $S_B=n(n+a_B)$. Anchored reachability pruning leaves $H_A$ and $H_B$ active vertices, including source and terminal. A product configuration chooses one vertex from each graph, giving the finite bounds
\begin{align}
B_{\mathrm{input}}^{2F}&=mn\,[m+2m(m-1)]\,[n+2n(n-1)],\\
B_{\mathrm{extended}}^{2F}&=mn(m+a_A)(n+a_B),\\
B_{\mathrm{active}}&=H_AH_B.
\end{align}
For CPS-3F there are no auxiliary locations and $B_{\mathrm{input}}^{3F}=m^2n^2$. Source and terminal require no additional states outside these products. These are combinatorial counts from the construction, not evaluations of an $O(\cdot)$ expression with an assumed constant of one. The theoretical CPS-2F $O(m^3n^3)$ expression counts potential configurations; see the thesis section 5.1, printed p.~26 (PDF p.~34), and manuscript section 4.4, PDF pp.~15--16~\cite{thesis,manuscript}.

Let $D$ be the number of discovered configurations, including the source. Figure~\ref{fig:graph_size_levels} compares $D$, $H_AH_B$, and $B_{\mathrm{input}}$ and names the denominator for each percentage. Active products can still include coupling-invalid states. A completed search need not enumerate all possible, valid or reachable configurations. For interrupted runs, $D$ is only observed before termination. Neither $100D/B_{\mathrm{input}}$ nor $100D/(H_AH_B)$ is an execution-progress estimate. Missing historical active-state counters remain missing and are not backfilled with measurements from timing repetitions.

Active-product bounds are available for \ActiveAvailable{} of the \GraphComparisons{} original CPS observations. Among the recorded CPS-2F active-product comparisons, discovered-state occupancy ranges from \ActiveTwoMin\% to \ActiveTwoMax\%, including the interrupted observation. These percentages describe the explored subset relative to a particular combinatorial denominator, not a universal sparsity ratio or a proof that all valid states were searched.

Single-graph edge counts include one wait self-transition per active vertex. Their finite bound is $E_A\leq H_A(H_A+1)/2$ (and analogously for B), because each non-wait edge strictly advances the topological index and entries are unique. Where both completed graphs have counters, candidate transition checks are bounded by $E_AE_B-H_AH_B$: each outgoing-entry pair defines one product transition candidate, simultaneous waits are excluded and a source is popped at most once. The implemented counter increments after the coupling filter; a budget-triggering check may be counted before relaxation. This does not establish that checks are stored unique graph edges. The product graph is implicit and no observed unique-edge count is available. Dynamic-programming labels are another quantity: at most $m$ A-hop labels per configuration in this implementation, giving a finite $mH_AH_B$ bound without an observed label count. Full formulas, assumptions and missing-data reasons are in the graph-bound table.
@@FIG:graph_size_levels@@

\section{Limitations}
The cohort is small, selected for a bounded pilot, and uses a pilot-informed alpha. Shared protein structure and identical protein sizes prevent claims of independent biological replication or protein scaling. Storm reuse and coordinate preparation are recorded; spatial coupling is not evidence of temporal correspondence. Three repeated executions of an input do not create three independent input instances or support a claim of statistical significance. Outcome variation must remain visible, and failed-run times must not enter successful-runtime averages.

Continuous and discrete fidelity are different scientific constraints. Feasibility checks do not independently certify optimality; resource exits do not certify infeasibility. External budgets are shared but internal work limits are method-specific. Historical instrumentation differs, historical missing counters remain missing, and the timing study is kept separate from original single observations. Process-aware readiness cannot guarantee an otherwise idle host throughout each measured phase. No general complexity law, exact packedness value, topological-preservation conclusion or publication-scale performance guarantee follows from these results. The provisional experimental prose in the supplied documents is not treated as measured evidence.

\section*{Supplementary saved-result tables and figure}
\input{tables/timing.tex}
\input{tables/quality.tex}
@@FIG:supplement_memory_phases@@
\begin{thebibliography}{9}
\bibitem{thesis} G. Gozoltzani. Supplied MSc thesis, \emph{M\_Sc\_\_Thesis\_\_\_Galit\_Gozoltzani.pdf}. Construction references: Observation 4 and section 5.1. The file hash is recorded in the accompanying provenance.
\bibitem{manuscript} \emph{Chain Pair Simplification under the Continuous Fr\'echet Distance}. Supplied manuscript, sections 3.1, 4.1 and 4.4. Experimental text is provisional; the file hash is recorded in the accompanying provenance.
\end{thebibliography}
'''
    for name,caption in CAPTIONS.items():
        figure='\n'.join([r'\begin{figure}[p]',r'\centering',r'\includegraphics[width=\textwidth]{figures/'+name+'.pdf}',r'\caption{'+tex(caption)+'}',r'\label{fig:'+name+'}',r'\end{figure}'])
        chapters=chapters.replace('@@FIG:'+name+'@@',figure)
    assert '@@FIG' not in chapters
    (OUT/'experiments.tex').write_text(chapters,encoding='utf-8')
    (OUT/'main.tex').write_text(r'''\documentclass[11pt]{report}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[a4paper,margin=22mm]{geometry}
\usepackage{amsmath,amssymb,booktabs,longtable,graphicx,hyperref}
\title{CPS experimental evaluation: evidence-backed draft}
\author{}
\date{Saved-data snapshot}
\begin{document}
\maketitle
\input{experiments.tex}
\end{document}
''',encoding='utf-8')
    (OUT/'figure_captions.md').write_text('\n\n'.join(f'## {k}\n\n{v}' for k,v in CAPTIONS.items()),encoding='utf-8')
    report=f'''# Experimental results for thesis and manuscript

## Status
Original cohort: {int(original.validated.sum())}/{len(original)} validated configurations. Fresh timing study: **{finished}/{len(timing)} recorded**, {valid} validated, {resource} resource-limited, {pending} pending. Three rounds are scheduled; inherited observations are not repetitions. This draft is provisional while timing slots remain pending.

## Supported quality findings (original cohort)
- Continuous-fidelity comparison: {int(cont.delta_k.eq(0).sum())}/{int(cont.delta_k.notna().sum())} matched cases have equal k; {int(cont.delta_k.gt(0).sum())} has a larger joint k.
- Discrete-fidelity comparison: all {int(disc.delta_k.eq(0).sum())} matched cases have equal k.
- {int(cont.equal_k_coupling_repair.sum())} continuous-fidelity cases repair coupling without increasing k. Independent coupling violations are permitted by their problem definition.
- The Atlantic illustrative example has equal k=5 with different coupling feasibility. EP111954/EP021984 illustrates CPS-2F k=16 versus CPS-3F k=79 under different fidelity constraints.
- {int(timing.objective_consistency.eq('unexpected_objective_difference').sum())} unexpected objective differences among recorded validated timing outputs. Pending slots support no timing conclusion.

## Interpretation and limits
The state figure distinguishes input-only finite bounds, active single-graph products and discovered states. Missing historical counters remain missing. Occupancy is not search completion. Transition checks are not unique edges. Alpha was chosen using pilot evidence; all protein inputs have 22 vertices and share 1o7j. These data do not establish general scaling, independent biological replication, packedness, topological preservation or statistical significance.

## Files and traceability
- `main.tex` includes the new `experiments.tex`; the user's original file is untouched.
- `tables/timing_observations.csv`: all 120 planned repetitions, individual phases, raw-result paths and outcomes.
- `tables/timing_summary.csv`: 40 configurations, all three statuses, success-only median/range/count, resource-exit times separately.
- `tables/claim_traceability.csv`: generated numerical chapter claims and source filters; table entries derive directly from CSVs.
- `tables/original_cohort.csv` and `tables/graph_bounds.csv`: preserved primary evidence, copied for this report.
- `figures/`: PDF and PNG; `figure_captions.md`: self-contained captions.
- `provenance/`: schedule/source hashes and build audit.

## Execution policy
The frozen schedule interleaves methods across three rounds. One worker at a time, unchanged warm-up/cache/thread policy, certificate=False, no algorithm changes, no replacement repetitions. Shared external time/memory caps have disclosed method-specific work limits. The schedule is resumable and completed slots are never duplicated.
'''
    (OUT/'findings.md').write_text(report,encoding='utf-8')
    progress=read(TIMING/'progress.json')
    with (OUT/'findings.md').open('a',encoding='utf-8') as f:
        f.write('\n## Recorded batch state\n'+f"Phase: {progress['phase']}. Last progress update: {progress['utc']}. "+str(progress.get('pause_reason',''))+'\n')
        if pending:
            f.write('Pending configurations (no replacement attempts or additional rounds):\n\n')
            for _,r in timing[timing.origin.eq('unexecuted')].iterrows():f.write(f'- {r.pair}, {r.method}, repetition {int(r.repetition)}.\n')
        f.write('\nThe readiness records preserve temporary PyCharm contention and all pauses. No unrelated processes were terminated and no system settings were changed.\n')
    (OUT/'status.md').write_text(f'# Publication preparation status\n\nOriginal cohort: 39/40 validated. Timing: {finished}/120 recorded, {valid} validated, {resource} resource exits, {pending} pending.\n\n'+
        ('Blocked by sustained competing host activity; the saved schedule can resume without duplicating completed repetitions.\n' if pending else 'All three scheduled rounds complete; no additional repetitions authorized.\n')+
        '\nResults, figures, generated tables, LaTeX draft and saved-results notebook are prepared from the current snapshot. No optimization algorithm was changed.\n',encoding='utf-8')
    totals={key:float(timing[key].sum()) for key in ['warmup_seconds','solver_seconds','validation_seconds','total_seconds']}
    dump(OUT/'tables/phase_resource_summary.json',dict(recorded=finished,phase_sums_seconds=totals,
        phase_sums_note='Includes recorded resource termination times; not a successful-runtime statistic.',
        sampled_peak_tree_rss_max_mib=float(timing.peak_rss_mib.max()),memory_includes_runtime_overhead=True,
        progress_snapshot=progress))
    sources=read(SOURCE/'source_references.json')
    original_tex=Path(sources[0]['path']).parent/'experiments.tex'
    if original_tex.exists():sources.append(dict(path=str(original_tex),sha256=sha(original_tex),purpose='Notation and context only; provisional numerical claims not inherited'))
    dump(OUT/'provenance/research_sources.json',sources)
    source_paths=[SOURCE/'results.csv',SOURCE/'graph_bounds.csv',TIMING/'results.csv',TIMING/'schedule.json',TIMING/'progress.json',Path(__file__)]
    dump(OUT/'provenance/build.json',dict(utc=datetime.now(timezone.utc).isoformat(),sources={str(p.relative_to(ROOT)):sha(p) for p in source_paths},schedule_environment=schedule['environment'],execution_code=schedule['execution_code'],optimization_executed_by_reporting=False))
    shutil.copyfile(SOURCE/'graph_bounds.md',OUT/'tables/graph_bound_derivation.md')
    evidence_audit(original,timing,schedule);structural_check();update_notebook();render_qa()

def evidence_audit(original,timing,schedule):
    records=[]
    assert len(timing)==120 and timing.configuration_id.is_unique
    assert timing.groupby('base_configuration_id').size().eq(3).all()
    assert sha(SOURCE/'results.csv')==schedule['original_results_sha256']
    for name,h in schedule['execution_code'].items():assert sha(ROOT/'notebooks_paper/graph_only_v3'/name)==h
    for _,r in timing.iterrows():
        if r.origin=='unexecuted':continue
        d=ROOT/r.attempt_path;config=read(d/'configuration.json');planned=next(x['configuration'] for x in schedule['rows'] if x['configuration_id']==r.configuration_id)
        assert {k:v for k,v in config.items() if k not in ['attempt_id','code_version']}==planned
        assert config['certificate'] is False and config['alpha']==.5
        assert sha(ROOT/config['input_file'])==config['input_sha256']
        for p,h in config['algorithms'].items():assert sha(ROOT/p)==h
        assert len(list(d.parent.glob('*/configuration.json')))==1
        raw=d/'raw_solver.json';val=d/'validation.json';before=None
        if r.validated:
            assert raw.exists() and val.exists();a=read(raw);b=read(val)
            assert a['certificate_enabled'] is False and b['feasibility_passed'] and b['status']=='validated'
            assert b['raw_checkpoint_sha256']==sha(raw)
            assert b['k']==r.k and b['kA']==r.kA and b['kB']==r.kB
            before=a['unix']<=b['unix'];assert before
            if r.method.startswith('CPS'):
                assert a['solver_status']=='optimal_configuration_graph'
                assert a['observation']['configuration_graph_entered'] is True and b['coupling_feasible']
        else:assert pd.isna(r.k) and pd.isna(r.kA) and pd.isna(r.kB)
        records.append(dict(configuration_id=r.configuration_id,attempt_path=r.attempt_path,configuration_hash=sha(d/'configuration.json'),
            input_hash_verified=True,code_hashes_verified=True,one_attempt=True,raw_saved_before_validation=before,
            raw_hash=sha(raw) if raw.exists() else None,validation_hash=sha(val) if val.exists() else None,status=r.status))
    assert not timing.objective_consistency.eq('unexpected_objective_difference').any()
    pd.DataFrame(records).to_csv(OUT/'provenance/attempt_audit.csv',index=False)
    dump(OUT/'provenance/audit_summary.json',dict(recorded=len(records),planned=120,original_results_hash_unchanged=True,
        three_scheduled_per_configuration=True,no_duplicate_attempts=True,algorithm_hashes_unchanged=True,
        validated=int(timing.validated.sum()),unexpected_objective_differences=0,
        feasibility_validation_not_independent_optimality_proof=True))

def structural_check():
    files=list(OUT.glob('*.tex'))+list((OUT/'tables').glob('*.tex'))
    # Follow generated table inclusions outside the chapter directory as well.
    seen={p.resolve() for p in files};i=0
    while i<len(files):
        for target in re.findall(r'\\input\{([^}]+)\}',files[i].read_text(encoding='utf-8')):
            p=OUT/target
            if not p.suffix:p=p.with_suffix('.tex')
            if p.exists() and p.resolve() not in seen:files.append(p);seen.add(p.resolve())
        i+=1
    alltext='\n'.join(p.read_text(encoding='utf-8') for p in files)
    missing=[]
    for command,target in re.findall(r'\\(input|includegraphics)(?:\[[^\]]*\])?\{([^}]+)\}',alltext):
        p=OUT/target
        if command=='input' and not p.suffix:p=p.with_suffix('.tex')
        if not p.exists():missing.append(str(p))
    labels=set(re.findall(r'\\label\{([^}]+)\}',alltext));refs=set(re.findall(r'\\ref\{([^}]+)\}',alltext))
    env_errors=[];brace_errors=[]
    for p in files:
        stack=[]
        body=p.read_text(encoding='utf-8')
        for command,env in re.findall(r'\\(begin|end)\{([^}]+)\}',body):
            if command=='begin':stack.append(env)
            elif not stack or stack.pop()!=env:env_errors.append(str(p))
        if stack:env_errors.append(str(p))
        level=0
        for token in re.findall(r'(?<!\\)[{}]',body):
            level+=1 if token=='{' else -1
            if level<0:brace_errors.append(str(p));break
        if level:brace_errors.append(str(p))
    tools={name:shutil.which(name) for name in ['pdflatex','xelatex','lualatex','tectonic']}
    citations=set(re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}',alltext))
    citations={c for group in citations for c in group.split(',')}
    bib=set(re.findall(r'\\bibitem\{([^}]+)\}',alltext))
    result=dict(missing_files=missing,unresolved_references=sorted(refs-labels),unresolved_citations=sorted(citations-bib),environment_errors=env_errors,brace_errors=brace_errors,latex_tools=tools,compiled=False)
    assert not missing and not refs-labels and not env_errors and not brace_errors and not citations-bib,result
    dump(OUT/'qa/latex_structural_check.json',result)
    with (OUT/'findings.md').open('a',encoding='utf-8') as f:f.write('\n## LaTeX checks\nFile inclusions, references and environment structure passed. No LaTeX compiler was found on PATH; the draft has not been compiled. PDF figures are rendered separately for visual inspection.\n')

def update_notebook():
    import nbformat
    p=ROOT/'notebooks_paper/expanded_comparison_saved_results.ipynb'
    backup=OUT/'provenance/expanded_comparison_notebook_before_publication.ipynb'
    if not backup.exists():shutil.copyfile(p,backup)
    nb=nbformat.read(p,as_version=4)
    # Reuse this notebook as the concise presentation; the prior reader is archived above.
    cells=[]
    def add(kind,text):
        cell=nbformat.v4.new_markdown_cell(text) if kind=='markdown' else nbformat.v4.new_code_cell(text)
        cell.metadata['publication_saved_results']=True;cells.append(cell)
    add('markdown','# Thesis/manuscript presentation — saved evidence only\n\nResearch questions: coupling cost, fidelity/compression trade-off, explored graph size and timing stability. The original cohort and fresh timing repetitions remain separate. No cell below launches optimization. Full results and limitations: `output/paper_v1/experimental_chapter/findings.md`.')
    add('code',"from pathlib import Path\nimport pandas as pd, json\nfrom IPython.display import display, Image, Markdown\nPROJECT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p/'output/paper_v1/experimental_chapter').exists())\nPUBLICATION = PROJECT/'output/paper_v1/experimental_chapter'  # central display parameter\ndisplay(Markdown((PUBLICATION/'findings.md').read_text(encoding='utf-8')))\ndisplay(pd.read_csv(PUBLICATION/'tables/completion.csv'))")
    add('code',"display(pd.read_csv(PUBLICATION/'tables/timing_summary.csv')[['pair_alias','method','successful','r1_status','r2_status','r3_status','successful_solver_median','successful_solver_min','successful_solver_max']])")
    for name,caption in CAPTIONS.items():
        add('markdown',f'## {name.replace("_"," ").capitalize()}\n\n{caption}')
        add('code',f"display(Image(filename=str(PUBLICATION/'figures/{name}.png')))")
    add('markdown','Full table links: [timing](../output/paper_v1/experimental_chapter/tables/timing_observations.csv), [quality](../output/paper_v1/experimental_chapter/tables/quality_details.csv), [graph bounds](../output/paper_v1/experimental_chapter/tables/graph_bounds.csv). Feasibility validation and algorithm optimality claims remain distinct.')
    nb.cells=cells;nbformat.write(nb,p)

def render_qa():
    import fitz
    from PIL import Image,ImageOps,ImageDraw
    thumbs=[]
    for p in sorted((OUT/'figures').glob('*.pdf')):
        doc=fitz.open(p);pix=doc[0].get_pixmap(matrix=fitz.Matrix(1,1),alpha=False)
        path=OUT/'qa'/f'{p.stem}.png';pix.save(str(path));doc.close()
        im=Image.open(path).convert('RGB');im.thumbnail((600,440));tile=Image.new('RGB',(620,470),'white');tile.paste(im,((620-im.width)//2,25));ImageDraw.Draw(tile).text((10,5),p.stem,fill='black');thumbs.append(tile)
    sheet=Image.new('RGB',(1240,470*((len(thumbs)+1)//2)),'white')
    for i,im in enumerate(thumbs):sheet.paste(im,((i%2)*620,(i//2)*470))
    sheet.save(OUT/'qa/contact_sheet.png')
    print('Saved publication draft and figures. Reporting did not run optimization.',flush=True)

def main():
    for d in ['tables','figures','provenance','qa']: (OUT/d).mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.facecolor':'white','figure.facecolor':'white','text.color':'black','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    original,timing,summary,paired,inputs,bounds,aliases,order=summaries()
    coupling_figure(paired,aliases,order);bound_figure(bounds,aliases,order)
    timing_figures(original,timing,summary,aliases,order);geometry(original)
    write_chapter(original,timing,summary,paired,inputs,bounds)

def verify_notebook():
    import nbformat
    from nbclient import NotebookClient
    p=ROOT/'notebooks_paper/expanded_comparison_saved_results.ipynb'
    nb=nbformat.read(p,as_version=4)
    for c in nb.cells:
        if c.cell_type=='code':assert not any(word in c.source for word in ['subprocess','solve_cps','continuous_decision','timing_execute'])
    NotebookClient(nb,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    nbformat.write(nb,p)
    dump(OUT/'qa/notebook_execution.json',dict(passed=True,code_cells=sum(c.cell_type=='code' for c in nb.cells),optimization_executed=False,notebook_sha256=sha(p)))
    build=read(OUT/'provenance/build.json');build['sources'][str(Path(__file__).relative_to(ROOT))]=sha(__file__);dump(OUT/'provenance/build.json',build)
    dump(OUT/'provenance/artifact_hashes.json',{str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_hashes.json'})
    print('Saved-results notebook executed successfully; no optimization calls.',flush=True)

def coverage_report():
    """Combine original single observations and extension; keep timing separate."""
    import sys,nbformat
    sys.path[:0]=[str(ROOT/'notebooks_paper'),str(ROOT/'notebooks_paper/graph_only_v3')]
    import coverage_extension as cx
    import expanded_comparison as shared
    global OUT
    chapter=OUT;ext=cx.OUT
    (ext/'qa').mkdir(exist_ok=True)
    old=pd.read_csv(SOURCE/'results.csv');fresh=cx.save()
    added=ROOT/'output/paper_v1/hurricane_extension_3';added_inputs=[]
    if (added/'configuration_manifest.json').exists():
        cx.OUT=added
        try:more=cx.save()
        finally:cx.OUT=ext
        more['observation_cohort']='three-hurricane extension single observation'
        fresh=pd.concat([fresh,more],ignore_index=True,sort=False)
        added_inputs=read(added/'input_manifest.json')['inputs']
    final=ROOT/'output/paper_v1/hurricane_final_20';final_inputs=[]
    if (final/'configuration_manifest.json').exists():
        cx.OUT=final
        try:final_rows=cx.save()
        finally:cx.OUT=ext
        final_rows['observation_cohort']='stratified final20 single observation'
        fresh=pd.concat([fresh,final_rows],ignore_index=True,sort=False)
        final_inputs=read(final/'input_manifest.json')['inputs']
    old['source_origin']=old.origin;old['origin']='reused';old['observation_cohort']='original alpha0.5 single observation'
    old['eligible']=True;fresh['eligible']=~fresh.origin.eq('blocked')
    f=pd.concat([old,fresh],ignore_index=True,sort=False);f['validated']=f.validated.fillna(False).astype(bool);f['resource_limited']=f.resource_limited.fillna(False).astype(bool)
    assert len(f)==96+4*len(added_inputs)+4*len(final_inputs) and f.configuration_id.is_unique and f.alpha.eq(.5).all()
    f['attempt_path']=f.attempt_path.where(f.attempt_path.notna(),None)
    f['observed']=f.origin.isin(['new','reused'])
    oldinputs=read(SOURCE/'input_manifest.json')['inputs'];newinputs=read(ext/'input_manifest.json')['inputs']
    inputs=[dict(r,eligible=True,coverage_cohort='original') for r in oldinputs]+[dict(r,coverage_cohort='extension') for r in newinputs]+[dict(r,coverage_cohort='three-hurricane extension') for r in added_inputs]+[dict(r,coverage_cohort='stratified final20') for r in final_inputs]
    aliases={};ordered=[]
    for domain,prefix in [('protein','P'),('hurricane','H')]:
        for i,r in enumerate([x for x in inputs if x['domain']==domain],1):aliases[r['pair']]=f'{prefix}{i}';ordered.append(r['pair'])
    f['pair_alias']=f.pair.map(aliases)
    gallery_table=ext/'attempt_gallery/per_attempt_metrics.csv'
    if gallery_table.exists():
        gm=pd.read_csv(gallery_table);lookup={str(r.attempt_path).replace('\\','/'):r for r in gm.itertuples()}
        for field in ['figure_png','figure_pdf']:
            f[field]=[getattr(lookup.get(str(p).replace('\\','/')),field,None) for p in f.attempt_path]
        f['attempt_metrics_json']=[str((ext/'attempt_gallery/metrics'/f'{lookup[str(p).replace(chr(92),chr(47))].attempt_id}.json').relative_to(ROOT)) if str(p).replace('\\','/') in lookup else None for p in f.attempt_path]
    f.to_csv(ext/'combined_results.csv',index=False)
    input_frame=pd.DataFrame(inputs);input_frame['pair_alias']=input_frame.pair.map(aliases);input_frame.to_csv(ext/'combined_input_manifest.csv',index=False)
    reuse=[]
    for domain in ['protein','hurricane']:
        members={}
        for r in inputs:
            if r['domain']!=domain:continue
            for identity in r['individuals']:members.setdefault(identity,[]).append(r['pair'])
        for identity,pairs in members.items():reuse.append(dict(domain=domain,identity=identity,identity_type='PDB structure' if domain=='protein' else 'storm',pair_count=len(pairs),pairs=json.dumps(pairs)))
    pd.DataFrame(reuse).to_csv(ext/'identity_reuse.csv',index=False)
    saved_out=shared.OUT;shared.OUT=ext
    bounds=shared.graph_tables(f[f.eligible].copy());shared.OUT=saved_out
    bounds['pair_alias']=bounds.pair.map(aliases);bounds.to_csv(ext/'graph_bounds.csv',index=False)
    comparisons=[]
    for pair in ordered:
        sub=f[f.pair.eq(pair)].set_index('method')
        for fidelity,base,joint in [('continuous',METHODS[0],METHODS[2]),('discrete',METHODS[1],METHODS[3])]:
            a,b=sub.loc[base],sub.loc[joint];both=bool(a.validated and b.validated)
            comparisons.append(dict(pair=pair,pair_alias=aliases[pair],domain=a.domain,fidelity=fidelity,cohort=a.observation_cohort,
                baseline_k=a.k if a.validated else None,joint_k=b.k if b.validated else None,delta_k=b.k-a.k if both else None,
                baseline_coupling_ratio=a.coupling_distance/a.delta3 if a.validated else None,
                joint_coupling_ratio=b.coupling_distance/b.delta3 if b.validated else None,
                equal_k_coupling_repair=bool(both and not a.coupling_feasible and b.coupling_feasible and a.k==b.k),
                baseline_status=a.status,joint_status=b.status,baseline_source=a.attempt_path,joint_source=b.attempt_path))
    paired=pd.DataFrame(comparisons);assert (paired.delta_k.dropna()>=0).all();paired.to_csv(ext/'matched_comparisons.csv',index=False)
    completion=f.groupby(['domain','method','termination_reason'],dropna=False).size().rename('count').reset_index()
    completion['planned_denominator']=completion.domain.map(f.groupby('domain').pair.nunique())
    completion.to_csv(ext/'completion.csv',index=False)
    OUT=ext;CAPTIONS.clear()
    plt.rcParams.update({'font.size':10,'axes.facecolor':'white','figure.facecolor':'white','text.color':'black','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    statuses=['validated']+sorted(set(f.termination_reason.dropna())-{'validated'})
    statuscolors={'validated':'#4daf4a','transition_budget_limit':'#e69f00','state_budget_limit':'#984ea3','single_chain_edge_budget_limit':'#377eb8','preparation_blocked':'.65','unexecuted':'.92','validation_incomplete':'#cc79a7','solver_timeout':'#d55e00','memory_limit':'#a65628'}
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    for ax,domain in zip(axes,['protein','hurricane']):
        sub=f[f.domain.eq(domain)];bottom=np.zeros(4);denom=len(sub)//4
        for status in statuses:
            vals=np.array([int((sub.method.eq(m)&sub.termination_reason.eq(status)).sum()) for m in METHODS])
            ax.bar(range(4),vals,bottom=bottom,label=status.replace('_',' '),color=statuscolors.get(status,'#e41a1c'),hatch='//' if status=='preparation_blocked' else None);bottom+=vals
        for i,m in enumerate(METHODS):ax.text(i,denom+.2,f'{int((sub.method.eq(m)&sub.validated).sum())}/{denom}',ha='center',fontsize=9)
        ax.set_ylim(0,denom*1.15);ax.set_xticks(range(4),[SHORT[m] for m in METHODS]);ax.set_title(f'{domain.capitalize()}: {denom} planned pairs');ax.set_ylabel('Configuration slots')
    handles,labels=axes[1].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=8,frameon=False)
    savefig(fig,'coverage_completion','Combined single-observation coverage with denominators of 14 planned protein pairs and 10 hurricane pairs per method. Hatched slots are preparation-blocked, not experiments or infeasible instances. Resource limits retain their exact reasons. The original 40 observations retain their provenance; the 120 timing repetitions are not pooled into this figure.')
    fig,axes=plt.subplots(2,1,figsize=(11,6),layout='constrained')
    for ax,domain in zip(axes,['protein','hurricane']):
        pairs=[p for p in ordered if f[f.pair.eq(p)].iloc[0].domain==domain]
        for j,m in enumerate(METHODS):
            s=f[f.domain.eq(domain)&f.method.eq(m)].set_index('pair').loc[pairs];x=np.arange(len(pairs))+(j-1.5)*.17
            good=s.validated.to_numpy();ax.scatter(x[good],s.k.to_numpy()[good],s=30,color=COLORS[m],label=m)
            bad=good&~s.coupling_feasible.fillna(False).to_numpy(dtype=bool)
            if m.startswith('Independent'):ax.scatter(x[bad],s.k.to_numpy()[bad],marker='x',color='black',s=75)
        for i,p in enumerate(pairs):
            s=f[f.pair.eq(p)]
            if not s.eligible.any():ax.axvspan(i-.45,i+.45,color='.9',hatch='//',zorder=-2)
        ax.set_xticks(range(len(pairs)),[aliases[p] for p in pairs],rotation=90);ax.set_title(domain.capitalize());ax.set_ylabel(r'$k=\max(k_A,k_B)$')
    h,l=axes[0].get_legend_handles_labels();h.append(Line2D([],[],marker='x',ls='',color='black'));l.append('Independent coupling exceeds bound')
    fig.legend(h,l,loc='outside lower center',ncol=3,fontsize=8,frameon=False)
    savefig(fig,'coverage_quality','Matched simplification sizes from validated outputs, with missing outcomes omitted rather than zero. Black crosses identify independent outputs that violate the coupling bound; this is permitted for those baselines. Hatched protein slots lack an eligible prepared input. CPS-2F and CPS-3F have different fidelity feasible sets. Pair aliases map to the combined input manifest.')
    fig,axes=plt.subplots(2,2,figsize=(12,6),sharex=True,layout='constrained')
    for col,fidelity in enumerate(['continuous','discrete']):
        d=paired[paired.fidelity.eq(fidelity)].set_index('pair').loc[ordered];x=np.arange(len(d));base=METHODS[col];joint=METHODS[col+2]
        axes[0,col].axhline(1,color='black',ls='--');axes[0,col].scatter(x-.12,d.baseline_coupling_ratio,color=COLORS[base],s=20,label=base);axes[0,col].scatter(x+.12,d.joint_coupling_ratio,color=COLORS[joint],s=20,label=joint)
        axes[0,col].set_ylabel(r'$d_{dF}(A^\prime,B^\prime)/\delta_3$');axes[0,col].set_title(fidelity.capitalize()+' fidelity');axes[0,col].legend(fontsize=8)
        axes[1,col].axhline(0,color='.6',lw=.7);axes[1,col].scatter(x,d.delta_k,color=COLORS[joint],s=24)
        repair=d.equal_k_coupling_repair.to_numpy();axes[1,col].scatter(x[repair],d.delta_k.to_numpy()[repair],marker='*',color='gold',edgecolor='black',s=100,label='Equal k; coupling satisfied')
        axes[1,col].set_ylabel(r'$k_{\rm CPS}-k_{\rm independent}$');axes[1,col].set_xticks(x,[aliases[p] for p in ordered],rotation=90,fontsize=8)
        axes[1,col].legend(fontsize=8)
    savefig(fig,'coverage_coupling_cost','Does coupling improve feasibility at a vertex cost? Normalized output coupling (top) and objective change against the independent baseline with matching fidelity (bottom). Stars mark equal-k coupling repairs. Missing points are missing comparisons, not zero cost. Independent coupling violations are not baseline failures. Both curves use identical saved thresholds within each pair.')
    fig,axes=plt.subplots(2,2,figsize=(10.5,6.3),layout='constrained')
    for i,domain in enumerate(['protein','hurricane']):
        for r in f[f.domain.eq(domain)].itertuples():
            if pd.isna(r.solver_seconds):continue
            for ax,y in [(axes[i,0],r.solver_seconds),(axes[i,1],r.peak_rss_mib)]:
                if pd.isna(y):continue
                if r.validated:ax.scatter(r.nA+r.nB,y,s=40,marker='o',edgecolor=COLORS[r.method],facecolor='none' if r.origin=='reused' else COLORS[r.method],alpha=.75)
                else:ax.scatter(r.nA+r.nB,y,s=60,marker='x',color=COLORS[r.method])
        axes[i,0].set_yscale('log');axes[i,0].set_ylabel(domain.capitalize()+' solver time (s)');axes[i,1].set_ylabel('Sampled peak tree RSS (MiB)')
        for ax in axes[i]:ax.set_xlabel('Input vertices m+n')
    method_legend(axes[0,0]);axes[0,1].legend(handles=[Line2D([],[],marker='o',ls='',mfc='none',color='black',label='Inherited single observation'),Line2D([],[],marker='o',ls='',color='black',label='New validated observation'),Line2D([],[],marker='x',ls='',color='black',label='Time/memory at resource exit')],fontsize=8)
    savefig(fig,'coverage_runtime_memory','Time and sampled process-tree memory versus input vertex count. Hollow circles preserve inherited single observations; filled circles are new validated observations; crosses are resource exits, whose times are time to termination. Repeated timing measurements remain separate. Protein preparation regimes differ, so these descriptive points do not identify a causal size effect or complexity law. RSS includes runtime overhead.')
    eligible_order=[p for p in ordered if f[f.pair.eq(p)].eligible.any()]
    bound_figure(bounds,aliases,eligible_order)
    CAPTIONS['graph_size_levels']='Combined eligible single observations: input-only finite bound, active single-graph product H_A H_B, and discovered states D. Percentages name the denominator; neither ratio is progress. Crosses mark interrupted searches. Missing historical counters stay missing. Blocked inputs have no constructed graph and are excluded from this graph-count comparison. Finite bounds follow the saved derivation, not a hidden constant of one in asymptotic notation.'
    fig,axes=plt.subplots(1,3,figsize=(11,3.6),layout='constrained')
    for ax,quantity,label in zip(axes,['auxiliary','single_active','single_edges_with_waits'],['Auxiliary matching locations','Active single-curve vertices','Single-curve edges incl. waits']):
        for m in METHODS[2:]:
            for side in ['A','B']:
                g=bounds[bounds.method.eq(m)&bounds.quantity.eq(quantity+'_'+side)]
                x=g['nA' if side=='A' else 'nB'];good=g.observed.notna()
                done=g.construction_completed.eq(True)
                ax.scatter(x[good&done],g.observed[good&done],color=COLORS[m],marker='o' if side=='A' else '^',s=24,alpha=.65)
                ax.scatter(x[good&~done],g.observed[good&~done],color=COLORS[m],marker='x',s=40)
        ax.set_xlabel('Original vertices in that curve');ax.set_ylabel(label);ax.set_yscale('symlog',linthresh=1)
    axes[0].legend(handles=[Line2D([],[],color=COLORS[m],marker='o',ls='',label=m) for m in METHODS[2:]],fontsize=8)
    savefig(fig,'coverage_single_graphs','Supplement: each point describes one single-curve graph (circles A, triangles B), not a pooled configuration count. Auxiliary matching locations are not selectable output vertices. CPS-3F has no auxiliaries; zero is shown on a symmetric-log axis. Edges include waits after anchored pruning; the construction budget applies earlier during enumeration. Missing counters are omitted; crosses identify counts from interrupted construction where available.')
    (ext/'figure_captions.md').write_text('\n\n'.join('## '+n+'\n\n'+c for n,c in CAPTIONS.items()),encoding='utf-8')
    coverage_text(f,fresh,paired,input_frame,bounds,chapter,ext,aliases)
    OUT=chapter
    structural_check()
    caption_path=chapter/'figure_captions.md'
    old_captions=caption_path.read_text(encoding='utf-8').split('\n# Coverage extension figures')[0]
    caption_path.write_text(old_captions+'\n# Coverage extension figures\n\n'+(ext/'figure_captions.md').read_text(encoding='utf-8'),encoding='utf-8')
    dump(ext/'artifact_hashes.json',{str(p.relative_to(ext)):sha(p) for p in ext.rglob('*') if p.is_file() and p.name!='artifact_hashes.json' and 'cache' not in p.parts})

def coverage_text(f,fresh,paired,inputs,bounds,chapter,ext,aliases):
    import nbformat
    eligible=int(f.eligible.sum());observed=int(f.observed.sum());valid=int(f.validated.sum());resource=int(f.resource_limited.sum());blocked=int((~f.eligible).sum());pending=int(f.origin.eq('unexecuted').sum())
    continuous=paired[paired.fidelity.eq('continuous')];discrete=paired[paired.fidelity.eq('discrete')]
    numbers=dict(planned=len(f),eligible=eligible,observed=observed,validated=valid,resource_limited=resource,blocked=blocked,pending=pending,
        new_observed=int(fresh.origin.eq('new').sum()),new_validated=int(fresh.validated.sum()),new_resource_limited=int(fresh.resource_limited.sum()),
        continuous_matched=int(continuous.delta_k.notna().sum()),continuous_equal=int(continuous.delta_k.eq(0).sum()),continuous_extra=int(continuous.delta_k.gt(0).sum()),
        discrete_matched=int(discrete.delta_k.notna().sum()),discrete_equal=int(discrete.delta_k.eq(0).sum()),equal_k_repairs=int(continuous.equal_k_coupling_repair.sum()))
    dump(ext/'combined_counts.json',numbers)
    completion=[]
    for (domain,method),g in f.groupby(['domain','method'],sort=False):completion.append(dict(Domain=domain,Method=SHORT[method],Planned=len(g),Observed=int(g.observed.sum()),Validated=int(g.validated.sum()),Resource=int(g.resource_limited.sum()),Blocked=int((~g.eligible).sum()),Pending=int(g.origin.eq('unexecuted').sum())))
    ct=pd.DataFrame(completion);ct.to_csv(ext/'combined_completion_summary.csv',index=False)
    latex_table(ct,ext/'combined_completion.tex','Combined coverage: original and new single observations; repeated timings excluded. Preparation-blocked slots have no experimental outcome.','tab:coverage-completion')
    ip=inputs[['pair_alias','pair','eligible','nA','nB','coverage_cohort']].copy();ip.loc[~ip.eligible,['nA','nB']]=np.nan
    latex_table(ip,ext/'combined_inputs.tex','Pair aliases and frozen sizes. Blocked pairs have no frozen executable input; source lengths and missing-residue audits are retained separately.','tab:coverage-inputs')
    protein=inputs[inputs.domain.eq('protein')];hurricane=inputs[inputs.domain.eq('hurricane')]
    pr=ext/'protein_audit/identifiers.csv';audit=pd.read_csv(pr)
    paper_rows=[]
    for r in protein.to_dict('records'):
        b=r['pair'].split('__')[1];a=audit[audit.verified.eq(b)].iloc[0]
        paper_rows.append(dict(pair=r['pair'],printed_identifier=a.printed,verified_identifier=b,paper_length=a.reported_length,local_selected_CA_length=a.prepared_ca_length,deposited_CA_length=a.deposited_ca_length,fragments=a.fragments,model=a.model,author_chain=a.author_chain,label_chain=a.label_chain,eligible=r['eligible'],cohort=r['coverage_cohort'],interpretation=a.interpretation))
    pd.DataFrame(paper_rows).to_csv(ext/'protein_paper_coverage.csv',index=False)
    newh=hurricane[~hurricane.coverage_cohort.eq('original')]
    hlist='\n'.join(f"- {r.pair}: {int(r.nA)}/{int(r.nB)} vertices." for r in newh.itertuples())
    report=f'''# Coverage extension: saved results

## Old timing study
120/120 scheduled repetitions recorded: 117 validated, 3 transition-budget exits. The four previously pending repetitions completed successfully. The three CPS-2F AL072006__AL011900 exits are preserved; no extra timing rounds or replacement attempts were launched.

## Protein-paper coverage
Tables 1–3 of the local paper were extracted and checked (see provenance/paper_tables.json for PDF pages and hash). All 14 requested B identities are inventoried; {int(protein.eligible.sum())} pairs have eligible inputs. Two pairs, 3ntx.a and 3jq0.a, are preparation-blocked because the existing policy splits their internally broken chains. No cropping, invented residue matching, or new gap bridging was introduced.

The original five protein inputs remain unchanged. 4eca.b/c reuse supplied aligned w=16 arrays with documented residue omissions and 0.001-angstrom coordinate storage. The five eligible varying-length chains reuse existing US-align structural_fit on all selected observed residues before independent endpoint-preserving w=16 sampling. Their B input lengths vary from 15 to 29; A retains 22 vertices. This preparation differs from historical index-Kabsch inputs, so cross-pair changes cannot be attributed to size alone.

Identifier corrections: printed 107j.a maps to verified 1o7j.a; 4cea.b/d to 4eca.b/d; 1toh to author chain A. Paper 1d9q.d length 297 differs from local selected/deposited C-alpha length 325. Other length discrepancies, author versus label chains, alternate conformers, model selection and gaps are recorded in protein_paper_coverage.csv and protein_audit/. Identity coverage is not exact reproduction of the paper's coordinates, thresholds, or results.

## New hurricane selection
Earlier additions used size-quantile bands within the prior maximum-track-size range, unused-storm preference, and geometry diversity. The final20 cohort instead uses seeded, preference-constrained stratified random sampling from the entire eligible supplied-pair pool, without a size cap; see its separate report and selection_pool.json. No solver outcomes entered either selection procedure. All supplied vertices and the shared basin frame are retained; temporal synchronization is not inferred. Cohort identity is retained in the combined input and result tables. Unweighted final20 percentages describe that stratified cohort, not the archive.

{hlist}

## Completion
- New observations: {numbers['new_observed']}/{eligible-40} executable slots; {numbers['new_validated']} validated; {numbers['new_resource_limited']} resource exits; {blocked} preparation-blocked slots; {pending} pending executable slots.
- Combined: {observed} observed of {len(f)} intended slots ({eligible} eligible); {valid} validated; {resource} resource exits. Blocked/unexecuted slots are not observations.
- Original 40 single observations, new extension observations and 120 timing repetitions retain separate provenance. No selection of favorable attempts.

## Quality and computation
- Continuous matched comparisons: {numbers['continuous_equal']}/{numbers['continuous_matched']} have equal independent/CPS objectives; {numbers['continuous_extra']} have a higher coupled objective.
- Discrete matched comparisons: {numbers['discrete_equal']}/{numbers['discrete_matched']} have equal objectives.
- {numbers['equal_k_repairs']} continuous-fidelity cases satisfy coupling at unchanged k where the saved independent output violates it. Such a violation is permitted for independent simplification.
- Missing solution sizes stay missing. Computational plots distinguish time to solution from time to resource termination. New runs are single observations, not timing repetitions; inherited times are marked separately.
- Graph figures show finite input bounds, active products and discovered states. They do not claim that the full graph was enumerated. Occupancy is not progress; transition checks are not unique edges. Source formulas and missing-counter reasons remain in graph_bounds.csv.

## Limitations and stopping point
All protein pairs share 1o7j.a, with further reuse of 4eca chains. Subsampling and differing alignment regimes limit biological independence and size comparisons. Alpha=0.5 was pilot-informed. Shared external time/memory budgets coexist with documented method-specific internal limits. Validation proves feasibility, not independent optimality. Resource exits do not prove infeasibility. No further inputs, thresholds, budget increases or timing repetitions were run.

## Reproduction
Frozen configuration_manifest.json / run_manifest.csv and input_manifest.json; raw outputs and phase/resource logs in slots/. Resume only genuinely unexecuted slots with the existing coverage_extension.py execute command after readiness passes. Reporting uses paper_publication.py coverage_report and performs no optimization.
'''
    active_cohort=ROOT/'output/paper_v1/hurricane_extension_3'
    if (ROOT/'output/paper_v1/hurricane_final_20/configuration_manifest.json').exists():active_cohort=ROOT/'output/paper_v1/hurricane_final_20'
    if not (active_cohort/'configuration_manifest.json').exists():active_cohort=ext
    progress=read(active_cohort/'progress.json') if (active_cohort/'progress.json').exists() else {}
    pause=''
    if progress.get('phase')=='paused':
        pause=('\n## Execution pause\n\nLast update: '+str(progress.get('utc'))+'. '+str(progress.get('pause_reason'))+
               '\n\nNext pending configuration: '+str(progress.get('pair'))+', '+str(progress.get('method'))+
               '. Readiness evidence is retained in the slot readiness JSON files. No solver attempt was launched for the failed readiness check.'+
               '\n\nResume command, after competing activity subsides:\n```powershell\n.venv/Scripts/python.exe notebooks_paper/coverage_extension.py execute --cohort '+active_cohort.name+'\n```\n')
        report+=pause
    (ext/'report.md').write_text(report,encoding='utf-8');(ext/'status.md').write_text(report.split('## Quality')[0]+pause,encoding='utf-8')
    # Preserve the completed original chapter before adding the extension section.
    archive=chapter/'provenance/experiments_before_coverage_extension.tex'
    if not archive.exists():shutil.copyfile(chapter/'experiments.tex',archive)
    base=archive.read_text(encoding='utf-8')
    replacements={'completion':'coverage_completion','coupling_vertex_cost':'coverage_coupling_cost','graph_size_levels':'graph_size_levels'}
    def retain_figure(match):
        return '' if any(r'\label{fig:'+name+'}' in match.group(0) for name in replacements) else match.group(0)
    base=re.sub(r'\\begin\{figure\}\[p\].*?\\end\{figure\}',retain_figure,base,flags=re.S)
    for before,after in replacements.items():base=base.replace(r'\ref{fig:'+before+'}',r'\ref{fig:extension-'+after+'}')
    base=base.replace('The completed primary cohort and the fresh timing study are distinct evidence sets.',
        'The completed original cohort, its fresh timing study, and the subsequent coverage extension are distinct evidence sets. Combined coverage and graph figures include the extension described in Section~\\ref{sec:coverage-extension}; original-cohort counts below retain their original denominators.')
    section=r'\section{Coverage extension at fixed $\alpha=0.5$}'+'\n'+r'\label{sec:coverage-extension}'+'\n'
    section+=f'The completed timing study contains 120 recorded repetitions: 117 validated outputs and three retained transition-budget exits. The four previously pending repetitions were completed without repeating other slots. The following extension consists of single observations, not additional timing rounds.\n\n'
    section+=f'The combined intended coverage is {len(protein)} protein pairs and {len(hurricane)} hurricane pairs, or {len(f)} method configurations. Of these, {eligible} have eligible inputs and {observed} have recorded observations; {valid} outputs passed validation, {resource} terminated at a resource limit, {blocked} are preparation-blocked, and {pending} eligible slots remain pending. Table~\\ref{{tab:coverage-completion}} exposes the denominators.\n\n'
    section+=r'\input{../coverage_extension/combined_completion.tex}'+'\n'+r'\input{../coverage_extension/combined_inputs.tex}'+'\n'
    section+=r'''The paper's Tables 1--3 were checked against the local PDF. The known printed identifiers \texttt{107j.a} and \texttt{4cea.b/d} are reconciled to \texttt{1o7j.a} and \texttt{4eca.b/d}; \texttt{1toh} denotes author chain A. In particular, the reported \texttt{1d9q.d} length of 297 differs from the 325 selected local C-alpha coordinates. Full identifier, residue, alternate-conformer, model and chain-break audits accompany the results. Coverage of identities does not reproduce the paper's original coordinate preparation, thresholds or results.

The original protein inputs were not rebuilt. The supplied aligned \texttt{4eca.b/c} inputs were reused, including their documented omissions. Five eligible varying-length additions use the existing rigid sequential US-align preparation on their entire observed chains, without guessed same-index residue correspondences. Endpoint-preserving stride sampling is applied independently with $w=16$: A has 22 vertices and the new B curves have 15--29 vertices. Unlike the original cohort, protein sizes are not all equal, but shared structure and different historical alignment regimes still prevent general scaling or independent biological-replication claims. The internally broken \texttt{3ntx.a} and \texttt{3jq0.a} chains are blocked: the available preparation splits fragments, and neither cropping nor a new whole-chain gap policy is introduced.

Earlier additional hurricane pairs used size-quantile coverage within the prior range, unused-storm preference and input-geometric diversity. The separately reported final stratified cohort uses a seeded random procedure and no input-size cap. Every supplied vertex and the shared basin frame is retained. No temporal synchronization follows from coordinates alone. Pair aliases and exact input hashes link the figures to the frozen manifest; full-precision thresholds retain the previous rule and all resource limits are unchanged.

'''
    section+=f"Among validated matched comparisons, continuous fidelity gives equal objectives in {numbers['continuous_equal']} of {numbers['continuous_matched']} cases and a larger coupled objective in {numbers['continuous_extra']}. Discrete fidelity gives equal objectives in {numbers['discrete_equal']} of {numbers['discrete_matched']} cases. There are {numbers['equal_k_repairs']} equal-objective continuous-fidelity coupling repairs. Independent coupling violations remain descriptive, not baseline failures. These counts include inherited single observations and never count timing repetitions as independent inputs.\n\n"
    figure_statements={'coverage_completion':'Completion includes the blocked and resource-limited cases rather than selecting successful pairs.',
        'coverage_quality':'The combined solution-size plot retains matched inputs and marks independent coupling violations.',
        'coverage_coupling_cost':'Normalized coupling and matched-fidelity objective differences isolate the coupling requirement.',
        'coverage_runtime_memory':'Time and memory versus input size are descriptive; resource-exit times and inherited single observations are distinguished.',
        'graph_size_levels':'Discovered states remain distinct from the finite input bound and active-graph product; interrupted counts are not final graph enumeration.',
        'coverage_single_graphs':'Auxiliary matching locations and single-curve graphs are different counting quantities, shown separately.'}
    for name,caption in CAPTIONS.items():
        section+=figure_statements[name]+r' See Figure~\ref{fig:extension-'+name+'}.\n'
        section+='\n'.join([r'\begin{figure}[p]',r'\centering',r'\includegraphics[width=\textwidth]{../coverage_extension/figures/'+name+'.pdf}',r'\caption{'+tex(caption)+'}',r'\label{fig:extension-'+name+'}',r'\end{figure}'])+'\n'
    section+='\n'+r'The original Atlantic and Pacific geometric illustrations are retained (Figure~\ref{fig:geometric_examples}); no new case was selected for an attractive outcome.'+'\n\n'
    section+=r'''The accompanying saved-results gallery provides a geometric view for every available attempt, with source coordinates, selected original vertices, endpoints, thresholds and validation status. Protein panels use static equal-aspect 3D views; hurricane panels retain the common projected basin frame. Input-only panels identify attempts without validated outputs. Identical timing-repeat geometries may share a figure, but timing outcomes remain separate. Figure generation is outside solver and validation measurements. The per-attempt metric table links raw checkpoints, validation recoveries, phase/resource evidence and finite-bound definitions. Missing historical counters remain missing; transition checks do not substitute for unique graph edges or dynamic-programming labels.

'''
    base=base.replace('Every protein curve has 22 vertices','Every protein curve in the original cohort has 22 vertices')
    base=base.replace('Shared protein structure and identical protein sizes prevent claims of independent biological replication or protein scaling.', 'Shared protein structure prevents claims of independent biological replication. Original protein sizes were identical; the extension varies sizes but also preparation regimes, preventing a controlled scaling claim.')
    base=base.replace(r'\section{Limitations}',section+r'\section{Limitations}',1)
    (chapter/'experiments.tex').write_text(base,encoding='utf-8')
    trace=[]
    for k,v in numbers.items():trace.append(dict(claim=k,value=v,source='combined_results.csv' if k not in ['continuous_matched','continuous_equal','continuous_extra','discrete_matched','discrete_equal','equal_k_repairs'] else 'matched_comparisons.csv'))
    pd.DataFrame(trace).to_csv(ext/'claim_traceability.csv',index=False)
    (chapter/'coverage_extension.md').write_text(report,encoding='utf-8')
    original_findings=(chapter/'findings.md').read_text(encoding='utf-8').split('\n## Coverage extension update')[0]
    (chapter/'findings.md').write_text(original_findings+f'\n## Coverage extension update\n\n{valid}/{observed} recorded combined single observations validated; {resource} resource exits; {blocked} preparation-blocked slots. See [extension report](../coverage_extension/report.md) for all {len(f)} intended slots and preparation caveats.\n',encoding='utf-8')
    notebook=ROOT/'notebooks_paper/expanded_comparison_saved_results.ipynb';nb=nbformat.read(notebook,4)
    # Preserve the original examples/timing; replace now-redundant original-cohort quality figures in the reader.
    remove=['coupling_vertex_cost','graph_size_levels','completion.png','work_indicators','supplement_memory_phases']
    kept=[]
    for c in nb.cells:
        if c.metadata.get('coverage_extension'):continue
        if c.cell_type=='code' and any(x in c.source for x in remove):
            if kept and kept[-1].cell_type=='markdown' and kept[-1].metadata.get('publication_saved_results'):kept.pop()
            continue
        kept.append(c)
    def add(kind,source):
        c=nbformat.v4.new_markdown_cell(source) if kind=='md' else nbformat.v4.new_code_cell(source);c.metadata['coverage_extension']=True;kept.append(c)
    add('md','# Coverage extension: fixed protocol, saved evidence\n\nOriginal cohort, fresh timing repetitions and new single observations remain distinct. Preparation-blocked slots are not observations.')
    add('code',"EXTENSION = PROJECT/'output/paper_v1/coverage_extension'\ndisplay(Markdown((EXTENSION/'report.md').read_text(encoding='utf-8')))\ndisplay(pd.read_csv(EXTENSION/'combined_completion_summary.csv'))")
    for name,caption in CAPTIONS.items():
        add('md','## '+name.replace('_',' ').capitalize()+'\n\n'+caption)
        add('code',f"display(Image(filename=str(EXTENSION/'figures/{name}.png')))")
    add('md','Full data: [combined results](../output/paper_v1/coverage_extension/combined_results.csv), [paper identity audit](../output/paper_v1/coverage_extension/protein_paper_coverage.csv), [graph bounds](../output/paper_v1/coverage_extension/graph_bounds.csv), [per-attempt geometric gallery](../output/paper_v1/coverage_extension/attempt_gallery/index.html). No cell launches optimization.')
    nb.cells=kept;nbformat.write(nb,notebook)
    audit_coverage(f,fresh,ext)
    added=ROOT/'output/paper_v1/hurricane_extension_3'
    if (added/'configuration_manifest.json').exists():
        newpairs={r['pair'] for r in read(added/'input_manifest.json')['inputs']}
        subset=f[f.pair.isin(newpairs)]
        audit_coverage(f,subset,added)
        summary=('# Three-pair hurricane extension\n\n'+subset[['pair','method','status','k','solver_seconds','validation_seconds','peak_rss_mib']].to_csv(index=False)+
                 '\nAll observations are single attempts under the frozen protocol; no timing repetitions. See ../coverage_extension/report.md and ../coverage_extension/attempt_gallery/index.html for combined comparisons and per-attempt evidence.\n'+pause)
        (added/'status.md').write_text(summary,encoding='utf-8')
    # Structural validation must also read the extension's generated LaTeX tables.
    dump(ext/'report_provenance.json',dict(utc=datetime.now(timezone.utc).isoformat(),source_hashes={str(p.relative_to(ROOT)):sha(p) for p in [SOURCE/'results.csv',TIMING/'results.csv',ext/'configuration_manifest.json',ext/'results.csv',Path(__file__)]},optimization_executed=False))
    render_qa()
    print('Coverage report saved:',numbers,flush=True)

def audit_coverage(combined,fresh,ext):
    m=read(ext/'configuration_manifest.json');checks=[]
    for r in m['rows']:
        if not r['eligible']:continue
        c=r['configuration'];assert sha(ROOT/c['input_file'])==c['input_sha256'] and c['certificate'] is False
        for p,h in c['algorithms'].items():assert sha(ROOT/p)==h
        slot=ext/'slots'/f"{r['plan_order']:02d}";attempts=list((slot/'attempts').glob('*/configuration.json'));assert len(attempts)<=1
        x=fresh[fresh.configuration_id.eq(r['configuration_id'])].iloc[0]
        if x.validated:
            d=ROOT/x.attempt_path;raw=read(d/'raw_solver.json');val=read(d/'validation.json')
            assert val['status']=='validated' and val['raw_checkpoint_sha256']==sha(d/'raw_solver.json') and raw['unix']<=val['unix']
            assert val['k']==x.k and val['kA']==x.kA and val['kB']==x.kB
            if r['method'].startswith('CPS'):assert raw['certificate_enabled'] is False and raw['solver_status']=='optimal_configuration_graph' and raw['observation']['configuration_graph_entered']
        else:assert pd.isna(x.k)
        checks.append(dict(configuration_id=r['configuration_id'],attempt_count=len(attempts),status=x.status,validated=x.validated,hashes_match=True))
    pd.DataFrame(checks).to_csv(ext/'execution_audit.csv',index=False)
    assert combined.configuration_id.is_unique

def attempt_gallery():
    """Consolidate saved checkpoints only. Never solve, validate anew, or generate auxiliaries."""
    import os,time,html,sys
    base=ROOT/'output/paper_v1';target=base/'coverage_extension/attempt_gallery'
    target.mkdir(exist_ok=True);(target/'figures').mkdir(exist_ok=True);(target/'metrics').mkdir(exist_ok=True)
    def optional(p):return read(p) if p.exists() else {}
    def rel(p):return str(Path(p).relative_to(ROOT)).replace('\\','/')
    def link(p):return os.path.relpath(p,target).replace('\\','/')
    def atomic_progress(x):
        x.update(utc=datetime.now(timezone.utc).isoformat(),optimization_executed=False)
        temp=target/'progress.tmp';dump(temp,x);os.replace(temp,target/'progress.json')
    # Paths, rather than configuration IDs, distinguish repetitions and diagnostics.
    folders=sorted({p.parent for name in ['raw_solver.json','solver_exit.json','supervisor.json'] for p in base.rglob(name) if 'attempts' in p.parts})
    rows=[];excluded=[];recoveries={};uses={}
    for source in [SOURCE/'results.csv',TIMING/'results.csv',base/'coverage_extension/results.csv',base/'hurricane_extension_3/results.csv',base/'hurricane_final_20/results.csv',base/'unified_pilot_results/results.csv']:
        if source.exists():
            for x in pd.read_csv(source).to_dict('records'):
                p=x.get('attempt_path')
                if isinstance(p,str) and p:uses.setdefault(str((ROOT/p).resolve()),[]).append(dict(table=rel(source),origin=x.get('origin'),repetition=x.get('repetition'),configuration_id=x.get('configuration_id')))
    for d in folders:
        c=optional(d/'configuration.json')
        if c.get('validation_replay_raw'):
            v=optional(d/'validation.json')
            recoveries.setdefault(str((ROOT/c['validation_replay_raw']).resolve()),[]).append((d,v))
    for d in folders:
        c=optional(d/'configuration.json')
        if not c or c.get('test_mode') or c.get('validation_replay_raw'):
            excluded.append(dict(path=rel(d),reason='validation recovery (linked to original)' if c.get('validation_replay_raw') else 'regression/test or missing configuration'));continue
        raw=optional(d/'raw_solver.json');payload=raw.get('raw_result',{});v=optional(d/'validation.json');vpath=d/'validation.json'
        original_validation=v.copy();recovery_records=recoveries.get(str((d/'raw_solver.json').resolve()),[])
        for rd,rv in sorted(recovery_records,key=lambda z:str(z[0])):
            if rv.get('status')=='validated':v=rv;vpath=rd/'validation.json';break
        sup=optional(d/'supervisor.json');exit_record=optional(d/'solver_exit.json');measure=optional(d/'measurements.json');warm=optional(d/'warmup.json')
        events=[json.loads(line) for line in (d/'phases.jsonl').read_text().splitlines()] if (d/'phases.jsonl').exists() else []
        obs={}
        for e in events:
            if 'solver_subphase' in e:obs.update(e)
        obs.update(raw.get('observation',{}));obs.update(exit_record.get('counters',{}))
        valid=v.get('status')=='validated'
        if valid and v.get('raw_checkpoint_sha256') and v['raw_checkpoint_sha256']!=sha(d/'raw_solver.json'):raise AssertionError(('Validation checkpoint mismatch',d))
        status='validated' if valid else exit_record.get('status') or v.get('status') or sup.get('status','unavailable')
        if status=='worker_exited':status='validation_incomplete' if raw else 'no_saved_solution'
        ident=hashlib.sha256(rel(d).encode()).hexdigest()[:18]
        selectedA=payload.get('indices_A');selectedB=payload.get('indices_B')
        r={k:c.get(k) for k in ['pair','domain','method','alpha','delta1','delta2','delta3','units','nA','nB','input_file','input_sha256','certificate','configuration_id']}
        r.update(attempt_id=ident,attempt_path=rel(d),plan_order=len(rows)+1,origin='historical attempt',status=status,validated=valid,
            evidence_cohort=rel(d).split('/')[3],observation_type='timing repetition' if 'timing_study' in d.parts else 'single observation / historical diagnostic',
            usage=json.dumps(uses.get(str(d.resolve()),[])),configuration_sha256=sha(d/'configuration.json'),code_provenance=json.dumps(c.get('code_version',raw.get('code_version',{}))),algorithm_hashes=json.dumps(c.get('algorithms',{})),
            selected_indices_A=json.dumps(selectedA) if selectedA is not None else None,selected_indices_B=json.dumps(selectedB) if selectedB is not None else None,
            kA=len(selectedA) if valid and selectedA is not None else None,kB=len(selectedB) if valid and selectedB is not None else None,k=v.get('k') if valid else None,
            raw_objective=payload.get('k'),validation_status=v.get('status','unavailable'),original_validation_status=original_validation.get('status','unavailable'),validation_recovered=vpath.parent!=d,
            coupling_distance=v.get('coupling_distance',payload.get('dDF')),coupling_feasible=v.get('coupling_feasible'),coupling_required=str(c.get('method','')).startswith('CPS'),
            solver_status=raw.get('solver_status'),solution_route=raw.get('solution_route'),configuration_graph_entered=obs.get('configuration_graph_entered'),
            search_completed=raw.get('solver_status')=='optimal_configuration_graph',states=obs.get('states_discovered',payload.get('states')),states_dequeued=obs.get('states_dequeued'),
            processed_states_lower_bound=obs.get('completed_states_lower_bound'),queue_size=obs.get('queue_size'),transitions=obs.get('product_transitions',payload.get('transitions')),
            chain_graphs_completed=obs.get('chain_graphs_completed'),counter_timestamp=obs.get('utc',raw.get('utc')),
            unique_configuration_edges=None,dynamic_programming_labels=None,
            warmup_seconds=warm.get('seconds'),solver_seconds=raw.get('solver_seconds',exit_record.get('solver_seconds')),validation_seconds=v.get('validation_seconds'),optional_measurement_seconds=measure.get('seconds'),
            total_process_seconds=sup.get('process_seconds'),peak_rss_mib=sup.get('approximate_sampled_tree_rss_peak_mib'),termination_phase=sup.get('limiting_phase') or ('solver' if exit_record else None),
            termination_reason=exit_record.get('failure_reason') or (status if not valid else None),
            raw_result_path=rel(d/'raw_solver.json') if raw else None,validation_path=rel(vpath) if vpath.exists() else None,resource_limits=json.dumps(dict(graph=c.get('graph_limits'),phases=c.get('phase_budgets'),memory_mib=c.get('memory_mib'))),
            measurement_status=measure.get('status','unavailable'),fidelity_evidence=json.dumps(v.get('checks',{})),fidelity_tolerances='Historical tolerance provenance unavailable in this schema; inspect linked code hashes')
        versions=c.get('code_version',raw.get('code_version',{}))
        validator=ROOT/'notebooks_paper/graph_only_v3/common.py'
        if any(str(k).replace('\\','/')=='notebooks_paper/graph_only_v3/common.py' and h==sha(validator) for k,h in versions.items()):
            r['fidelity_tolerances']='Discrete fidelity and coupling: delta + 1e-9; coupling path cost check: abs difference < 1e-7; continuous decision unchanged at delta, geometric TOL=1e-10 (algorithm hash linked)'
        r['external_limiting_phase_elapsed_seconds']=sup.get('phase_elapsed_seconds') if sup.get('limiting_phase') else None
        r['memory_scope']=sup.get('memory_scope','Historical monitor scope unavailable')
        r['solver_internal_checks']=raw.get('solver_internal_validation')
        r['validation_recovery_process_seconds']=optional(vpath.parent/'supervisor.json').get('process_seconds') if vpath.parent!=d else None
        if 'timing_study' in d.parts:r['observation_type']='timing repetition'
        elif 'hurricane_final_20' in d.parts:r['observation_type']='stratified final20 single observation'
        elif any(name in d.parts for name in ['coverage_extension','hurricane_extension_3']):r['observation_type']='new single observation'
        elif uses.get(str(d.resolve())):r['observation_type']='reused cohort observation'
        else:r['observation_type']='historical diagnostic'
        for side in ['A','B']:
            for endpoint in ['lower','upper']:r[f'continuous_{side}_{endpoint}']=measure.get(f'continuous_{side}_{endpoint}')
            g=obs.get('single_graph_'+side,{})
            for field in ['auxiliary','extended_vertices','dense_states','active_states','edges_including_waits','construction_completed']:
                value=g.get(field)
                if field=='auxiliary' and value is None:value=payload.get('auxiliary_'+side)
                r[f'{field}_{side}']=value
            r['retained_fraction_'+side]=r['k'+side]/r['n'+side] if valid and r['k'+side] is not None and r['n'+side] else None
        r['retained_fraction_total']=(r['kA']+r['kB'])/(r['nA']+r['nB']) if valid and r['kA'] is not None and r['kB'] is not None else None
        r['edge_count_convention']='Active/pruned outgoing entries including one wait per active vertex' if any(r['edges_including_waits_'+s] is not None for s in ['A','B']) else 'unavailable'
        r['count_scope']='completed search (not full enumeration)' if r['search_completed'] else 'observed before termination / phase not reached'
        r['metric_availability']={k:('not_applicable' if not r['coupling_required'] and k in ['states','states_dequeued','queue_size','transitions','unique_configuration_edges','dynamic_programming_labels'] else 'missing' if value is None else 'available') for k,value in r.items()}
        evidence=dict(configuration=c,raw_checkpoint=raw,validation=v,original_validation=original_validation,validation_recoveries=[dict(path=rel(rd),record=rv) for rd,rv in recovery_records],solver_exit=exit_record,supervisor=sup,measurements=measure,observation=obs)
        dump(target/'metrics'/f'{ident}.json',dict(metrics=r,evidence=evidence));rows.append(r)
    # Legacy workers stored integrated validation and timing in *.result.json.
    # Preserve their weaker phase provenance, rather than inventing raw checkpoints.
    for config_path in sorted(base.glob('**/workers/*.json')):
        c=optional(config_path)
        if not isinstance(c.get('param'),dict):continue
        result_path=Path(c.get('result_path',str(config_path.with_name(config_path.stem+'.result.json'))))
        p=c['param'];saved=optional(result_path);evidence_path=result_path if result_path.exists() else config_path
        legacy_table=None
        for candidate in [config_path.parent.parent/name for name in ['pilot_results.csv','protein_results.csv','hurricane_results.csv']]:
            if not candidate.exists():continue
            table=pd.read_csv(candidate)
            if 'case' not in table or 'method' not in table:continue
            match=table[table['case'].eq(p.get('case')) & table.method.eq(c.get('method'))]
            if len(match)==1:
                saved={**match.replace({np.nan:None}).iloc[0].to_dict(),**saved};legacy_table=candidate;break
        ident=hashlib.sha256(rel(evidence_path).encode()).hexdigest()[:18]
        r={key:None for key in rows[0]} if rows else {}
        r.update({key:p.get(key) for key in ['domain','pair','alpha','delta1','delta2','delta3','nA','nB','units','input_file']})
        input_path=ROOT/(p.get('input_file') or '__missing__');valid=saved.get('constraints_valid') is True
        r.update(attempt_id=ident,attempt_path=rel(evidence_path),configuration_id=p.get('case','')+'|'+c.get('method',''),plan_order=len(rows)+1,
            method=c.get('method'),origin='historical legacy attempt',observation_type='legacy integrated-validation observation',evidence_cohort=rel(result_path).split('/')[3],
            status='legacy_validated' if valid else saved.get('status','unavailable'),validated=valid,validation_status='legacy integrated constraints_valid' if valid else 'unavailable',
            input_sha256=sha(input_path) if input_path.exists() else None,input_hash_provenance='Calculated during evidence audit; historical run hash not asserted',
            configuration_sha256=sha(config_path),code_provenance='Unavailable in worker result; consult historical cohort provenance',algorithm_hashes=None,
            selected_indices_A=json.dumps(saved['indices_A']) if 'indices_A' in saved else None,selected_indices_B=json.dumps(saved['indices_B']) if 'indices_B' in saved else None,
            kA=saved.get('kA') if valid else None,kB=saved.get('kB') if valid else None,k=saved.get('k') if valid else None,
            raw_objective=saved.get('k'),certificate=saved.get('certificate_enabled'),coupling_distance=saved.get('coupling_distance'),coupling_feasible=saved.get('coupling_feasible'),coupling_required=c.get('method','').startswith('CPS'),
            solver_status=saved.get('raw_status'),solution_route=saved.get('route'),configuration_graph_entered=saved.get('configuration_graph_entered'),search_completed=saved.get('raw_status')=='optimal_configuration_graph',
            states=saved.get('states_discovered',saved.get('states')),transitions=saved.get('transitions_examined',saved.get('transitions')),chain_graphs_completed=saved.get('chain_graphs_completed'),
            warmup_seconds=saved.get('warmup_seconds'),solver_seconds=saved.get('solve_seconds'),total_process_seconds=saved.get('process_wall_seconds'),peak_rss_mib=saved.get('peak_process_rss_mib'),termination_phase=saved.get('phase') if not valid else None,termination_reason=saved.get('failure_reason',saved.get('status')) if not valid else None,
            validation_path=rel(result_path) if valid else None,legacy_result_path=rel(result_path),fidelity_evidence='Integrated constraints_valid flag; no separate external-validation timing/checkpoint',
            fidelity_tolerances='Historical hashed code/protocol; not inferred from present implementation',edge_count_convention='unavailable',count_scope='legacy saved counters; full enumeration not asserted')
        r['memory_scope']='Legacy peak_process_rss_mib field; process-tree scope not asserted by this adapter'
        for side in ['A','B']:
            r['auxiliary_'+side]=saved.get('auxiliary_'+side)
            r['discrete_fidelity_'+side]=saved.get('discrete_fidelity_'+side)
            for endpoint in ['lower','upper']:r[f'continuous_{side}_{endpoint}']=saved.get(f'continuous_{side}_{endpoint}')
            r['retained_fraction_'+side]=r['k'+side]/r['n'+side] if r['k'+side] is not None else None
        r['retained_fraction_total']=(r['kA']+r['kB'])/(r['nA']+r['nB']) if valid else None
        r['metric_availability']={k:'missing' if value is None else 'available' for k,value in r.items()}
        dump(target/'metrics'/f'{ident}.json',dict(metrics=r,evidence=dict(configuration=c,legacy_result=saved,legacy_summary_source=rel(legacy_table) if legacy_table else None)));rows.append(r)
    f=pd.DataFrame(rows);pd.DataFrame(excluded).to_csv(target/'excluded_or_recovery_records.csv',index=False)
    # Reuse the established finite-bound derivation without optimizer execution.
    sys.path.insert(0,str(ROOT/'notebooks_paper/graph_only_v3'))
    import expanded_comparison as shared
    previous=shared.OUT;shared.OUT=target
    try:bounds=shared.graph_tables(f[f.method.isin(['CPS-2F','CPS-3F']) & f.nA.notna() & f.nB.notna()])
    finally:shared.OUT=previous
    # Legacy result files are not raw_solver.json folders. Carry their saved
    # auxiliary counts into the same finite formulas, without rebuilding graphs.
    bypath={r['attempt_path']:r for r in rows}
    for i,b in bounds.iterrows():
        r=bypath[b.attempt_path];q=b.quantity
        count=None;denominator=b.finite_upper_bound
        for side in ['A','B']:
            aux=r.get('auxiliary_'+side);n=r.get('n'+side)
            if aux is not None and n is not None:
                if q=='auxiliary_'+side:count=aux
                if q=='single_dense_'+side:count=n*(n+aux)
                if q=='single_active_'+side and pd.isna(denominator):denominator=n*(n+aux)
        if q=='configuration_states_extended_bound' and all(r.get('auxiliary_'+s) is not None for s in ['A','B']):
            denominator=r['nA']*r['nB']*(r['nA']+r['auxiliary_A'])*(r['nB']+r['auxiliary_B'])
        if count is not None and pd.isna(b.observed):bounds.at[i,'observed']=count;bounds.at[i,'missing_reason']=''
        if pd.notna(denominator):bounds.at[i,'finite_upper_bound']=denominator
        observed=bounds.at[i,'observed']
        if pd.notna(observed) and pd.notna(denominator) and denominator>0:bounds.at[i,'percent_of_bound']=100*observed/denominator
    bounds.to_csv(target/'per_attempt_bounds.csv',index=False)
    for r in rows:
        r['bounds']=bounds[bounds.attempt_path.eq(r['attempt_path'])].replace({np.nan:None}).to_dict('records')
        for b in r['bounds']:
            if b['quantity'] in ['configuration_states_input_bound','configuration_states_extended_bound','configuration_states_active_bound','candidate_transition_checks']:
                r[b['quantity']+'_upper']=b['finite_upper_bound'];r[b['quantity']+'_percent']=b['percent_of_bound'];r[b['quantity']+'_formula']=b['formula']
        r['bounds_table_path']=rel(target/'per_attempt_bounds.csv')
    # Rendering occurs only after measured execution has stopped.
    completed=0;last=time.monotonic();cards=[];canonical=set()
    for r in rows:
        p=ROOT/(r.get('input_file') or '__missing__');key=None;figure=None;figure_note='Input unavailable in saved configuration'
        if p.exists() and p.suffix=='.npz':
            if r.get('input_sha256') and sha(p)!=r['input_sha256']:raise AssertionError(('Input hash mismatch',p))
            arrays=np.load(p);A,B=arrays['A'],arrays['B'];valid=r['validated'] and r['selected_indices_A'] is not None and r['selected_indices_B'] is not None
            ia=json.loads(r['selected_indices_A']) if valid else [];ib=json.loads(r['selected_indices_B']) if valid else []
            if valid:
                for ids,P in [(ia,A),(ib,B)]:assert ids[0]==0 and ids[-1]==len(P)-1 and all(a<b for a,b in zip(ids,ids[1:]))
            fingerprint=dict(input=sha(p),indices=[ia,ib],method=r['method'],thresholds=[r['delta1'],r['delta2'],r['delta3']],status=r['status'],coupling=r['coupling_distance'])
            key=hashlib.sha256(json.dumps(fingerprint,sort_keys=True).encode()).hexdigest()[:22];figure=target/'figures'/key
            canonical.add(key)
            if not figure.with_suffix('.png').exists() or not figure.with_suffix('.pdf').exists():
                three=A.shape[1]==3;fig=plt.figure(figsize=(9,7.4));ax=fig.add_subplot(111,projection='3d' if three else None)
                for P,indices,col,label in [(A,ia,'#2166ac','A'),(B,ib,'#b2182b','B')]:
                    ax.plot(*P.T,color=col,alpha=.35,lw=1.2,label=label+' original')
                    if valid:
                        Q=P[indices];ax.plot(*Q.T,color=col,lw=2,marker='o',ms=3.5,label=label+' simplified / selected')
                    ax.scatter(*P[0],marker='s',s=65,color=col,edgecolors='black',label=label+' start')
                    ax.scatter(*P[-1],marker='*',s=115,color=col,edgecolors='black',label=label+' end')
                allpoints=np.vstack([A,B]);mid=(allpoints.min(0)+allpoints.max(0))/2;span=max(np.ptp(allpoints,axis=0).max()/2,1e-9)*1.07
                ax.set_xlim(mid[0]-span,mid[0]+span);ax.set_ylim(mid[1]-span,mid[1]+span)
                if three:ax.set_zlim(mid[2]-span,mid[2]+span);ax.set_box_aspect((1,1,1));ax.set_zlabel('z ('+str(r['units'])+')')
                else:ax.set_aspect('equal',adjustable='box')
                ax.set_xlabel('x ('+str(r['units'])+')');ax.set_ylabel('y ('+str(r['units'])+')')
                output=f"output {r['kA']}/{r['kB']}, k={r['k']}" if valid else 'INPUT ONLY; no externally validated output shown'
                ax.set_title(f"{r['pair']} | {r['method']}\ninput {len(A)}/{len(B)}; {output}\nalpha={r['alpha']}; delta=({r['delta1']:.6g}, {r['delta2']:.6g}, {r['delta3']:.6g}) {r['units']}",fontsize=10,pad=18)
                note=f"{r['status']} | "+('Static 3D view; distances evaluated in 3D.' if three else 'Common basin projected coordinates; all supplied vertices.')
                if valid and r['coupling_feasible'] is False:note+=f"\nCOUPLING {r['coupling_distance']:.6g} > delta3 {r['delta3']:.6g} — valid independent solution."
                elif valid and r['coupling_distance'] is not None:note+=f"\nDiscrete output coupling {r['coupling_distance']:.6g}; threshold {r['delta3']:.6g}."
                ax.legend(fontsize=7,loc='upper left');fig.text(.5,.035,note,ha='center',fontsize=9,color='#9c2020' if valid and r['coupling_feasible'] is False else 'black')
                fig.subplots_adjust(top=.79,bottom=.14);fig.savefig(figure.with_suffix('.png'),dpi=140,facecolor='white');fig.savefig(figure.with_suffix('.pdf'),facecolor='white');plt.close(fig)
            figure_note='Saved original coordinates and selected indices; exact thresholds in linked metrics'
        r['figure_png']=rel(figure.with_suffix('.png')) if figure else None;r['figure_pdf']=rel(figure.with_suffix('.pdf')) if figure else None;r['figure_note']=figure_note
        metric=target/'metrics'/f"{r['attempt_id']}.json";detail=read(metric);detail['metrics']=r;dump(metric,detail)
        esc=html.escape;links=[f'<a href="{link(metric)}">metrics + bounds + provenance</a>',f'<a href="{link(ROOT/r["attempt_path"])}">attempt files</a>']
        for label,field in [('raw result','raw_result_path'),('validation','validation_path'),('PDF','figure_pdf')]:
            if r.get(field):links.append(f'<a href="{link(ROOT/r[field])}">{label}</a>')
        image=f'<a href="{link(figure.with_suffix(".png"))}"><img loading="lazy" src="{link(figure.with_suffix(".png"))}" alt="{esc(r["pair"])}"></a>' if figure else '<p>Input geometry unavailable.</p>'
        tags=[r.get(k,'') for k in ['domain','pair','method','status','observation_type']]
        cards.append('<article data-tags="'+esc(' | '.join(map(str,tags)),quote=True)+'"><h3>'+esc(str(r['pair'])+' — '+str(r['method']))+'</h3><p>'+esc(' | '.join(map(str,tags)))+'</p>'+image+'<p>'+' · '.join(links)+'</p><small>'+esc(r['attempt_path'])+'</small></article>')
        completed+=1
        if time.monotonic()-last>=5 or completed==len(rows):
            atomic_progress(dict(phase='figure_generation',finished_attempt_views=completed,planned_attempt_views=len(rows),canonical_figures=len(canonical)));print(f'FIGURES {completed}/{len(rows)} attempt views; {len(canonical)} canonical geometries',flush=True);last=time.monotonic()
    for r in rows:
        r['metric_availability']=json.dumps(r['metric_availability']);r.pop('bounds',None)
    pd.DataFrame(rows).to_csv(target/'per_attempt_metrics.csv',index=False)
    available=[]
    for key in ['selected_indices_A','selected_indices_B','coupling_distance','continuous_A_lower','continuous_B_upper','auxiliary_A','auxiliary_B','extended_vertices_A','dense_states_A','active_states_A','edges_including_waits_A','construction_completed_A','states','states_dequeued','processed_states_lower_bound','queue_size','transitions','unique_configuration_edges','dynamic_programming_labels','warmup_seconds','solver_seconds','validation_seconds','optional_measurement_seconds','total_process_seconds','peak_rss_mib']:
        graph_metric=key in ['states','states_dequeued','processed_states_lower_bound','queue_size','transitions','unique_configuration_edges','dynamic_programming_labels'] or key.startswith(('auxiliary_','extended_vertices_','dense_states_','active_states_','edges_including_waits_','construction_completed_'))
        applicable=[r for r in rows if r['coupling_required']] if graph_metric else rows
        available.append(dict(metric=key,available=sum(r.get(key) is not None for r in applicable),applicable_attempts=len(applicable),missing=sum(r.get(key) is None for r in applicable),not_applicable=len(rows)-len(applicable)))
    pd.DataFrame(available).to_csv(target/'metric_availability.csv',index=False)
    controls=''
    for field in ['domain','pair','method','status','observation_type']:
        options=sorted({str(r.get(field,'')) for r in rows});controls+=f'<label>{field} <select><option value="">All</option>'+''.join('<option>'+html.escape(v)+'</option>' for v in options)+'</select></label> '
    (target/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>CPS saved-attempt gallery</title><style>body{font:15px Arial;margin:2em;background:white;color:#111}header{position:sticky;top:0;background:white;padding:1em;border-bottom:1px solid #aaa}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:20px}article{border:1px solid #ddd;padding:12px;overflow-wrap:anywhere}img{width:100%}select{max-width:220px}small{color:#555}</style><h1>Saved CPS attempts</h1><p>Historical diagnostics, reused evidence, new single observations and timing repetitions remain separate. A figure may be shared when geometry and status agree. No optimization runs in this gallery. Blank metrics mean unavailable; zero is a recorded value; applicability is explicit in JSON. Independent coupling violations are not baseline failures.</p><p><a href="per_attempt_metrics.csv">All metrics CSV</a> · <a href="per_attempt_bounds.csv">Finite bounds CSV</a></p><header>'+controls+' <b id="count"></b></header><main>'+''.join(cards)+'</main><script>function filter(){let values=[...document.querySelectorAll("select")].map(s=>s.value),n=0;document.querySelectorAll("article").forEach(a=>{let tags=a.dataset.tags.split(" | ");a.hidden=!values.every((v,i)=>!v||tags[i]===v);if(!a.hidden)n++});document.getElementById("count").textContent=n+" attempts shown"}document.querySelectorAll("select").forEach(s=>s.onchange=filter);filter()</script>',encoding='utf-8')
    dump(target/'summary.json',dict(attempts=len(rows),validated=sum(r['validated'] for r in rows),canonical_figures=len(canonical),figure_links=sum(r['figure_png'] is not None for r in rows),excluded_or_recovery_records=len(excluded),optimization_executed=False,source_code_sha256=sha(Path(__file__))))
    atomic_progress(dict(phase='completed',finished_attempt_views=len(rows),planned_attempt_views=len(rows),canonical_figures=len(canonical)))
    (target/'README.md').write_text('# Saved-attempt evidence\n\nOne row per phase-separated attempt, deduplicated by attempt directory. Reuse locations are in usage; validation recoveries link to the original raw checkpoint and retain original validation status. Regression injections are excluded and listed. Older pre-checkpoint formats remain in historical result tables and are not silently asserted to have external validation.\n\nGeometry is generated after measured execution. All plotted coordinates come directly from saved NPZ inputs; no smoothing, alignment, auxiliary generation or optimization occurs. Auxiliary locations are matching positions, not selectable output vertices.\n\nGraph edge counts include waits and post-pruning active entries where recorded. Dequeued heap states and processed-state lower bounds remain distinct. Unique configuration edges and DP-label counts are unavailable. Historical missing counters remain blank. Bound occupancy is not execution progress. Fidelity intervals are optional saved continuous-distance brackets; required feasibility predicates are recorded separately. Exact numerical tolerance provenance resides in the hashed implementation, not inferred from these brackets.\n',encoding='utf-8')
    print('Saved gallery:',target,flush=True)

def final20_report():
    """New stratified evidence, separately reported; saved artifacts only."""
    import sys,nbformat,platform
    sys.path[:0]=[str(ROOT/'notebooks_paper'),str(ROOT/'notebooks_paper/graph_only_v3')]
    import coverage_extension as cx
    import expanded_comparison as shared
    global OUT
    chapter=ROOT/'output/paper_v1/experimental_chapter';out=ROOT/'output/paper_v1/hurricane_final_20'
    for sub in ['figures','qa','tables','provenance']:(out/sub).mkdir(exist_ok=True)
    prior_out=cx.OUT;cx.OUT=out
    try:f=cx.save()
    finally:cx.OUT=prior_out
    inputs=pd.read_csv(out/'input_manifest.csv');m=read(out/'configuration_manifest.json');pool=read(out/'selection_pool.json')
    f=f.merge(inputs[['pair','stratum','size_quartile','ratio_group','selection_order']],on='pair',validate='many_to_one')
    f['cohort']='stratified final20';f['observed']=f.origin.eq('new');f.to_csv(out/'results_with_strata.csv',index=False)
    order=inputs.sort_values('selection_order').pair.tolist();aliases={p:f'S{i+1}' for i,p in enumerate(order)}
    f['pair_alias']=f.pair.map(aliases)
    f.to_csv(out/'results_with_strata.csv',index=False)
    alias_table=inputs[['pair','nA','nB','stratum','selection_order','input_sha256','delta1','delta2','delta3']].copy();alias_table.insert(0,'pair_alias',alias_table.pair.map(aliases));alias_table.to_csv(out/'tables/pair_aliases.csv',index=False)
    completion=[]
    for method,g in f.groupby('method',sort=False):completion.append(dict(method=method,planned=len(g),observed=int(g.observed.sum()),validated=int(g.validated.sum()),resource_limited=int(g.resource_limited.sum()),validation_incomplete=int(g.status.eq('validation_incomplete').sum()),pending=int(g.origin.eq('unexecuted').sum())))
    completion=pd.DataFrame(completion);completion.to_csv(out/'tables/completion_by_method.csv',index=False)
    strata=f.groupby(['domain','stratum','method','status'],dropna=False).size().rename('count').reset_index();strata.to_csv(out/'tables/completion_by_stratum.csv',index=False)
    f.groupby(['domain','method','status']).size().rename('count').reset_index().to_csv(out/'tables/termination_by_domain_method.csv',index=False)
    pairing=[]
    for pair,g in f.groupby('pair',sort=False):
        g=g.set_index('method')
        for fidelity,base,joint in [('continuous',METHODS[0],METHODS[2]),('discrete',METHODS[1],METHODS[3])]:
            a,b=g.loc[base],g.loc[joint];both=bool(a.validated and b.validated)
            pairing.append(dict(pair=pair,pair_alias=aliases[pair],fidelity=fidelity,baseline_k=a.k if a.validated else None,joint_k=b.k if b.validated else None,delta_k=b.k-a.k if both else None,baseline_coupling_ratio=a.coupling_distance/a.delta3 if a.validated else None,joint_coupling_ratio=b.coupling_distance/b.delta3 if b.validated else None,equal_k_coupling_repair=bool(both and a.k==b.k and not a.coupling_feasible and b.coupling_feasible),baseline_source=a.attempt_path,joint_source=b.attempt_path))
    paired=pd.DataFrame(pairing);paired.to_csv(out/'tables/matched_quality.csv',index=False)
    old_out=shared.OUT;shared.OUT=out
    try:bounds=shared.graph_tables(f)
    finally:shared.OUT=old_out
    bounds['pair_alias']=bounds.pair.map(aliases);bounds.to_csv(out/'graph_bounds.csv',index=False)
    audit_coverage(f,f,out)
    plot_out=OUT;OUT=out;CAPTIONS.clear()
    plt.rcParams.update({'font.size':9,'axes.facecolor':'white','figure.facecolor':'white','text.color':'black','pdf.fonttype':42})
    statuses=['validated']+sorted(set(f.status)-{'validated'});palette={'validated':'#4daf4a','unexecuted':'.9','transition_budget_limit':'#e69f00','state_budget_limit':'#984ea3','solver_timeout':'#d55e00','validation_incomplete':'#cc79a7','single_chain_edge_budget_limit':'#377eb8','memory_limit':'#a65628'}
    fig,axes=plt.subplots(1,2,figsize=(13,4.3),layout='constrained')
    bottom=np.zeros(4)
    for status in statuses:
        values=np.array([sum((f.method==method)&(f.status==status)) for method in METHODS]);axes[0].bar(range(4),values,bottom=bottom,label=status.replace('_',' '),color=palette.get(status,'.4'));bottom+=values
    axes[0].set_xticks(range(4),[SHORT[x] for x in METHODS]);axes[0].set_ylabel('Configurations / 20 per method');axes[0].legend(fontsize=7)
    slist=list(pool['allocation']);x=np.arange(len(slist));bottom=np.zeros(len(slist))
    for status in statuses:
        values=np.array([sum((f.stratum==s)&(f.status==status)) for s in slist]);axes[1].bar(x,values,bottom=bottom,color=palette.get(status,'.4'));bottom+=values
    axes[1].set_xticks(x,[s.replace('_','\n') for s in slist],fontsize=7);axes[1].set_ylabel('Configurations in stratum');axes[1].set_title('Empty strata retained; denominator = 4 × allocated pairs')
    savefig(fig,'final20_completion','All 80 scheduled configurations, including pending and resource outcomes. Left: 20 per method. Right: recorded allocation by length quartile and ratio; balanced means ratio <=1.5. Empty strata are shown. Completion is not validation success; resource exits are not infeasibility.')
    fig,ax=plt.subplots(figsize=(13,4.5),layout='constrained');positions={p:i for i,p in enumerate(order)}
    for i,method in enumerate(METHODS):
        g=f[f.method.eq(method)&f.validated];xx=g.pair.map(positions)+(i-1.5)*.16;ax.scatter(xx,g.k,color=COLORS[method],label=method)
        if method.startswith('Independent'):
            bad=g.coupling_feasible.eq(False);ax.scatter(xx[bad],g.loc[bad,'k'],marker='x',s=90,color='black')
    ax.set_xticks(range(len(order)),[aliases[p] for p in order]);ax.set_ylabel(r'$k=\max(k_A,k_B)$');ax.legend(ncol=4,fontsize=8)
    savefig(fig,'final20_quality','Matched input pairs and thresholds, only externally validated outputs. Black crosses mark independent coupling violations, not independent-method failures. Missing solutions are omitted, never assigned zero. S aliases follow frozen selection order.')
    coupling_figure(paired,aliases,order);bound_figure(bounds,aliases,order)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for method in METHODS:
        g=f[f.method.eq(method)&f.validated];axes[0].scatter(g.retained_total,g.solver_seconds,color=COLORS[method],label=method)
        failed=f[f.method.eq(method)&f.resource_limited];axes[1].scatter(failed.nA+failed.nB,failed.solver_seconds,color=COLORS[method],marker='x',s=70)
    axes[0].set_xlabel(r'Retained fraction $(k_A+k_B)/(n_A+n_B)$');axes[0].set_ylabel('Successful solver time (s)');axes[0].set_yscale('log');axes[0].legend(fontsize=8)
    axes[1].set_xlabel('Input vertices nA+nB');axes[1].set_ylabel('Resource time to termination (s)');axes[1].set_title('Output quality unknown for resource exits')
    if not f.resource_limited.any():axes[1].text(.5,.5,'No resource exits recorded',transform=axes[1].transAxes,ha='center')
    savefig(fig,'final20_compression_cost','Continuous and discrete fidelity define different feasible sets. Left: retained fraction versus solver time for validated solutions only. Right: resource-exit times, with no invented output quality. These are single attempts, not repeated timing estimates.')
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for method in METHODS:
        g=f[f.method.eq(method)]
        for ax,key in zip(axes,['input_size','states','transitions']):
            for valid,marker in [(True,'o'),(False,'x')]:
                z=g[g.validated.eq(valid)&g.solver_seconds.notna()];xx=z.nA+z.nB if key=='input_size' else z[key];keep=xx.notna();ax.scatter(xx[keep],z.loc[keep,'solver_seconds'],c=COLORS[method],marker=marker,s=30)
            ax.set_xlabel({'input_size':'Input vertices nA+nB','states':'Discovered states','transitions':'Candidate transition checks'}[key]);ax.set_yscale('log');ax.set_ylabel('Solver seconds')
        
    method_legend(axes[0]);savefig(fig,'final20_work','Descriptive work indicators only: input size, discovered states, and candidate transition checks. Circles are validated solution times; crosses are unsuccessful termination times. Missing counters are omitted. Transition checks are not unique edges or all operations. No complexity law is fitted.')
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    for method in METHODS:
        g=f[f.method.eq(method)&f.peak_rss_mib.notna()];ax.scatter(g.nA+g.nB,g.peak_rss_mib,color=COLORS[method],label=method)
        bad=g.resource_limited;ax.scatter((g.nA+g.nB)[bad],g.peak_rss_mib[bad],marker='x',color='black',s=80)
    ax.set_xlabel('Input vertices nA+nB');ax.set_ylabel('Sampled process-tree peak RSS (MiB)');ax.legend(fontsize=8)
    savefig(fig,'final20_memory','Supplementary sampled process-tree RSS, including runtime overhead; not exact algorithm memory. Black crosses indicate resource exits. Missing measurements are not zero.')
    captions=CAPTIONS.copy();(out/'figure_captions.md').write_text('\n\n'.join('## '+k+'\n'+v for k,v in captions.items()),encoding='utf-8')
    render_qa();OUT=plot_out
    numbers=dict(planned=80,observed=int(f.observed.sum()),validated=int(f.validated.sum()),resource_limited=int(f.resource_limited.sum()),validation_incomplete=int(f.status.eq('validation_incomplete').sum()),pending=int(f.origin.eq('unexecuted').sum()),eligible_pool=pool['eligible_pairs'],selected_pairs=len(inputs),reused_prior_storms=sum(len(r['reused_prior_storms']) for r in read(out/'input_manifest.json')['inputs']),reused_within_storms=sum(len(r['reused_current_storms']) for r in read(out/'input_manifest.json')['inputs']))
    for fidelity in ['continuous','discrete']:
        p=paired[paired.fidelity.eq(fidelity)];numbers[fidelity+'_matched']=int(p.delta_k.notna().sum());numbers[fidelity+'_equal']=int(p.delta_k.eq(0).sum());numbers[fidelity+'_extra']=int(p.delta_k.gt(0).sum());numbers[fidelity+'_equal_k_repairs']=int(p.equal_k_coupling_repair.sum())
    dump(out/'summary.json',numbers)
    earlier=pd.read_csv(ROOT/'output/paper_v1/coverage_extension/matched_comparisons.csv');earlier=earlier[~earlier.pair.isin(order)]
    layers=[]
    for layer,data in [('earlier pilot/development cohorts',earlier),('stratified final20',paired)]:
        for fidelity,g in data.groupby('fidelity'):
            layers.append(dict(layer=layer,fidelity=fidelity,matched_validated=int(g.delta_k.notna().sum()),equal_objective=int(g.delta_k.eq(0).sum()),higher_coupled_objective=int(g.delta_k.gt(0).sum()),equal_k_coupling_repairs=int(g.equal_k_coupling_repair.sum())))
    pd.DataFrame(layers).to_csv(out/'tables/quality_by_evidence_layer.csv',index=False)
    pd.DataFrame([dict(claim=k,value=v,source='tables/matched_quality.csv' if k.startswith(('continuous_','discrete_')) else 'selection_pool.json and input_manifest.json' if k in ['eligible_pool','selected_pairs','reused_prior_storms','reused_within_storms'] else 'results_with_strata.csv') for k,v in numbers.items()]).to_csv(out/'claim_traceability.csv',index=False)
    latex_table(completion,out/'tables/completion.tex','Stratified final20 cohort. Twenty planned configurations per method; pending slots are not observations.','tab:final20-completion')
    quality=pd.DataFrame([dict(Fidelity=q,Matched=numbers[q+'_matched'],Equal=numbers[q+'_equal'],Higher=numbers[q+'_extra'],EqualKRepairs=numbers[q+'_equal_k_repairs']) for q in ['continuous','discrete']]);latex_table(quality,out/'tables/quality.tex','Within-fidelity matched objectives; independent coupling violations are permitted.','tab:final20-quality')
    progress=read(out/'progress.json') if (out/'progress.json').exists() else {};pause=progress.get('pause_reason','')
    report=f'''# Stratified final hurricane cohort

## Scope and selection
Seed 20260916; {pool['eligible_pairs']} eligible supplied unordered pairs; full track sizes {pool['track_size_range']}. Quartile cutpoints: {pool['quartile_cutpoints']}; ratio split at 1.5. Allocation: {pool['allocation']}. All candidates, exclusions, empty strata and deterministic redistribution are saved. No upper-size cap. The selected tracks span {int(inputs[['nA','nB']].min().min())}–{int(inputs[['nA','nB']].max().max())} vertices. Prior storm reuses: {numbers['reused_prior_storms']}; within-cohort reuses: {numbers['reused_within_storms']}.

Preference-constrained stratified random sampling is not simple random sampling. Unweighted percentages describe this cohort, not the hurricane archive. Alpha=.5 was chosen using earlier pilot evidence. The sampling frame is the supplied pairs.csv, not every combinatorial pairing of all storms.

## Execution
{numbers['observed']}/80 observed; {numbers['validated']} validated; {numbers['resource_limited']} resource-limited; {numbers['validation_incomplete']} validation-incomplete; {numbers['pending']} pending. Current phase: {progress.get('phase')}. Pause reason: {pause or 'none'}.
Resource exits are completed observations and never retried. Raw checkpoints precede external validation. Independent coupling is descriptive. Validation checks feasibility, not independent optimality. Common external 60-second/512-MiB budgets coexist with disclosed method-specific internal work limits. No optimization algorithm, endpoint rule, tolerance, warmup, cache policy or thread setting was changed.

## Matched quality
Continuous: {numbers['continuous_equal']}/{numbers['continuous_matched']} matched objectives equal; {numbers['continuous_extra']} higher coupled objectives; {numbers['continuous_equal_k_repairs']} equal-k coupling repairs.
Discrete: {numbers['discrete_equal']}/{numbers['discrete_matched']} matched objectives equal; {numbers['discrete_extra']} higher coupled objectives.
Continuous/discrete comparisons have different fidelity feasible sets. Unvalidated quality stays missing. Timing and memory are single-attempt observations; no scaling law or significance claim is supported.

## Evidence layers and files
This new cohort is reported separately in results_with_strata.csv, tables/, graph_bounds.csv, and figures/. Earlier cohorts remain in ../coverage_extension/combined_results.csv with explicit cohort labels. The separate 120-repetition study is unchanged. Original geometric illustrations are retained; the gallery covers every saved attempt rather than selecting only attractive outcomes.
[Gallery](../coverage_extension/attempt_gallery/index.html), [metrics](../coverage_extension/attempt_gallery/per_attempt_metrics.csv). Auxiliary locations are matching positions, not selectable output vertices. Transition checks are not unique configuration edges; unique edges and DP labels remain unavailable. Finite bounds count like objects; interrupted counts are observed before termination, not progress percentages.

## Resume / stopping point
If pending slots remain: `.venv/Scripts/python.exe notebooks_paper/coverage_extension.py execute --cohort hurricane_final_20` after readiness passes. Recorded attempts are skipped. No further datasets, thresholds, budget increases or timing rounds are authorized.
'''
    (out/'report.md').write_text(report,encoding='utf-8');(out/'status.md').write_text(report,encoding='utf-8')
    chapter_path=chapter/'experiments.tex';body=chapter_path.read_text(encoding='utf-8');body=re.sub(r'\\section\{Final stratified hurricane cohort\}.*?(?=\\section\{Limitations\})','',body,flags=re.S)
    section=r'\section{Final stratified hurricane cohort}\label{sec:final20}'+'\n'
    section+=f"A separately frozen cohort uses seed 20260916 to select 20 pairs from {pool['eligible_pairs']} eligible supplied pairs. All source vertices and the common basin frame are retained. The eligible track range is {pool['track_size_range'][0]}--{pool['track_size_range'][1]} vertices, without an upper-size cap. Maximum-length quartile cutpoints are {', '.join(map(str,pool['quartile_cutpoints']))}; each quartile is subdivided by length ratio at 1.5. Empty strata receive zero slots; round-robin allocation over nonempty strata redistributes shortages deterministically. Within each stratum, seeded random choice among candidates minimizing prior or current storm reuse gives a preference-constrained stratified sample. It is not a simple random sample, and unweighted cohort percentages are not archive-wide estimates. Alpha was chosen using pilot evidence.\n\n"
    section+=f"The frozen rotated method schedule has {numbers['observed']} of 80 recorded observations, {numbers['validated']} validated solutions, {numbers['resource_limited']} resource exits and {numbers['pending']} pending slots. Resource termination does not establish infeasibility. No outcome-based replacements or retries occur. These single observations remain separate from the original 120 timing repetitions.\n\n"
    section+=r'\input{../hurricane_final_20/tables/completion.tex}'+'\n'+r'\input{../hurricane_final_20/tables/quality.tex}'+'\n'
    section+=f"Under continuous fidelity, equal independent and coupled objectives occur in {numbers['continuous_equal']} of {numbers['continuous_matched']} matched validated cases; {numbers['continuous_extra']} coupled objectives are higher, and {numbers['continuous_equal_k_repairs']} cases restore coupling at equal objective. Under discrete fidelity, {numbers['discrete_equal']} of {numbers['discrete_matched']} matched objectives are equal. These within-fidelity comparisons isolate coupling; CPS-2F versus CPS-3F instead compares different fidelity feasible sets.\n\n"
    for name,caption in captions.items():
        section+=r'Figure~\ref{fig:final20-'+name+'} presents the saved evidence for this cohort.\n'
        section+='\n'.join([r'\begin{figure}[p]',r'\centering',r'\includegraphics[width=\textwidth]{../hurricane_final_20/figures/'+name+'.pdf}',r'\caption{'+tex(caption)+'}',r'\label{fig:final20-'+name+'}',r'\end{figure}'])+'\n'
    body=body.replace(r'\section{Limitations}',section+r'\section{Limitations}',1);chapter_path.write_text(body,encoding='utf-8')
    notebook=ROOT/'notebooks_paper/expanded_comparison_saved_results.ipynb';nb=nbformat.read(notebook,4);nb.cells=[c for c in nb.cells if not c.metadata.get('final20')]
    def add(kind,text):
        c=nbformat.v4.new_markdown_cell(text) if kind=='md' else nbformat.v4.new_code_cell(text);c.metadata['final20']=True;nb.cells.append(c)
    add('md','# Separately sampled final20 cohort\nSeeded stratification; single attempts; earlier pilot and timing evidence remain distinct.')
    add('code',"FINAL20 = PROJECT/'output/paper_v1/hurricane_final_20'\ndisplay(Markdown((FINAL20/'report.md').read_text(encoding='utf-8')))\ndisplay(pd.read_csv(FINAL20/'tables/completion_by_method.csv'))")
    for name,caption in captions.items():add('md',caption);add('code',f"display(Image(filename=str(FINAL20/'figures/{name}.png')))")
    nbformat.write(nb,notebook)
    OUT=chapter;structural_check();OUT=plot_out
    dump(out/'provenance/report_build.json',dict(utc=datetime.now(timezone.utc).isoformat(),reporting_python=sys.executable,reporting_python_version=sys.version,platform=platform.platform(),source_sha256=sha(Path(__file__)),source_results_sha256=sha(out/'results.csv'),optimization_executed=False))
    print('Final20 saved report',numbers,flush=True)

def final20_audit():
    import sys,platform,psutil,numba,os
    out=ROOT/'output/paper_v1/hurricane_final_20';m=read(out/'configuration_manifest.json');pool=read(out/'selection_pool.json');inputs=read(out/'input_manifest.json')['inputs']
    workers=[p.info for p in psutil.process_iter(['pid','cmdline']) if any('worker.py' in str(a) for a in (p.info['cmdline'] or []))];assert not workers
    assert len(m['rows'])==80 and len(inputs)==20 and len({tuple(sorted(r['individuals'])) for r in inputs})==20
    prior={tuple(p) for p in pool['prior_pairs']};assert not any(tuple(sorted(r['individuals'])) in prior for r in inputs)
    used={z for p in prior for z in p};strata={s:[p for p in pool['pool'] if p['stratum']==s] for s in pool['allocation']};taken={s:0 for s in strata};rng=np.random.default_rng(20260916);chosen=[]
    while len(chosen)<20:
        for s in strata:
            if taken[s]>=pool['allocation'][s]:continue
            options=[p for p in strata[s] if p['pair'] not in chosen];cost=min(sum(z in used for z in p['individuals']) for p in options);tied=[p for p in options if sum(z in used for z in p['individuals'])==cost];pick=tied[int(rng.integers(len(tied)))];chosen.append(pick['pair']);used.update(pick['individuals']);taken[s]+=1
    assert chosen==[p['pair'] for p in inputs]
    template=read(ROOT/'output/paper_v1/hurricane_extension_3/configuration_manifest.json')['rows'][0]['configuration'];records=[]
    for i,p in enumerate(inputs):
        arrays=np.load(ROOT/p['input_file']);assert sha(ROOT/p['input_file'])==p['input_sha256']
        for side,identity in zip(['A','B'],p['individuals']):assert np.array_equal(arrays[side],np.loadtxt(ROOT/'data/hurricanes'/f'{identity}.txt'))
        assert p['delta1']==.5*np.linalg.norm(np.diff(arrays['A'],axis=0),axis=1).mean() and p['delta2']==.5*np.linalg.norm(np.diff(arrays['B'],axis=0),axis=1).mean()
        methods=[r['method'] for r in m['rows'][i*4:i*4+4]];shift=i%4;assert methods==METHODS[shift:]+METHODS[:shift]
    for r in m['rows']:
        c=r['configuration'];assert c['certificate'] is False and c['alpha']==.5
        for key in ['graph_limits','phase_budgets','memory_mib','algorithms']:assert c[key]==template[key]
        for p,h in c['algorithms'].items():assert sha(ROOT/p)==h
        attempts=list((out/'slots'/f"{r['plan_order']:02d}"/'attempts').glob('*/configuration.json'));assert len(attempts)<=1
        row=dict(configuration_id=r['configuration_id'],attempts=len(attempts),input_sha256=c['input_sha256'],configuration_sha256=sha(out/'slots'/f"{r['plan_order']:02d}"/'configuration.json'),status='pending')
        if attempts:
            d=attempts[0].parent;raw=read(d/'raw_solver.json') if (d/'raw_solver.json').exists() else {};v=read(d/'validation.json') if (d/'validation.json').exists() else {}
            row['status']=v.get('status','resource exit / validation unavailable');row['raw_sha256']=sha(d/'raw_solver.json') if raw else None
            if v.get('status')=='validated':
                assert v['raw_checkpoint_sha256']==row['raw_sha256'] and raw['unix']<=v['unix']
                if c['method'].startswith('CPS'):assert raw['certificate_enabled'] is False and raw['solver_status']=='optimal_configuration_graph' and raw['observation']['configuration_graph_entered']
        records.append(row)
    pd.DataFrame(records).to_csv(out/'provenance/execution_audit.csv',index=False)
    gallery=ROOT/'output/paper_v1/coverage_extension/attempt_gallery';metrics=pd.read_csv(gallery/'per_attempt_metrics.csv');own=metrics[metrics.attempt_path.str.contains('hurricane_final_20')];own.to_csv(out/'per_attempt_metrics.csv',index=False)
    bounds=pd.read_csv(gallery/'per_attempt_bounds.csv');bounds[bounds.attempt_path.str.contains('hurricane_final_20')].to_csv(out/'per_attempt_bounds.csv',index=False)
    for r in own.itertuples():assert (ROOT/r.figure_pdf).exists() and (ROOT/r.figure_png).exists()
    timestamp=datetime.now(timezone.utc).isoformat()
    dump(out/'provenance/environment_audit.json',dict(utc=timestamp,interpreter=sys.executable,interpreter_sha256=sha(Path(sys.executable)),python=sys.version,numpy=np.__version__,numba=numba.__version__,platform=platform.platform(),memory_gib=psutil.virtual_memory().total/2**30,logical_cpus=psutil.cpu_count(),scope='Audit-time environment; unchanged supervisor launches sys.executable. Per-attempt code/configuration hashes remain authoritative.',worker_source_sha256=sha(ROOT/'notebooks_paper/graph_only_v3/worker.py'),supervisor_source_sha256=sha(ROOT/'notebooks_paper/graph_only_v3/supervisor.py'),thread_and_cache_policy='Existing supervisor/common policy; output-local caches, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1; no settings changed'))
    result=dict(utc=timestamp,active_workers=workers,selection_reproduced=True,rotated_order_verified=True,exact_tracks_verified=True,unchanged_resource_policy=True,algorithm_hashes_verified=True,scheduled=80,attempted=sum(r['attempts'] for r in records),remaining=sum(r['attempts']==0 for r in records),gallery_links=len(own),optimization_executed=False)
    dump(out/'provenance/final_audit.json',result)
    dump(out/'artifact_hashes.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.name!='artifact_hashes.json' and 'cache' not in p.parts})
    print(result,flush=True)

def verify_gallery():
    import pymupdf,psutil
    target=ROOT/'output/paper_v1/coverage_extension/attempt_gallery';qa=target/'qa';qa.mkdir(exist_ok=True)
    workers=[p.info for p in psutil.process_iter(['pid','cmdline']) if any('worker.py' in str(a) for a in (p.info['cmdline'] or []))]
    assert not workers,workers
    table=pd.read_csv(target/'per_attempt_metrics.csv');combined=pd.read_csv(target.parent/'combined_results.csv')
    assert table.attempt_id.is_unique and combined.configuration_id.is_unique
    for row in table.itertuples():
        assert (ROOT/row.figure_png).exists() and (ROOT/row.figure_pdf).exists()
        details=read(target/'metrics'/f'{row.attempt_id}.json');assert details['metrics']['attempt_path']==row.attempt_path
    for row in combined[combined.observed].itertuples():
        assert (ROOT/row.figure_png).exists() and (ROOT/row.figure_pdf).exists() and (ROOT/row.attempt_metrics_json).exists()
    new=ROOT/'output/paper_v1/hurricane_extension_3';m=read(new/'configuration_manifest.json')
    prior=read(target.parent/'configuration_manifest.json');template=next(r['configuration'] for r in prior['rows'] if r['eligible'])
    for row in m['rows']:
        c=row['configuration'];assert c['certificate'] is False and c['alpha']==.5
        for key in ['graph_limits','phase_budgets','memory_mib','algorithms']:assert c[key]==template[key]
        assert sha(ROOT/c['input_file'])==c['input_sha256']
        for p,h in c['algorithms'].items():assert sha(ROOT/p)==h
        arrays=np.load(ROOT/c['input_file'])
        for side,identity in zip(['A','B'],c['pair'].split('__')):assert np.array_equal(arrays[side],np.loadtxt(ROOT/'data/hurricanes'/f'{identity}.txt'))
        assert c['delta1']==.5*np.linalg.norm(np.diff(arrays['A'],axis=0),axis=1).mean()
        assert c['delta2']==.5*np.linalg.norm(np.diff(arrays['B'],axis=0),axis=1).mean()
        assert len(list((new/'slots'/f"{row['plan_order']:02d}"/'attempts').glob('*/configuration.json')))==1
    bounds=pd.read_csv(target/'per_attempt_bounds.csv');comparable=bounds.dropna(subset=['observed','finite_upper_bound'])
    assert (comparable.observed<=comparable.finite_upper_bound+1e-6).all()
    for p in sorted(set(table.figure_pdf)):
        with pymupdf.open(ROOT/p) as doc:assert len(doc)==1 and doc[0].rect.width>0
    chosen=[table[(table.domain=='protein') & table.validated].iloc[0],table[(table.domain=='hurricane') & table.validated & table.coupling_feasible.eq(False)].iloc[0],table[table.status.eq('transition_budget_limit')].iloc[0],table[table.pair.eq('EP061986__EP131991') & table.method.eq('CPS-2F')].iloc[0]]
    for i,row in enumerate(chosen):
        with pymupdf.open(ROOT/row.figure_pdf) as doc:doc[0].get_pixmap(matrix=pymupdf.Matrix(1,1)).save(qa/f'example_{i}.png')
    audit=dict(active_workers=workers,attempt_records=len(table),canonical_figures=table.figure_pdf.nunique(),primary_observations=int(combined.observed.sum()),primary_validated=int(combined.validated.sum()),primary_figure_links=int(combined[combined.observed].figure_png.notna().sum()),new_attempts=12,one_attempt_each=True,source_tracks_exact=True,resource_policy_unchanged=True,algorithm_hashes_unchanged=True,comparable_bounds_checked=len(comparable),pdfs_opened=table.figure_pdf.nunique(),optimization_executed=False)
    dump(target/'delivery_audit.json',audit)
    latest=pd.read_csv(new/'results.csv');matched=pd.read_csv(target.parent/'matched_comparisons.csv')
    lines=['# Completed continuation and three-pair extension','',
        'The previously pending 27 eligible coverage-extension configurations completed once each. The coverage extension now has 48/48 validated eligible observations; eight preparation-blocked slots remain unexecuted. Historical resource exits were not retried.',
        '',f'The three-pair extension has {len(latest)} recorded single observations, {int(latest.validated.sum())} validated, {int(latest.resource_limited.sum())} resource exits. No timing repetitions were added.',
        '', '| Pair | nA/nB | Independent continuous k | CPS-2F k | Independent discrete k | CPS-3F k |','|---|---:|---:|---:|---:|---:|']
    for pair,g in latest.groupby('pair',sort=False):
        methods=g.set_index('method');first=g.iloc[0]
        values=[num(methods.loc[method,'k']) for method in [METHODS[0],METHODS[2],METHODS[1],METHODS[3]]]
        lines.append('| '+pair+f" | {int(first.nA)}/{int(first.nB)} | "+' | '.join(values)+' |')
    lines+=['',f'Combined primary coverage: {int(combined.observed.sum())}/{len(combined)} intended configurations observed, {int(combined.validated.sum())} validated, {int(combined.resource_limited.sum())} resource exits, {int((~combined.eligible).sum())} preparation-blocked. Repetitions and historical diagnostics are excluded from this denominator.',
        '', '## Scientific comparison']
    for fidelity in ['continuous','discrete']:
        g=matched[matched.fidelity.eq(fidelity)];lines.append(f'- {fidelity.capitalize()} fidelity: {int(g.delta_k.eq(0).sum())}/{int(g.delta_k.notna().sum())} validated matched objectives are equal; {int(g.delta_k.gt(0).sum())} coupled objectives are larger.')
    lines+=['- Independent coupling violations are permitted. Continuous/discrete methods have different feasible sets. Validation establishes feasibility, not independent proof of optimality.',
        '- Shared protein structure, subsampling, differing historical preparation regimes, pilot-informed alpha, and single new observations limit generalization. Resource exits do not establish infeasibility.',
        '', '## Evidence and presentation',f'- [Gallery](../coverage_extension/attempt_gallery/index.html): {len(table)} attempt records, {table.figure_pdf.nunique()} canonical PNG/PDF figures. Includes legacy integrated validation, phase-separated validation and timing repetitions with distinct provenance; these counts are not unique-configuration coverage.',
        '- [Metrics](../coverage_extension/attempt_gallery/per_attempt_metrics.csv), [availability](../coverage_extension/attempt_gallery/metric_availability.csv), and [finite bounds](../coverage_extension/attempt_gallery/per_attempt_bounds.csv). Exact thresholds, selected indices, hashes, recovery links, graph counters and phase/resource records are retained where available.',
        '- Unique configuration-graph edge and DP-label counts are unavailable. Missing historical measurements remain missing; zeros and not-applicable fields are distinguished. Transition checks are not unique edges. Bound occupancy is not search progress.',
        '- Figures were generated after timed execution. Protein views are static 3D; hurricanes retain the common projected frame. Unvalidated attempts show inputs only. Identical timing-repeat geometries may share a figure while their measurements remain separate.',
        '- [Combined results](../coverage_extension/combined_results.csv), [chapter](../experimental_chapter/experiments.tex), and existing saved-results notebook were updated. The LaTeX draft passed structural checks; no compiler was available.',
        '', '## Implementation and stopping point',
        'Only orchestration output/schedule parameterization and saved-evidence reporting were extended. Worker, validator and optimization algorithms are unchanged. The pre-extension orchestrator source is retained as provenance text. One attempt per new slot, no budget increases, no certificate shortcuts, no new timing rounds. No experiment worker remains active.']
    (new/'completion_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');print(audit,flush=True)

if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='verify_notebook':verify_notebook()
    elif len(sys.argv)>1 and sys.argv[1]=='coverage_report':coverage_report()
    elif len(sys.argv)>1 and sys.argv[1]=='attempt_gallery':attempt_gallery()
    elif len(sys.argv)>1 and sys.argv[1]=='verify_gallery':verify_gallery()
    elif len(sys.argv)>1 and sys.argv[1]=='final20_report':final20_report()
    elif len(sys.argv)>1 and sys.argv[1]=='final20_audit':final20_audit()
    else:main()
