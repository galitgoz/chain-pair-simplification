"""Read-only CPS-2F evidence extraction and plotting. Never imports optimizers."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/paper_v1/plots_auxiliary_bound_percentages'
SOURCE = ROOT / 'output/paper_v1/coverage_extension/combined_results.csv'
METRICS = ROOT / 'output/paper_v1/coverage_extension/attempt_gallery/per_attempt_metrics.csv'
BLUE, GREEN = '#2675B8', '#388443'

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    return json.loads(p.read_text(encoding='utf-8'))

def collect():
    combined = pd.read_csv(SOURCE)
    primary = combined[combined.method.eq('CPS-2F') & combined.alpha.eq(.5) & combined.observed.eq(True)].copy()
    metrics = pd.read_csv(METRICS)
    metrics['normalized_path'] = metrics.attempt_path.str.replace('\\', '/', regex=False)
    assert primary.pair.is_unique and len(primary) == 45
    records = []
    for r in primary.itertuples():
        attempt = ROOT / r.attempt_path
        config = read(attempt/'configuration.json')
        assert config['certificate'] is False and config['alpha'] == .5
        assert config['input_sha256'] == digest(ROOT/config['input_file'])
        rawpath = attempt/'raw_solver.json'
        progresspath = attempt/'progress.json'
        raw = read(rawpath) if rawpath.exists() else {}
        progress = read(progresspath) if progresspath.exists() else {}
        observation = raw.get('observation', progress)
        assert observation.get('chain_graphs_completed') == 2
        assert observation.get('configuration_graph_entered') is True
        limited = r.status == 'transition_budget_limit'
        assert limited or r.status == 'validated'
        if limited:
            assert read(attempt/'solver_exit.json')['status'] == 'transition_budget_limit'
        matching = metrics[metrics.normalized_path.eq(attempt.relative_to(ROOT).as_posix())]
        assert len(matching) == 1
        metric = matching.iloc[0]
        record = dict(pair=r.pair, alias=r.pair_alias, domain=r.domain, alpha=.5,
                      nA=int(r.nA), nB=int(r.nB), delta1=r.delta1, delta2=r.delta2, delta3=r.delta3,
                      units=r.units, input_sha256=config['input_sha256'],
                      attempt_path=attempt.relative_to(ROOT).as_posix(), status=r.status,
                      subsequent_search_resource_limited=limited,
                      common_comparison_set=bool(combined[combined.pair.eq(r.pair) & combined.alpha.eq(.5)].validated.eq(True).all()),
                      preprocessing='w=16; endpoint-preserving protein subsampling' if r.domain=='protein' else 'Full supplied hurricane tracks; no additional subsampling',
                      configuration_sha256=digest(attempt/'configuration.json'))
        for side in ['A','B']:
            graph = observation.get('single_graph_'+side)
            n = record['n'+side]
            if graph is not None:
                assert graph['construction_completed'] is True
                assert graph['original_vertices'] == n
                auxiliary = int(graph['auxiliary'])
                assert graph['extended_vertices'] - n == auxiliary
                evidence = rawpath if raw else progresspath
                field = ('observation.' if raw else '')+'single_graph_'+side+'.auxiliary'
                completion = 'Explicit completed single-curve graph; auxiliary generation finished before graph return'
            else:
                # Older successful checkpoints recorded counts but no per-graph flags.
                assert not limited and raw['solver_status'] == 'optimal_configuration_graph'
                auxiliary = int(raw['raw_result']['auxiliary_'+side])
                evidence = rawpath
                field = 'raw_result.auxiliary_'+side
                completion = 'Successful configuration-graph checkpoint plus chain_graphs_completed=2; per-curve flag not recorded'
            assert auxiliary == metric['auxiliary_'+side]
            bound = 2*n*(n-1)
            assert 0 <= auxiliary <= bound
            record.update({
                'a'+side: auxiliary, 'bound_'+side: bound,
                'percentage_'+side: 100*auxiliary/bound,
                'exact_fraction_'+side: f'{100*auxiliary}/{bound}',
                'construction_completed_'+side: True,
                'construction_evidence_'+side: completion,
                'count_source_'+side: evidence.relative_to(ROOT).as_posix(),
                'count_field_'+side: field, 'count_source_sha256_'+side: digest(evidence)})
        records.append(record)
    df = pd.DataFrame(records)
    df['alias_number'] = df.alias.str[1:].astype(int)
    df = df.sort_values(['domain','alias_number']).reset_index(drop=True)
    assert df.domain.value_counts().to_dict() == {'hurricane':33,'protein':12}
    assert df[df.common_comparison_set].domain.value_counts().to_dict() == {'hurricane':31,'protein':12}
    assert set(df[df.subsequent_search_resource_limited].alias) == {'H4','H18'}
    return df

def draw(df, name, include_proteins=False, compact=False):
    hurricanes = df[df.domain.eq('hurricane')]
    groups = []
    if include_proteins:
        groups.append(('Proteins | existing w=16 subsampling', df[df.domain.eq('protein')]))
    for i, indices in enumerate(np.array_split(np.arange(len(hurricanes)), 3)):
        groups.append((f'Hurricanes | panel {i+1}/3 | no additional subsampling', hurricanes.iloc[indices]))
    count = len(groups)
    height = 4.7*count+1.
    fig, axes = plt.subplots(count, 1, figsize=(7.6,10.8) if compact else (18,height), squeeze=False)
    fig.subplots_adjust(left=.11 if compact else .07, right=.985, top=.91, bottom=.105 if compact else 3./height,
                        hspace=.75 if compact else 1.15)
    upper = 13. if include_proteins else 5.
    for idx, (title, data) in enumerate(groups):
        ax = axes[idx,0]; x = np.arange(len(data)); width=.36
        for side, offset, color, hatch in [('A',-width/2,BLUE,None),('B',width/2,GREEN,'///')]:
            bars=ax.bar(x+offset,data['percentage_'+side],width,color=color,edgecolor='black',linewidth=.45,hatch=hatch,zorder=3)
            ax.bar_label(bars,labels=[f'{v:.1f}%' for v in data['percentage_'+side]],padding=3,
                         fontsize=7 if compact else 9,rotation=90 if compact else 0)
        labels = [(r.alias if compact else r.pair)+('*' if r.subsequent_search_resource_limited else '') for r in data.itertuples()]
        ax.set_xticks(x,labels,rotation=0 if compact else 90,fontsize=8 if compact else 10)
        ax.set_ylim(0,upper)
        ax.set_yticks([0,3,6,9,12] if include_proteins else np.arange(6))
        ax.set_ylabel('Own auxiliary bound (%)',fontsize=9 if compact else 12)
        ax.set_title(f'({chr(97+idx)}) {title}',loc='left',fontsize=10 if compact else 13,pad=12)
        ax.grid(axis='y',color='#DDDDDD',linewidth=.6,zorder=0)
        ax.spines[['top','right']].set_visible(False)
        ax.margins(x=.025)
    scope = 'Common comparison set: 31 hurricane + 12 protein pairs' if compact else ('All tested pairs: 33 hurricane + 12 protein pairs' if include_proteins else 'All 33 tested hurricane pairs')
    fig.suptitle('CPS-2F added auxiliary locations | '+r'$\alpha=0.5$'+'\n'+scope,fontsize=12 if compact else 18,y=.985)
    fig.legend(handles=[Patch(facecolor=BLUE,edgecolor='black',label='Curve A (solid)'),
                        Patch(facecolor=GREEN,edgecolor='black',hatch='///',label='Curve B (hatched)')],
               loc='upper center',bbox_to_anchor=(.5,.955 if compact else .955),ncol=2,frameon=False,fontsize=9 if compact else 12)
    note = ('Each bar: 100 × added locations / [2n(n−1)], using its own curve size n; original vertices excluded.\n'
            'Proteins: existing w=16; hurricanes: no additional subsampling. ' if include_proteins else
            'Each bar: 100 × added locations / [2n(n−1)], using its own curve size n; original vertices excluded.\nHurricanes: full supplied tracks, no additional subsampling. ')
    note += '\nNormalized percentages do not compare absolute counts or isolate a domain effect.'
    note += ('\nCommon set excludes H4 and H18; stable aliases are mapped in pair_aliases.csv.' if compact else
             '\n* Auxiliary and single-curve construction completed; the subsequent CPS search reached its transition limit.')
    fig.text(.5,.018,note,ha='center',va='bottom',fontsize=7.5 if compact else 11,linespacing=1.5)
    for ext in ['pdf','png']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=300,facecolor='white',metadata={'Title':name} if ext=='pdf' else None)
    plt.close(fig)
    df.drop(columns=['alias_number']).to_csv(OUT/f'{name}_values.csv',index=False,float_format='%.17g')

def build():
    OUT.mkdir(parents=True,exist_ok=True)
    df=collect()
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,'text.color':'black',
                         'axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','figure.facecolor':'white','axes.facecolor':'white'})
    draw(df[df.domain.eq('hurricane')],'hurricanes_all33')
    draw(df,'combined_all45',include_proteins=True)
    draw(df[df.common_comparison_set],'publication_common43',include_proteins=True,compact=True)
    df[['alias','pair','domain','nA','nB','common_comparison_set','subsequent_search_resource_limited','preprocessing']].to_csv(OUT/'pair_aliases.csv',index=False)
    captions='''# English figure captions

## hurricanes_all33
Added CPS-2F auxiliary matching locations as percentages of each curve's finite construction bound, for all 33 tested hurricane pairs at alpha=0.5. Blue solid bars denote curve A and green hatched bars denote curve B. All supplied track vertices are retained, without additional subsampling. The three panels show complete pair identifiers. Asterisks mark AL072006__AL011900 (H4) and AL061921__AL101926 (H18): both curves' auxiliary and single-curve graph construction completed, with recorded counts, before the subsequent CPS configuration search reached its transition-check limit. Their construction counts are complete, although no simplified solution was returned.

## combined_all45
Added CPS-2F auxiliary-location utilization for all 45 tested pairs (33 hurricane and 12 protein pairs), at alpha=0.5. The protein panel uses the existing endpoint-preserving w=16 subsampling; the hurricane panels use full supplied tracks with no additional subsampling. All panels use the same 0–13% y-axis scale. Full pair identifiers and percentage labels are shown. Asterisks have the construction-complete/search-resource-limited meaning defined above.

## publication_common43
Added CPS-2F auxiliary-location utilization on the common validated four-method comparison set: 31 hurricane and 12 protein pairs, at alpha=0.5. The two transition-limited hurricane cases are excluded from this comparison set, not from the all-pair construction figures. Domains occupy separate panels with the same 0–13% scale. Stable aliases preserve the existing combined-cohort identifiers; gaps are intentional and no aliases are renumbered. The complete mapping is exported in pair_aliases.csv. Proteins use existing w=16 subsampling; hurricanes have no additional subsampling.

## Definition and interpretation (applies to every figure)
For X=A or B, each bar is 100 a_X/[2 n_X(n_X−1)], where n_X is the number of original vertices in the saved prepared input and a_X counts added auxiliary matching locations only, excluding original vertices. The finite bound follows the implemented construction: n_X vertex-centered spheres, n_X−1 input segments, and at most two interior intersections for each sphere–segment pair. Endpoint exclusion and merging of coincident curve parameters can reduce the count. Zero-length segments add no isolated intersections in this implementation. This is a finite combinatorial auxiliary-location bound, not an asymptotic O-expression evaluated with constant one, and not a configuration-state or edge bound.

Each curve uses its own denominator. These normalized percentages neither directly compare absolute auxiliary counts nor isolate the effect of data domain: input sizes, geometry, thresholds and preprocessing differ. Auxiliary locations are matching positions, not selectable simplified-output vertices. CPS-3F is omitted because it adds no auxiliary locations. Percentages are not search-completion indicators. Bar labels round to one decimal place; the CSV files retain integer numerators and denominators, exact rational percentage strings and full-precision numeric percentages.

Counts come from the single primary alpha=0.5 observation chosen in combined_results.csv, joined by exact attempt path to saved per-attempt evidence; timing repetitions and earlier resource diagnostics are not substituted. Older successful checkpoints with no per-curve completion flags are supported by their recorded counts, chain_graphs_completed=2, and successful configuration-graph route. Per-row count sources, JSON fields, hashes and completion evidence accompany the plotted values. No optimization, auxiliary-coordinate generation or validation was rerun.
'''
    (OUT/'captions.md').write_text(captions,encoding='utf-8')
    audit=dict(pairs=45,hurricane_pairs=33,protein_pairs=12,common_hurricane_pairs=31,common_protein_pairs=12,
               limited_pairs=df[df.subsequent_search_resource_limited][['pair','aA','aB','nA','nB','percentage_A','percentage_B','construction_completed_A','construction_completed_B']].to_dict('records'),
               optimization_executed=False,source_sha256={str(p.relative_to(ROOT)):digest(p) for p in [SOURCE,METRICS,ROOT/'cps_paper_algorithms.py',Path(__file__)]})
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2))
    return df

def coupling_size():
    """Add IC/CPS-2F plots; leave every existing auxiliary figure unchanged."""
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator
    import nbformat
    from nbclient import NotebookClient
    import shutil
    combined=pd.read_csv(SOURCE,float_precision='round_trip')
    metrics=pd.read_csv(METRICS,float_precision='round_trip')
    metrics['normalized_path']=metrics.attempt_path.str.replace('\\','/',regex=False)
    primary=combined[combined.alpha.eq(.5) & combined.observed.eq(True)]
    methods={'Independent continuous','Independent discrete','CPS-2F','CPS-3F'}
    common=[pair for pair,g in primary.groupby('pair') if len(g)==4 and set(g.method)==methods and g.validated.eq(True).all()]
    assert len(common)==43
    records=[]
    for pair in common:
        g=primary[primary.pair.eq(pair)]
        first=g.iloc[0]
        row=dict(pair=pair,alias=first.pair_alias,domain=first.domain,alpha=.5,nA=int(first.nA),nB=int(first.nB),
                 units=first.units,input_sha256=first.input_sha256,delta1=first.delta1,delta2=first.delta2,delta3=first.delta3,
                 absolute_coupling_tolerance=1e-9)
        assert g.input_sha256.nunique()==1 and g.delta3.nunique()==1 and first.delta3>0
        for method,tag in [('Independent continuous','IC'),('CPS-2F','CPS2F')]:
            r=g[g.method.eq(method)].iloc[0]
            key=r.attempt_path.replace('\\','/')
            matched=metrics[metrics.normalized_path.eq(key)]
            assert len(matched)==1
            metric=matched.iloc[0]
            vpath=ROOT/metric.validation_path
            validation=read(vpath)
            assert validation['status']=='validated'
            # Saved validation explicitly measures coupling of selected OUTPUT curves.
            distance=float(validation['coupling_distance'])
            assert math.isclose(distance,r.coupling_distance,rel_tol=1e-14,abs_tol=1e-12)
            assert math.isclose(distance,metric.coupling_distance,rel_tol=1e-14,abs_tol=1e-12)
            kA,kB=int(validation['kA']),int(validation['kB'])
            K=max(kA,kB)
            assert K==int(r.k)==int(metric.k)==validation['k']
            assert kA==int(r.kA)==int(metric.kA) and kB==int(r.kB)==int(metric.kB)
            feasible=distance<=row['delta3']+1e-9
            assert feasible==validation['coupling_feasible']==bool(r.coupling_feasible)
            assert validation['coupling_constraint_required']==(tag=='CPS2F')
            row.update({f'kA_{tag}':kA,f'kB_{tag}':kB,f'K_{tag}':K,
                        f'coupling_distance_{tag}':distance,f'normalized_coupling_{tag}':distance/row['delta3'],
                        f'coupling_satisfied_{tag}':feasible,f'attempt_path_{tag}':key,
                        f'validation_path_{tag}':vpath.relative_to(ROOT).as_posix(),
                        f'validation_sha256_{tag}':digest(vpath),f'status_{tag}':r.status})
        row['additional_K']=row['K_CPS2F']-row['K_IC']
        records.append(row)
    df=pd.DataFrame(records)
    df['domain_order']=df.domain.map({'protein':0,'hurricane':1})
    df['alias_number']=df.alias.str[1:].astype(int)
    df=df.sort_values(['domain_order','alias_number']).reset_index(drop=True)
    assert df.domain.value_counts().to_dict()=={'hurricane':31,'protein':12}
    actual=dict(IC_satisfies=int(df.coupling_satisfied_IC.sum()),IC_violates=int((~df.coupling_satisfied_IC).sum()),
                CPS2F_satisfies=int(df.coupling_satisfied_CPS2F.sum()),
                violating_equal_K=int((~df.coupling_satisfied_IC & df.additional_K.eq(0)).sum()),
                violating_plus_one=int((~df.coupling_satisfied_IC & df.additional_K.eq(1)).sum()),equal_K_overall=int(df.additional_K.eq(0).sum()))
    expected=dict(IC_satisfies=14,IC_violates=29,CPS2F_satisfies=43,violating_equal_K=24,violating_plus_one=5,equal_K_overall=38)
    discrepancies={k:dict(expected=expected[k],observed=v) for k,v in actual.items() if v!=expected[k]}
    print('COMMON-SET VALIDATION',actual,'DISCREPANCIES',discrepancies,flush=True)
    df.drop(columns=['domain_order','alias_number']).to_csv(OUT/'coupling_size_values.csv',index=False,float_format='%.17g')
    df[['alias','pair','domain','nA','nB']].to_csv(OUT/'coupling_size_aliases.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,'figure.facecolor':'white','axes.facecolor':'white'})
    ratio_max=float(df[['normalized_coupling_IC','normalized_coupling_CPS2F']].max().max())
    ycoupling=math.ceil(ratio_max*1.15*2)/2
    ysize=int(df[['K_IC','K_CPS2F']].max().max())+3
    handles=[Patch(facecolor=BLUE,label='IC (independent continuous)'),Patch(facecolor=GREEN,label='CPS-2F'),
             Line2D([0],[0],color='black',linestyle='--',label='Coupling threshold')]
    footer='Common validated four-method set. IC above the threshold remains a valid independent solution.\nCost means output size K, not runtime or graph states. Equal K does not imply identical selected vertices.'
    def rowplots(axes,data,compact):
        x=np.arange(len(data));w=.38
        for tag,offset,color in [('IC',-w/2,BLUE),('CPS2F',w/2,GREEN)]:
            axes[0].bar(x+offset,data['normalized_coupling_'+tag],w,color=color,zorder=3)
            axes[1].bar(x+offset,data['K_'+tag],w,color=color,zorder=3)
        axes[0].axhline(1,color='black',ls='--',lw=1.3,zorder=4)
        axes[0].set(title='Achieved coupling',ylabel='Output coupling distance / '+r'$\delta_3$',ylim=(0,ycoupling))
        axes[1].set(title='Output size',ylabel=r'$K=\max(|A^\prime|,|B^\prime|)$',ylim=(0,ysize))
        axes[1].yaxis.set_major_locator(MaxNLocator(integer=True,nbins=6))
        for i,r in enumerate(data.itertuples()):
            if r.additional_K!=0:axes[1].annotate(f'{r.additional_K:+d}',(i+w/2,r.K_CPS2F),xytext=(0,4),textcoords='offset points',ha='center',fontsize=9 if compact else 11)
        for ax in axes:
            labels=data.alias if compact else data.pair
            ax.set_xticks(x,labels,rotation=0 if compact else 90,fontsize=9 if compact else 10)
            ax.set_axisbelow(True);ax.grid(axis='y',color='#DDDDDD',linewidth=.6)
            ax.spines[['top','right']].set_visible(False);ax.margins(x=.015)
            ax.title.set_fontsize(12 if compact else 17)
    # Large format: precisely two side-by-side panels, complete pair names.
    fig,axes=plt.subplots(1,2,figsize=(40,10))
    fig.subplots_adjust(left=.035,right=.99,bottom=.29,top=.80,wspace=.12)
    rowplots(axes,df,False)
    for ax in axes:
        ax.axvspan(-.5,11.5,color='#EFEFEF',zorder=0)
        ax.axvline(11.5,color='#777777',lw=1)
        for x,label in [(5.5,'Proteins (12 pairs)'),(27,'Hurricanes (31 pairs)')]:
            ax.text(x,1.02,label,transform=ax.get_xaxis_transform(),ha='center',fontsize=12)
        ax.set_title(ax.get_title(),pad=35)
    fig.suptitle('IC vs. CPS-2F: achieved output coupling and output size',fontsize=23,y=.98)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.935),ncol=3,frameon=False,fontsize=15)
    for ext in ['pdf','png']:fig.savefig(OUT/f'coupling_size_common43_large.{ext}',dpi=300)
    plt.close(fig)
    # Publication layout: protein row followed by three manageable hurricane rows.
    h=df[df.domain.eq('hurricane')]
    groups=[('Proteins | 12 pairs',df[df.domain.eq('protein')])]
    groups += [(f'Hurricanes | group {i+1}/3',h.iloc[ix]) for i,ix in enumerate(np.array_split(np.arange(len(h)),3))]
    fig,axes=plt.subplots(4,2,figsize=(14,13))
    fig.subplots_adjust(left=.075,right=.985,bottom=.09,top=.86,hspace=.78,wspace=.22)
    for i,(title,data) in enumerate(groups):
        rowplots(axes[i],data,True)
        fig.text(.53,axes[i,0].get_position().y1+.037,title,ha='center',fontsize=12,weight='bold')
    fig.suptitle('IC vs. CPS-2F: achieved output coupling and output size\nCommon comparison set: 31 hurricane + 12 protein pairs',fontsize=16,y=.985)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.944),ncol=3,frameon=False,fontsize=11)
    fig.text(.53,.025,'Stable pair aliases are mapped in coupling_size_aliases.csv.',ha='center',fontsize=9)
    for ext in ['pdf','png']:fig.savefig(OUT/f'coupling_size_common43_publication.{ext}',dpi=300)
    plt.close(fig)
    caption='''# IC versus CPS-2F: coupling and output size

**Caption (both layouts).** Achieved discrete output coupling and simplification size at alpha=0.5 for the common comparison set of 43 pairs (31 hurricane and 12 protein pairs) with validated primary outputs from all four methods. Blue bars denote independent continuous simplification (IC); green bars denote CPS-2F. Left panels show the recorded discrete Frechet distance between the two simplified OUTPUT curves, divided by delta3. The dashed line at one is the coupling threshold. IC values above one violate the joint coupling constraint but remain valid independent solutions. Right panels show actual K=max(|A'|,|B'|), with increases annotated above CPS-2F bars. Here cost means output size, not runtime or configuration-state count. Equal K does not imply identical vertex selections. All pairs are retained, including those where IC already satisfies coupling. Protein inputs use the existing w=16 subsampling; hurricane tracks have no additional subsampling. The publication layout keeps coupling left and size right in every row, with consistent scales for each metric and stable aliases mapped in coupling_size_aliases.csv.

**Data and interpretation.** One primary observation per pair and method is taken from combined_results.csv; per_attempt_metrics.csv is joined by exact normalized attempt path, then distances and counts are checked against the saved validation records. Timing repetitions, early diagnostics and resource-limited pairs are excluded. Coupling satisfaction uses the existing absolute tolerance d <= delta3 + 1e-9 in each input's coordinate units; bars show the unmodified ratio d/delta3 without clamping. Feasibility validation is distinct from independent proof of optimality. No optimization, geometry, coupling distances or fidelity predicates were recomputed. The comparison holds continuous input-fidelity constraints fixed and adds the output coupling constraint. It is not a runtime comparison.
'''
    caption+='\n**Saved-data check.** '+json.dumps(actual)+'.\n\nDiscrepancies against the requested expectations: '+(json.dumps(discrepancies) if discrepancies else 'none')+'.\n'
    (OUT/'coupling_size_captions.md').write_text(caption,encoding='utf-8')
    audit=dict(actual=actual,expected=expected,discrepancies=discrepancies,absolute_tolerance=1e-9,ratio_max=ratio_max,
               coupling_ylim=[0,ycoupling],size_ylim=[0,ysize],optimization_executed=False,
               source_sha256={str(p.relative_to(ROOT)):digest(p) for p in [SOURCE,METRICS,ROOT/'notebooks_paper/graph_only_v3/common.py']})
    (OUT/'coupling_size_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    notebook=ROOT/'notebooks_paper/auxiliary_bound_percentages.ipynb'
    backup=OUT/'qa/notebook_before_coupling_size.ipynb'
    if not backup.exists():shutil.copy2(notebook,backup)
    nb=nbformat.read(notebook,as_version=4)
    if not any(c.metadata.get('coupling_size_section') for c in nb.cells):
        cells=[nbformat.v4.new_markdown_cell('## IC versus CPS-2F: achieved coupling and output size\nCommon validated set: **31 hurricane + 12 protein pairs**, alpha=0.5. Every row retains coupling on the left and actual K on the right. These cells read saved charts and values only.'),
               nbformat.v4.new_code_cell("import json\ncoupling_values = pd.read_csv(PLOT_DIR/'coupling_size_values.csv')\ndisplay(pd.DataFrame([json.loads((PLOT_DIR/'coupling_size_audit.json').read_text())['actual']]))\ndisplay(Image(filename=str(PLOT_DIR/'coupling_size_common43_large.png'), width=1600))"),
               nbformat.v4.new_code_cell("display(Image(filename=str(PLOT_DIR/'coupling_size_common43_publication.png'), width=1100))\ndisplay(pd.read_csv(PLOT_DIR/'coupling_size_aliases.csv'))\ndisplay(Markdown((PLOT_DIR/'coupling_size_captions.md').read_text()))")]
        for cell in cells:cell.metadata['coupling_size_section']=True
        nb.cells.extend(cells)
    NotebookClient(nb,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    nbformat.write(nb,notebook)
    return df

if __name__=='__main__':
    import sys
    coupling_size() if '--coupling-size' in sys.argv else build()
