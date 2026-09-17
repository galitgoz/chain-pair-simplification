"""Stage 1 protocol and reporting adapters. Optimization stays in existing modules."""
from pathlib import Path
import os,sys,json,hashlib,time,itertools,traceback,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output/paper_v1'
sys.path.insert(0,str(ROOT))
os.environ['NUMBA_CACHE_DIR']=str(OUT/'cache')
os.environ['MPLCONFIGDIR']=str(OUT/'cache/matplotlib')
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import pandas as pd
METHODS=['Independent continuous','Independent discrete','CPS-2F','CPS-3F']
COLORS=dict(zip(METHODS,['#0072B2','#56B4E9','#009E73','#D55E00']))
LIMITS=dict(seconds=8.,max_states=30_000,max_transitions=500_000,max_chain_edges=200_000)
PROTOCOL=dict(alphas=[.5,1.],protein_existing_stride=16,endpoints='retain first and last input indices',
    hard_seconds=20.,startup_seconds=90.,memory_mib=512,solver_limits=LIMITS,repetitions=1,
    selection='two smallest aligned supplied protein pairs at fixed existing w=16; smallest hurricane pair per basin; size ties broken by identifiers',
    thresholds='delta1=alpha*mean_edge(A); delta2=alpha*mean_edge(B); delta3=unrounded discrete Frechet(A,B)',
    scope='Stage 1 pipeline pilot, not final paper parameter selection')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def prepare_protocol():
    from curve_algorithms import discrete_frechet
    folder=ROOT/'data/decimated protein backbones'
    curves=pd.read_csv(folder/'decimated_curves.csv');pairs=pd.read_csv(folder/'decimated_pairs.csv')
    eligible=pairs[pairs.frame.eq('aligned')&pairs.w.eq(PROTOCOL['protein_existing_stride'])].copy()
    eligible['size_max']=eligible[['n_A','n_B']].max(axis=1);eligible['size_sum']=eligible.n_A+eligible.n_B
    eligible=eligible.sort_values(['size_max','size_sum','chain_A','chain_B'])
    storms=pd.read_csv(ROOT/'data/hurricanes/pairs.csv')
    storms['size_max']=storms[['n_A','n_B']].max(axis=1);storms['size_sum']=storms.n_A+storms.n_B
    chosen_storms=storms.sort_values(['size_max','size_sum','curve_A','curve_B']).groupby('basin',sort=True).head(1)
    candidates=[];inputs=[]
    chainmeta=pd.read_csv(folder/'decimated_chains.csv')
    for r in eligible.itertuples():
        candidates.append(dict(domain='protein',pair=r.chain_A+'__'+r.chain_B,nA=r.n_A,nB=r.n_B,
            eligible=True,selection_rank=len(candidates)+1,selected=len(inputs)<2))
        if len(inputs)>=2:continue
        A=curves[curves.chain.eq(r.chain_A)&curves.frame.eq('aligned')&curves.w.eq(r.w)].sort_values('t')
        B=curves[curves.chain.eq(r.chain_B)&curves.frame.eq('aligned')&curves.w.eq(r.w)].sort_values('t')
        gaps=[]
        for chain in [r.chain_A,r.chain_B]:
            meta=chainmeta[chainmeta.chain.eq(chain)&chainmeta.frame.eq('deposited')&chainmeta.w.eq(1)].iloc[0]
            gaps.append(str(meta.gaps))
        inputs.append(dict(domain='protein',pair=r.chain_A+'__'+r.chain_B,id_A=r.chain_A,id_B=r.chain_B,
            A=A[['x','y','z']].to_numpy(float),B=B[['x','y','z']].to_numpy(float),units='angstrom',
            provenance='supplied decimated_curves.csv; coordinates reused exactly',
            subsampling=f'existing w={r.w}; no further sampling; stored rows '+','.join(map(str,A.source_row)),
            gaps_A=gaps[0],gaps_B=gaps[1],alignment='supplied common aligned frame (existing Kabsch); no new fit',
            extent_policy='entire supplied subsampled curve; not a fragment; not a full-resolution backbone',
            source_A=str((folder/'decimated_curves.csv').relative_to(ROOT)),source_B=str((folder/'decimated_curves.csv').relative_to(ROOT)),basin=''))
    keys=set(zip(chosen_storms.curve_A,chosen_storms.curve_B))
    for r in storms.itertuples():candidates.append(dict(domain='hurricane',pair=r.curve_A+'__'+r.curve_B,nA=r.n_A,nB=r.n_B,
        eligible=True,selection_rank=None,selected=(r.curve_A,r.curve_B) in keys))
    for r in chosen_storms.itertuples():
        sources=[ROOT/f'data/hurricanes/{x}.txt' for x in [r.curve_A,r.curve_B]]
        inputs.append(dict(domain='hurricane',pair=r.curve_A+'__'+r.curve_B,id_A=r.curve_A,id_B=r.curve_B,
            A=np.loadtxt(sources[0]),B=np.loadtxt(sources[1]),units='projected km',
            provenance='supplied HURDAT2-derived coordinates; data/hurricanes/README.md; all supplied track vertices',
            subsampling='none added; upstream filtering and consecutive-duplicate removal documented in dataset README',
            gaps_A='timestamps absent; temporal gaps unknown',gaps_B='timestamps absent; temporal gaps unknown',
            alignment='unchanged shared basin azimuthal-equidistant frame; no independent alignment',
            extent_policy='entire supplied track; not a fragment',source_A=str(sources[0].relative_to(ROOT)),
            source_B=str(sources[1].relative_to(ROOT)),basin=r.basin))
    manifests=[];parameters=[]
    for item in inputs:
        A=item.pop('A');B=item.pop('B');assert np.isfinite(A).all() and np.isfinite(B).all()
        path=OUT/'inputs'/f"{item['domain']}_{item['pair']}.npz"
        np.savez(path,A=A,B=B)
        item.update(nA=len(A),nB=len(B),dimensions=A.shape[1],input_file=str(path.relative_to(ROOT)),
            input_sha256=sha(path),source_A_sha256=sha(ROOT/item['source_A']),source_B_sha256=sha(ROOT/item['source_B']),
            endpoint_policy=PROTOCOL['endpoints'],pilot_selected_before_optimization=True)
        manifests.append(item)
        sa=float(np.linalg.norm(np.diff(A,axis=0),axis=1).mean());sb=float(np.linalg.norm(np.diff(B,axis=0),axis=1).mean())
        d3=float(discrete_frechet(A,B))
        for alpha in PROTOCOL['alphas']:
            parameters.append(dict(domain=item['domain'],pair=item['pair'],case=f"{item['domain']}:{item['pair']}:a{alpha:g}",
                nA=len(A),nB=len(B),units=item['units'],alpha=alpha,s_A=sa,s_B=sb,
                delta1=alpha*sa,delta2=alpha*sb,delta3=d3,input_file=item['input_file']))
    manifest=pd.DataFrame(manifests);params=pd.DataFrame(parameters)
    shared=[]
    for a,b in itertools.combinations(manifests,2):
        chains=set([a['id_A'],a['id_B']])&set([b['id_A'],b['id_B']])
        structures=set(x.split('.')[0] for x in [a['id_A'],a['id_B']])&set(x.split('.')[0] for x in [b['id_A'],b['id_B']]) if a['domain']==b['domain']=='protein' else set()
        if chains or structures:shared.append(dict(pair_1=a['pair'],pair_2=b['pair'],shared_ids=';'.join(sorted(chains)),shared_structures=';'.join(sorted(structures))))
    manifest.to_csv(OUT/'input_manifest.csv',index=False);params.to_csv(OUT/'pilot_parameters.csv',index=False)
    pd.DataFrame(candidates).to_csv(OUT/'selection_audit.csv',index=False)
    pd.DataFrame(shared,columns=['pair_1','pair_2','shared_ids','shared_structures']).to_csv(OUT/'shared_entities.csv',index=False)
    frozen=dict(**PROTOCOL,configurations=len(params)*len(METHODS),method_names=METHODS,
        parameter_sha256=sha(OUT/'pilot_parameters.csv'),input_manifest_sha256=sha(OUT/'input_manifest.csv'),
        algorithms={name:sha(ROOT/name) for name in ['curve_algorithms.py','cps_paper_algorithms.py']})
    (OUT/'protocol.json').write_text(json.dumps(frozen,indent=2))
    return manifest,params

def solve(A,B,d1,d2,d3,method,certificate=True):
    from curve_algorithms import shortest_independent
    from cps_paper_algorithms import independent_discrete,solve_cps
    if method.startswith('Independent'):
        fn=shortest_independent if method=='Independent continuous' else independent_discrete
        a,b=fn(A,d1),fn(B,d2)
        return dict(status='optimal',route='independent_optimization',indices_A=a,indices_B=b)
    raw=solve_cps(A,B,d1,d2,d3,method,certificate=certificate,**LIMITS)
    status=raw['status']
    raw['route']={'optimal_independent_certificate':'independent_certificate',
        'optimal_alternate_certificate':'alternate_certificate','optimal_configuration_graph':'configuration_graph'}.get(status,'none')
    raw['status']='optimal' if status.startswith('optimal') else 'infeasible' if status.startswith('infeasible') else status
    raw['raw_status']=status
    return raw

def measures(A,B,result,d1,d2,d3,method):
    from curve_algorithms import continuous_distance,continuous_decision,discrete_frechet
    from cps_paper_algorithms import discrete_coupling
    ia=np.asarray(result['indices_A'],int);ib=np.asarray(result['indices_B'],int)
    for P,ids in [(A,ia),(B,ib)]:
        assert ids[0]==0 and ids[-1]==len(P)-1 and (np.diff(ids)>0).all()
    a,b=A[ia],B[ib];coupling,d=discrete_coupling(a,b)
    assert tuple(coupling[0])==(0,0) and tuple(coupling[-1])==(len(a)-1,len(b)-1)
    moves=np.diff(coupling,axis=0);assert ((moves>=0)&(moves<=1)).all() and (moves.sum(axis=1)>=1).all()
    discrete=method in ['Independent discrete','CPS-3F']
    for P,Q,delta in [(A,a,d1),(B,b,d2)]:
        assert continuous_decision(P,Q,delta)
        if discrete:assert discrete_frechet(P,Q)<=delta+1e-9
    if method.startswith('CPS'):assert d<=d3+1e-9
    ca=continuous_distance(A,a);cb=continuous_distance(B,b)
    return dict(kA=len(a),kB=len(b),k=max(len(a),len(b)),compression_removed_pct=100*(1-(len(a)+len(b))/(len(A)+len(B))),
        continuous_A_lower=ca[0],continuous_A_upper=ca[1],continuous_B_lower=cb[0],continuous_B_upper=cb[1],
        discrete_fidelity_A=float(discrete_frechet(A,a)),discrete_fidelity_B=float(discrete_frechet(B,b)),
        coupling_distance=d,coupling_budget_ratio=d/d3 if d3 else np.nan,coupling_feasible=d<=d3+1e-9,
        unequal_output_lengths=len(a)!=len(b),repeated_coupling_index=bool((moves==0).any()),
        indices_A=ia.tolist(),indices_B=ib.tolist(),coupling=coupling.tolist(),constraints_valid=True)

def warmup(method):
    P=np.array([[0.,0.],[1,0],[2,0]])
    solve(P,P,.5,.5,1.,method,certificate=False)

def run_pilot(domain):
    import psutil
    params=pd.read_csv(OUT/'pilot_parameters.csv');protocol=json.loads((OUT/'protocol.json').read_text())
    assert sha(OUT/'pilot_parameters.csv')==protocol['parameter_sha256']
    assert sha(OUT/'input_manifest.csv')==protocol['input_manifest_sha256']
    for name,digest in protocol['algorithms'].items():assert sha(ROOT/name)==digest
    for r in pd.read_csv(OUT/'input_manifest.csv').itertuples():assert sha(ROOT/r.input_file)==r.input_sha256
    if (OUT/f'{domain}_results.csv').exists():
        raise RuntimeError('Saved pilot exists. Read saved results; do not overwrite or count cached results as repetitions.')
    gate=pd.read_csv(OUT/'validation.csv')
    blocked=set(gate.loc[~gate.outcome.eq('PASS'),'affected_method'])
    rows=[];paths={}
    for param in params[params.domain.eq(domain)].to_dict('records'):
        for method in METHODS:
            row={**param,'method':method,'repetition':1,'cached':False,'runtime_kind':'diagnostic single execution',
                'status':'not_run','route':'none','failure_reason':'','hard_seconds':PROTOCOL['hard_seconds'],'memory_limit_mib':PROTOCOL['memory_mib']}
            if method in blocked or 'all' in blocked:
                row.update(status='blocked_correctness',failure_reason='See validation.csv substantive check failure');rows.append(row);continue
            stem=param['case'].replace(':','_')+'_'+method.replace(' ','_')
            task=OUT/'workers'/f'{stem}.json';result_path=OUT/'workers'/f'{stem}.result.json';ready=OUT/'workers'/f'{stem}.ready'
            task.write_text(json.dumps(dict(param=param,method=method,result_path=str(result_path),ready=str(ready))))
            log=open(OUT/'logs'/f'{stem}.log','w',encoding='utf-8')
            start=time.perf_counter();process=subprocess.Popen([sys.executable,str(Path(__file__).with_name('pilot_worker.py')),str(task)],cwd=ROOT,
                stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
            peak=0;run_start=None;limit=None
            while process.poll() is None:
                try:
                    tree=psutil.Process(process.pid).children(recursive=True)+[psutil.Process(process.pid)]
                    rss=0
                    for child in tree:
                        try:rss+=child.memory_info().rss
                        except psutil.NoSuchProcess:pass
                    peak=max(peak,rss)
                except psutil.NoSuchProcess:pass
                if ready.exists() and run_start is None:run_start=time.perf_counter()
                if peak>PROTOCOL['memory_mib']*1024**2:limit='memory_limit'
                elif run_start is not None and time.perf_counter()-run_start>PROTOCOL['hard_seconds']:limit='time_limit'
                elif run_start is None and time.perf_counter()-start>PROTOCOL['startup_seconds']:limit='startup_time_limit'
                if limit:
                    for child in reversed(tree):
                        try:child.kill()
                        except psutil.NoSuchProcess:pass
                    process.wait();break
                time.sleep(.05)
            log.close();row['process_wall_seconds']=time.perf_counter()-start;row['peak_process_rss_mib']=peak/1024**2
            row['measured_execution_count']=1 if run_start is not None else 0
            if limit:row.update(status=limit,failure_reason='Isolated worker resource cap; not a proof of infeasibility')
            elif result_path.exists():row.update(json.loads(result_path.read_text()))
            else:row.update(status='software_error',failure_reason=f'Worker exit {process.returncode}; see {stem}.log')
            paths[param['case']+'|'+method]={k:row.pop(k) for k in ['indices_A','indices_B','coupling'] if k in row}
            rows.append(row)
            print(param['case'],method,row['status'],row.get('route'),flush=True)
            pd.DataFrame(rows).to_csv(OUT/f'{domain}_results.csv',index=False)
            (OUT/f'{domain}_paths.json').write_text(json.dumps(paths,indent=2))
            print(f"{param['pair']} alpha={param['alpha']} {method}: {row['status']} k={row.get('k','NA')}",flush=True)
    frame=pd.DataFrame(rows)
    for i,r in frame.iterrows():
        baseline='Independent continuous' if r.method=='CPS-2F' else 'Independent discrete' if r.method=='CPS-3F' else r.method
        base=frame[(frame.case==r.case)&(frame.method==baseline)]
        frame.loc[i,'baseline_method']=baseline
        frame.loc[i,'baseline_missing']=base.empty or base.iloc[0].status!='optimal'
        if r.status=='optimal' and not frame.loc[i,'baseline_missing']:
            frame.loc[i,'vertex_cost_vs_matching_baseline']=r.k-base.iloc[0].k
    frame.to_csv(OUT/f'{domain}_results.csv',index=False)
    return frame

def compact_results(frame):
    cols=['pair','alpha','method','kA','kB','k','vertex_cost_vs_matching_baseline','status','route','solve_seconds']
    return frame.reindex(columns=cols)

def plot_domain(domain):
    import matplotlib.pyplot as plt
    frame=pd.read_csv(OUT/f'{domain}_results.csv');params=pd.read_csv(OUT/'pilot_parameters.csv')
    folder=OUT/'figures'/('proteins' if domain=='protein' else 'hurricanes');folder.mkdir(exist_ok=True)
    with plt.rc_context({'figure.facecolor':'white','axes.facecolor':'white','text.color':'black','axes.labelcolor':'black','font.size':10}):
        cases=list(frame.case.drop_duplicates());x=np.arange(len(cases))
        fig,axes=plt.subplots(1,2,figsize=(12,4.7))
        for i,m in enumerate(METHODS):
            r=frame[frame.method.eq(m)].set_index('case').reindex(cases)
            for ax,field in zip(axes,['k','coupling_budget_ratio']):
                values=pd.to_numeric(r.get(field,pd.Series(np.nan,index=r.index)),errors='coerce')
                ax.plot(x+(i-1.5)*.04,values,'o-',color=COLORS[m],label=m,ms=5,lw=1)
                for j,v in enumerate(values):
                    if pd.isna(v):ax.text(j+(i-1.5)*.045,.03,'missing',rotation=90,color=COLORS[m],transform=ax.get_xaxis_transform(),fontsize=7)
        labels=[c.split(':')[1].replace('__',' / ')+'\nα='+c.split(':')[-1][1:] for c in cases]
        for ax in axes:ax.set_xticks(x,labels,rotation=20,ha='right');ax.grid(alpha=.15)
        axes[0].set_ylabel(r'$k=\max(|A^\prime|,|B^\prime|)$');axes[1].set_ylabel('Discrete coupling / δ₃')
        axes[1].axhline(1,color='black',ls='--',lw=1)
        for ax,field in zip(axes,['k','coupling_budget_ratio']):
            if field not in frame or pd.to_numeric(frame[field],errors='coerce').isna().all():
                ax.set_ylim(0,1.2);ax.set_yticks([])
                ax.text(.5,.55,'No verified outputs within this run\nSee failure reasons in the results table',ha='center',va='center',transform=ax.transAxes)
        fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',ncol=4,frameon=False)
        fig.tight_layout(rect=[0,0,1,.9]);fig.savefig(folder/'paired_comparison.png',dpi=180,bbox_inches='tight');fig.savefig(folder/'paired_comparison.svg',bbox_inches='tight')
        plt.show();plt.close(fig)
        paths=json.loads((OUT/f'{domain}_paths.json').read_text());r=frame.iloc[0]
        inp=np.load(ROOT/r.input_file);A,B=inp['A'],inp['B']
        fig=plt.figure(figsize=(12,3.6))
        allpoints=np.vstack([A,B]);mid=.5*(allpoints.min(0)+allpoints.max(0));half=.55*np.ptp(allpoints,axis=0).max()
        for i,m in enumerate(METHODS):
            ax=fig.add_subplot(1,4,i+1,projection='3d' if A.shape[1]==3 else None)
            for P,color,side in [(A,'#333333','A'),(B,'#999999','B')]:
                ax.plot(*P.T,color=color,alpha=.35,lw=1)
                p=paths.get(r.case+'|'+m,{});ids=p.get('indices_'+side)
                if ids is not None:ax.plot(*P[ids].T,'o',color=COLORS[m],ms=3,ls='-' if side=='A' else '--')
            ax.set_xlim(mid[0]-half,mid[0]+half);ax.set_ylim(mid[1]-half,mid[1]+half)
            if A.shape[1]==3:ax.set_zlim(mid[2]-half,mid[2]+half);ax.set_box_aspect([1,1,1])
            else:ax.set_aspect('equal')
            ax.set_title(m+(' — missing' if not p else ''),fontsize=9);ax.tick_params(labelsize=7)
            ax.set_xlabel('x ('+r.units+')',fontsize=7);ax.set_ylabel('y ('+r.units+')',fontsize=7)
        fig.suptitle(f"{r.pair.replace('__',' / ')}; α={r.alpha:g}; δ₁={r.delta1:.4g}, δ₂={r.delta2:.4g}, δ₃={r.delta3:.4g} {r.units}\nOriginal coordinates; faint inputs, solid A′ / dashed B′; pipeline pilot",fontsize=10)
        fig.tight_layout(rect=[0,0,1,.83]);fig.savefig(folder/'geometry_example.png',dpi=180,bbox_inches='tight');plt.show();plt.close(fig)
