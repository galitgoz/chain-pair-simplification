"""Expanded cohort selection and saved-result presentation; shared solvers are reused."""
from pathlib import Path
import sys, json, time, argparse, hashlib
import numpy as np
import pandas as pd
from common import ROOT, atomic, sha, stamp, BUDGETS
import controlled_followup as runner

OUT=ROOT/'output/paper_v1/expanded_comparison'
OLD=ROOT/'output/paper_v1/unified_pilot_results'
METHODS=['Independent continuous','Independent discrete','CPS-2F','CPS-3F']
COLORS=dict(zip(METHODS,['#377eb8','#984ea3','#e66101','#1b9e77']))
LIMITS=dict(seconds=55.,max_states=100000,max_transitions=20000000,max_chain_edges=800000)
RULE='''# Frozen selection rule

Use supplied aligned w=16 protein pairs listed in decimated_pairs.csv and supplied
hurricane pairs listed in hurricanes/pairs.csv; retain the four original pairs.
No solver outcomes enter selection. Compute features from input coordinates only:
log(max(nA,nB)), log(max(nA,nB)/min(nA,nB)), mean chord/path-length ratio,
and centroid separation divided by mean path length. Standardize each feature
over the eligible candidate pool (zero-variance features contribute zero).
Select three additions per domain by farthest-first distance to the retained and
already selected pairs. First minimize reuse of individual structures/storms;
then maximize the minimum feature distance; break ties by pair identifier.
Protein identity for reuse is PDB structure, not chain: all supplied aligned
protein pairs share 1o7j.a, so that reuse is unavoidable. Only one additional
4eca chain will be selected when unused structures remain. All supplied protein
curves have 325 prepared points before the existing w=16 sampling, yielding 22
points; this is not full-resolution backbone data. New hurricanes use every
supplied row, with no alignment, cropping, or subsampling. Pair orientation is
the supplied orientation. Duplicate unordered pairs are removed deterministically.
Freeze all arrays, source hashes, numerical thresholds and configurations before
any new optimization. No selection by success, compression, or runtime.
'''

def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def rows():return read(OUT/'configuration_manifest.json')['rows']
def configure():runner.COHORT=OUT
def feature(A,B):
    lengths=[np.linalg.norm(np.diff(P,axis=0),axis=1).sum() for P in [A,B]]
    return [np.log(max(len(A),len(B))),np.log(max(len(A),len(B))/min(len(A),len(B))),
        np.mean([np.linalg.norm(P[-1]-P[0])/max(l,1e-12) for P,l in zip([A,B],lengths)]),
        np.linalg.norm(A.mean(0)-B.mean(0))/max(np.mean(lengths),1e-12)]

def prepare():
    runner.cohort_active_check();OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'configuration_manifest.json').exists():
        print('Manifest already frozen; preserving it.');return
    # This rule is written before computing rankings or running any solver.
    (OUT/'selection_rule.md').write_text(RULE,encoding='utf-8')
    old=pd.read_csv(OLD/'results.csv');assert len(old)==32
    curves=pd.read_csv(ROOT/'data/decimated protein backbones/decimated_curves.csv')
    pp=pd.read_csv(ROOT/'data/decimated protein backbones/decimated_pairs.csv')
    hp=pd.read_csv(ROOT/'data/hurricanes/pairs.csv')
    pools={};source_files=[ROOT/'data/decimated protein backbones/decimated_curves.csv',ROOT/'data/decimated protein backbones/decimated_pairs.csv',ROOT/'data/hurricanes/pairs.csv']
    for domain in ['protein','hurricane']:
        pool=[];seen=set()
        source=pp[(pp.frame=='aligned')&(pp.w==16)] if domain=='protein' else hp
        for r in source.sort_values(['chain_A','chain_B'] if domain=='protein' else ['curve_A','curve_B']).to_dict('records'):
            a,b=(r['chain_A'],r['chain_B']) if domain=='protein' else (r['curve_A'],r['curve_B'])
            key=tuple(sorted([a,b]))
            if key in seen:continue
            seen.add(key)
            if domain=='protein':
                def load(c):return curves[(curves.chain==c)&(curves.frame=='aligned')&(curves.w==16)].sort_values('t')[['x','y','z']].to_numpy(float)
                A,B=load(a),load(b)
            else:
                A,B=[np.loadtxt(ROOT/'data/hurricanes'/f'{c}.txt') for c in [a,b]]
            assert A.ndim==B.ndim==2 and len(A)>1 and len(B)>1
            pool.append(dict(domain=domain,pair=a+'__'+b,A=A,B=B,individuals=[a.split('.')[0],b.split('.')[0]] if domain=='protein' else [a,b],
                features=feature(A,B),basin=r.get('basin'),source_record=r))
        pools[domain]=pool
    selected=[];audit=[]
    for domain,pool in pools.items():
        retained=set(old[old.domain==domain].pair)
        chosen=[p for p in pool if p['pair'] in retained];assert len(chosen)==2
        X=np.asarray([p['features'] for p in pool]);sd=X.std(0);sd[sd<1e-12]=1
        for p,x in zip(pool,(X-X.mean(0))/sd):p['standardized_features']=x.tolist()
        used={i for p in chosen for i in p['individuals']}
        for _ in range(3):
            candidates=[p for p in pool if p['pair'] not in {q['pair'] for q in chosen}]
            if not candidates:break
            def score(p):
                reuse=sum(i in used for i in p['individuals'])
                distance=min(np.linalg.norm(np.array(p['standardized_features'])-q['standardized_features']) for q in chosen)
                return reuse,-distance,p['pair']
            pick=min(candidates,key=score);pick['selection_score']=score(pick);chosen.append(pick);used.update(pick['individuals'])
        selected.extend(chosen)
        for p in pool:audit.append({k:v for k,v in p.items() if k not in ['A','B','source_record']}|dict(selected=p['pair'] in {q['pair'] for q in chosen}))
    atomic(OUT/'candidate_selection.json',dict(rule_sha256=sha(OUT/'selection_rule.md'),candidates=audit,**stamp()))
    algorithms=read(OLD/'cohort.json')['algorithms']
    for n,h in algorithms.items():assert sha(ROOT/n)==h
    template=read(ROOT/old.iloc[0].attempt_path/'configuration.json')
    manifest=[];inputs=[];(OUT/'inputs').mkdir(exist_ok=True)
    for r in old.to_dict('records'):
        c=read(ROOT/r['attempt_path']/'configuration.json')
        assert sha(ROOT/c['input_file'])==c['input_sha256']
        manifest.append(dict(**{k:r[k] for k in ['configuration_id','plan_order','domain','pair','nA','nB','units','input_file','input_sha256','method','alpha','delta1','delta2','delta3']},
            configuration=c,reuse_path=r['attempt_path'],original_origin=r['origin']))
    from curve_algorithms import discrete_frechet
    for p in selected:
        domain,pair,A,B=p['domain'],p['pair'],p['A'],p['B'];prior=old[old.pair==pair]
        if len(prior):
            c=read(ROOT/prior.iloc[0].attempt_path/'configuration.json');data=np.load(ROOT/c['input_file'])
            assert np.array_equal(A,data['A']) and np.array_equal(B,data['B'])
        else:
            path=OUT/'inputs'/f'{domain}_{pair}.npz';np.savez(path,A=A,B=B)
            c={k:v for k,v in template.items() if k not in ['code_version','configuration_id','diagnostic','case']}
            c.update(domain=domain,pair=pair,nA=len(A),nB=len(B),units='angstrom' if domain=='protein' else 'projected km',
                input_file=str(path.relative_to(ROOT)),input_sha256=sha(path),s_A=float(np.linalg.norm(np.diff(A,axis=0),axis=1).mean()),
                s_B=float(np.linalg.norm(np.diff(B,axis=0),axis=1).mean()),delta3=float(discrete_frechet(A,B)))
            for alpha in [.5,1.]:
                for method in METHODS:
                    order=len(manifest)+1;cid=f'{domain}:{pair}:a{alpha:g}|{method}'
                    config={**c,'alpha':alpha,'delta1':alpha*c['s_A'],'delta2':alpha*c['s_B'],'method':method,'certificate':False,
                        'graph_limits':LIMITS,'phase_budgets':BUDGETS,'memory_mib':512,'algorithms':algorithms,
                        'configuration_id':cid,'case':cid.split('|')[0],'diagnostic':f'expanded{order:02d}'}
                    manifest.append(dict(configuration_id=cid,plan_order=order,**{k:config[k] for k in ['domain','pair','nA','nB','units','input_file','input_sha256','method','alpha','delta1','delta2','delta3']},configuration=config,reuse_path=None,original_origin=None))
        inputrow=dict(domain=domain,pair=pair,nA=len(A),nB=len(B),input_file=c['input_file'],input_sha256=c['input_sha256'],
            inherited=bool(len(prior)),units=c['units'],subsampling='existing w=16; 325 supplied points -> 22' if domain=='protein' else 'none; full supplied track',
            alignment='existing Kabsch to 1o7j.a' if domain=='protein' else 'supplied common basin coordinates; no alignment',
            basin=p['basin'],individuals=p['individuals'],features=p['features'],source_record=p['source_record'])
        inputs.append(inputrow)
    atomic(OUT/'input_manifest.json',dict(inputs=inputs,source_hashes={str(p.relative_to(ROOT)):sha(p) for p in source_files},**stamp()))
    pd.DataFrame([{k:v for k,v in r.items() if k!='source_record'} for r in inputs]).to_csv(OUT/'input_manifest.csv',index=False)
    atomic(OUT/'configuration_manifest.json',dict(rows=manifest,algorithms=algorithms,rule_sha256=sha(OUT/'selection_rule.md'),**stamp()))
    pd.DataFrame([{k:v for k,v in r.items() if k!='configuration'}|{'configuration_json':json.dumps(r['configuration'])} for r in manifest]).to_csv(OUT/'run_manifest.csv',index=False)
    atomic(OUT/'provenance.json',dict(inherited=32,new_configurations=len(manifest)-32,algorithms=algorithms,
        code={n:sha(Path(__file__).with_name(n)) for n in ['common.py','worker.py','supervisor.py','controlled_followup.py','expanded_comparison.py']},
        instrumentation='At completion of each main single-curve graph, save dense and active state counts, edges including waits, extended and auxiliary counts. No inner-loop changes; small boundary observation overhead included in solver time.',
        policy='Shared external time/memory budgets, disclosed method-specific internal edge limits. Independent discrete and CPS-3F preprocessing: 2000000; main CPS graphs: 800000.',**stamp()))
    for r in manifest:
        slot=OUT/'slots'/f"{r['plan_order']:02d}";slot.mkdir(parents=True,exist_ok=True)
        atomic(slot/'slot.json',dict(configuration_id=r['configuration_id'],plan_order=r['plan_order'],reuse_path=r['reuse_path']))
        if not r['reuse_path']:atomic(slot/'configuration.json',r['configuration'])
    save();print('FROZEN',len(inputs),'pairs;',len(manifest)-32,'new configurations',flush=True)

def save():
    configure();frame=pd.DataFrame([runner.cohort_result(r) for r in rows()])
    # The shared legacy aggregator labels a validation timeout as incomplete;
    # retain that validation status while also exposing the limiting resource.
    frame['termination_reason']=frame['status']
    for i,r in frame.iterrows():
        if r.attempt_path:
            p=ROOT/r.attempt_path/'supervisor.json'
            sup=read(p) if p.exists() else {}
            if sup.get('limiting_phase'):
                frame.at[i,'termination_reason']=sup['status']
                frame.at[i,'resource_limited']=True
    frame['validation_recovery_path']=''
    recovery_indices=list(OUT.glob('validation_recovery*/index.json'))
    for recovery_index in recovery_indices:
        info=read(recovery_index);d=ROOT/info['attempt_directory'];v=read(d/'validation.json')
        mask=frame.configuration_id.eq(info['configuration_id'])
        if v['status']=='validated':
            i=frame.index[mask][0];source=ROOT/frame.at[i,'attempt_path']/'raw_solver.json'
            assert v['raw_checkpoint_sha256']==sha(source)
            frame.at[i,'original_validation_seconds']=frame.at[i,'validation_seconds']
            frame.at[i,'original_limiting_phase']=frame.at[i,'limiting_phase']
            frame.at[i,'limiting_phase']=None
            for key,value in dict(status='validated',validated=True,resource_limited=False,termination_reason='validated',
                validation_status='validated',validation_seconds=v['validation_seconds'],validation_recovery_path=str(d.relative_to(ROOT)),
                kA=v['kA'],kB=v['kB'],k=v['k'],coupling_distance=v['coupling_distance'],coupling_feasible=v['coupling_feasible'],
                retained_A=v['kA']/frame.at[i,'nA'],retained_B=v['kB']/frame.at[i,'nB'],
                retained_total=(v['kA']+v['kB'])/(frame.at[i,'nA']+frame.at[i,'nB'])).items():frame.at[i,key]=value
    frame['original_origin']=[r['original_origin'] for r in rows()]
    frame['observation_cohort']=np.where(frame.origin.eq('reused'),'inherited unified pilot','new expansion')
    for key in ['kA','kB','k']:
        frame['raw_'+key]=frame[key]
        frame.loc[~frame.validated,key]=np.nan
    if (OUT/'execution_manifest_alpha05.json').exists():
        frame.to_csv(OUT/'all_alpha_results.csv',index=False)
        frame[frame.alpha.eq(1.)].to_csv(OUT/'alpha1_pilot_evidence.csv',index=False)
        frame=frame[frame.alpha.eq(.5)].copy()
        frame['pair_cohort']=np.where(frame.origin.eq('reused'),'original pilot pair','newly selected pair')
    frame.to_csv(OUT/'results.csv',index=False)
    counts=dict(configurations=len(frame),inherited=int(frame.origin.eq('reused').sum()),new_executed=int(frame.origin.eq('new').sum()),
        validated=int(frame.validated.sum()),resource_limited=int(frame.resource_limited.sum()),unexecuted=int(frame.origin.eq('unexecuted').sum()))
    atomic(OUT/'counts.json',counts);return frame

class LiveProgress:
    """Five-second observation outside workers and all optimization loops."""
    def __init__(self,frame_loader=None,output=None,total_new=24,batch_started_unix=None):
        import threading
        self.loader=frame_loader or save;self.output=output or OUT;self.total_new=total_new
        self.start=time.monotonic();self.data={};self.graph={};self.frame=self.loader();self.lock=threading.RLock()
        self.batch_started_unix=batch_started_unix or read(OUT/'execution_manifest_alpha05.json')['unix']
        self.stop_event=threading.Event();self.thread=threading.Thread(target=self.loop,daemon=True)
    def __enter__(self):self.thread.start();return self
    def __exit__(self,*args):self.stop_event.set();self.thread.join();self.publish()
    def current(self,r,phase):
        with self.lock:
            self.graph={};self.data=dict(pair=r['pair'],method=r['method'],nA=r['nA'],nB=r['nB'],configuration_id=r['configuration_id'],phase=phase,phase_started=time.monotonic(),phase_budget_seconds=None)
            if 'repetition' in r:self.data['repetition']=r['repetition']
    def observe(self,event):
        with self.lock:
            self.data.update(event)
            self.data['phase_started']=time.monotonic()-event['phase_elapsed_seconds']
            p=event.get('worker_progress',{})
            keys=['states_discovered','states_dequeued','completed_states_lower_bound','product_transitions','queue_size','chain_graphs_completed']
            if any(k in p for k in keys):self.graph={k:p.get(k) for k in keys}|dict(counter_timestamp=p.get('utc'))
    def refresh(self):
        with self.lock:self.frame=self.loader()
    def finish(self,state,reason=None):
        with self.lock:self.data.update(phase=state,pause_reason=reason,phase_started=time.monotonic(),phase_budget_seconds=None)
    def loop(self):
        while not self.stop_event.wait(5):self.publish()
    def publish(self):
        with self.lock:
            f=self.frame;new=f[f.origin.ne('reused')];finished=int(new.origin.eq('new').sum());incomplete=int(f.status.eq('validation_incomplete').sum())
            counts=dict(new_slots_finished=finished,new_slots_total=self.total_new,new_slots_finished_percent=100*finished/self.total_new,
                reused_observations=int(f.origin.eq('reused').sum()),validated=int(f.validated.sum()),
                resource_limited=int((f.resource_limited&~f.status.eq('validation_incomplete')).sum()),validation_incomplete=incomplete,
                failed=int(f.status.isin(['software_error','validation_failed','attempt_incomplete']).sum()),pending=int(f.origin.eq('unexecuted').sum()),
                solver_completed=int(f.get('raw_status',pd.Series(dtype=str)).fillna('').str.startswith('optimal').sum()))
            state={**counts,**{k:v for k,v in self.data.items() if k not in ['worker_progress','phase_started']},**stamp(),
                'phase_elapsed_seconds':time.monotonic()-self.data.get('phase_started',time.monotonic()),
                'current_session_elapsed_seconds':time.monotonic()-self.start,
                'total_batch_elapsed_seconds':time.time()-self.batch_started_unix,
                'batch_elapsed_definition':'Wall time since alpha0.5 manifest freeze, including pauses and recovery; current session separately timed',
                'graph_counters':self.graph,'counter_note':'states_dequeued counts heap pops; completed_states_lower_bound counts prior expanded nonterminal states; not search completion'}
            state['budget_used_percent']={k:100*self.graph[k]/limit for k,limit in [('states_discovered',100000),('product_transitions',20000000)] if self.graph.get(k) is not None}
            if self.total_new==120:state.update(finished_repetitions=finished,planned_repetitions=self.total_new,finished_repetitions_percent=100*finished/self.total_new)
            atomic(self.output/'progress.json',state)
            with (self.output/'progress_history.jsonl').open('a',encoding='utf-8') as log:log.write(json.dumps(state)+'\n')
            print(f"PROGRESS {finished}/{self.total_new} ({100*finished/self.total_new:.1f}%) finished; reused={counts['reused_observations']} validated={counts['validated']} resource={counts['resource_limited']} validation_incomplete={incomplete} failed={counts['failed']} pending={counts['pending']} | {state.get('pair','')} {state.get('method','')} rep={state.get('repetition','-')} {state.get('phase','')} {state['phase_elapsed_seconds']:.1f}s/{state.get('phase_budget_seconds')}s | counters={self.graph}",flush=True)

def primary_prepare():
    import shutil
    runner.cohort_active_check()
    if not (OUT/'execution_manifest_alpha05.json').exists():
        archive=OUT/'presentation_before_alpha05';archive.mkdir(exist_ok=False)
        for name in ['results.csv','counts.json','report.md','status.md','graph_bounds.csv','graph_state_summary.csv','matched_comparisons.csv','completion.csv','final_audit.json','artifact_hashes.json']:
            if (OUT/name).exists():shutil.copy2(OUT/name,archive/name)
        shutil.copytree(OUT/'figures',archive/'figures')
        selected=[r for r in rows() if r['alpha']==.5];assert len(selected)==40
        atomic(OUT/'execution_manifest_alpha05.json',dict(rows=selected,original_manifest_sha256=sha(OUT/'configuration_manifest.json'),
            alpha_selection='alpha=0.5 chosen after pilot evidence; exploratory comparison, not independent confirmatory selection',
            policy='Preserve 16 inherited observations, recover slot33 validation only, execute untouched new alpha0.5 slots once. No alpha1 launches.',**stamp()))
        pd.DataFrame([{k:v for k,v in r.items() if k!='configuration'} for r in selected]).to_csv(OUT/'execution_manifest_alpha05.csv',index=False)
    save()

def recover(slot_number=33):
    from supervisor import attempt
    runner.cohort_active_check();primary_prepare();base=OUT/('validation_recovery' if slot_number==33 else f'validation_recovery_slot{slot_number}');base.mkdir(exist_ok=True)
    existing=list((base/'attempts').glob('*/configuration.json'))
    r=next(r for r in rows() if r['plan_order']==slot_number);original=next((OUT/'slots'/f'{slot_number:02d}'/'attempts').glob('*/raw_solver.json'))
    if existing:
        assert len(existing)==1;d=existing[0].parent
        if not (d/'validation.json').exists():raise RuntimeError('Recovery already started but not finalized; inspect, no retry')
    else:
        c={**read(original)['configuration'],'diagnostic':f'slot{slot_number}_validation_recovery','validation_replay_raw':str(original.relative_to(ROOT)),
            'validation_replay_sha256':sha(original)}
        with LiveProgress() as progress:
            progress.current(r,'validation_recovery_preparation')
            d=attempt(c,output_directory=base,monitor=progress.observe)
            progress.finish('recovery_finished')
    atomic(base/'index.json',dict(configuration_id=r['configuration_id'],attempt_directory=str(d.relative_to(ROOT)),original_raw_sha256=sha(original),optimization_reexecuted=False))
    save();print('VALIDATION RECOVERY',read(d/'validation.json'),flush=True)
    if read(d/'validation.json')['status']=='validation_failed':raise RuntimeError('Feasibility/validator failure; stop batch')

def primary_execute():
    from supervisor import attempt
    runner.cohort_active_check();primary_prepare();configure()
    recovery=OUT/'validation_recovery/index.json'
    assert recovery.exists(),'Recover saved validation first'
    assert read(ROOT/read(recovery)['attempt_directory']/'validation.json')['status']!='validation_failed'
    with LiveProgress() as progress:
        for r in read(OUT/'execution_manifest_alpha05.json')['rows']:
            if r['reuse_path']:continue
            slot=OUT/'slots'/f"{r['plan_order']:02d}"
            existing=list((slot/'attempts').glob('*/configuration.json'))
            if existing:
                outcome=save().set_index('configuration_id').loc[r['configuration_id']].to_dict()
                if outcome['status'] in ['software_error','validation_failed','attempt_incomplete']:
                    progress.finish('paused','Existing systemic error: '+r['configuration_id']);return
                continue
            if (slot/'launch.json').exists():raise RuntimeError('Incomplete launch marker; no retry')
            progress.current(r,'host_readiness')
            if not runner.cohort_readiness(slot):
                atomic(OUT/'primary_pause.json',dict(reason='host_readiness',configuration_id=r['configuration_id'],**stamp()))
                progress.finish('paused','Host readiness failed; see slot readiness record');return
            config=read(slot/'configuration.json');assert config==r['configuration'] and config['alpha']==.5
            atomic(slot/'launch.json',dict(configuration_sha256=sha(slot/'configuration.json'),automatic_retries=False,**stamp()))
            progress.current(r,'startup');d=attempt(config,output_directory=slot,monitor=progress.observe)
            atomic(slot/'attempt_index.json',dict(attempt_directory=str(d.relative_to(ROOT))))
            outcome=runner.cohort_result(r);progress.refresh()
            with progress.lock:
                progress.graph={k:outcome.get(v) for k,v in [('states_discovered','states'),('states_dequeued','states_dequeued'),('completed_states_lower_bound','completed_states_lower_bound'),('product_transitions','transitions'),('queue_size','queue_size'),('chain_graphs_completed','chain_graphs_completed')]}
                progress.graph['counter_timestamp']=read(d/'supervisor.json')['utc']
            print('END',r['plan_order'],outcome['status'],outcome.get('k'),flush=True)
            if outcome['status'] in ['software_error','validation_failed','attempt_incomplete'] or str(outcome.get('raw_status')).startswith('infeasible'):
                atomic(OUT/'primary_pause.json',dict(reason='systemic software/feasibility review',configuration_id=r['configuration_id'],status=outcome['status'],**stamp()))
                progress.finish('paused','Systemic software or feasibility failure');return
        progress.refresh();progress.finish('completed');atomic(OUT/'primary_completed.json',dict(**stamp(),counts=read(OUT/'counts.json')))

def execute():
    from supervisor import attempt
    configure();runner.cohort_active_check()
    for r in rows():
        if r['reuse_path']:continue
        slot=OUT/'slots'/f"{r['plan_order']:02d}"
        existing=list((slot/'attempts').glob('*/configuration.json'))
        if existing:
            status=runner.cohort_result(r)['status']
            if status in ['software_error','validation_failed','attempt_incomplete']:break
            continue
        if (slot/'launch.json').exists():raise RuntimeError('Incomplete launch marker; no retry')
        if not runner.cohort_readiness(slot):
            atomic(OUT/'pause.json',dict(configuration_id=r['configuration_id'],reason='host_readiness',**stamp()));break
        config=read(slot/'configuration.json');assert config==r['configuration']
        assert sha(ROOT/config['input_file'])==config['input_sha256']
        for n,h in config['algorithms'].items():assert sha(ROOT/n)==h
        atomic(slot/'launch.json',dict(configuration_sha256=sha(slot/'configuration.json'),automatic_retries=False,**stamp()))
        print('START',r['plan_order'],r['configuration_id'],flush=True)
        d=attempt(config,output_directory=slot);atomic(slot/'attempt_index.json',dict(attempt_directory=str(d.relative_to(ROOT))))
        result=runner.cohort_result(r);save();print('END',r['plan_order'],result['status'],result['k'],flush=True)
        if result['status'] in ['software_error','validation_failed','attempt_incomplete'] or str(result.get('raw_status')).startswith('infeasible'):
            atomic(OUT/'pause.json',dict(reason='systemic review required',configuration_id=r['configuration_id'],status=result['status'],**stamp()));break
    save();print('CHECKPOINT',read(OUT/'counts.json'),flush=True)

THEORY='''# Graph counting and numerical bounds

## Sources and assumptions
Gozoltzani thesis: Chapter 3, Observation 4, printed p.18 (PDF p.26);
Chapter 5, section 5.1, printed p.26 (PDF p.34). Manuscript:
sections 3.1, 4.1 and 4.4, PDF pp.7-8, 11-12, 15-16.
The experimental chapters are provisional and are not evidence for these results.
Both constructions use original output-vertex indices and matching positions on
the extended input curves. For a curve with n original vertices, each of n
spheres intersects each of n-1 nondegenerate straight edges in at most two
locations. Thus U(n)=2n(n-1) is a finite upper bound on added locations and
n*=n+a <= n+U(n). Assumption: no edge is contained in a sphere boundary.
The implementation skips degenerate edges and preserves their original endpoints,
counts distinct interior curve parameters (not distinct spatial coordinates),
and merges parameters within 1e-10. Endpoint intersections are not added.
The numerical tolerances do not create more than two candidates per edge/sphere.
No packedness assumption is made or tested here.

## Mapping to the implementation
For one curve, dense indices (i,p) cover S_A=n(n+aA) locations: owner i in
the extended curve and selectable original vertex p. CPS-3F has aA=0 and
S_A=n^2. These dense slots are NOT all valid single-curve graph vertices.
The code removes vertices not on an anchored source-to-terminal path; the
reported active count H_A is the remaining number. Source and terminal are
included. Every active vertex has one zero-cost self/wait transition.
E_A counts outgoing entries after this pruning, including these waits.
Auxiliary locations cannot be selected as simplified output vertices.

The finite input-only configuration-state bound is
B_input=n*m*(n+2n(n-1))*(m+2m(m-1)) for CPS-2F;
for CPS-3F it is n^2*m^2 (one bound, no duplicate auxiliary-based bound).
When completed single-curve graphs supply auxiliary counts,
B_extended=S_A*S_B is a tighter finite combinatorial bound. When active
counts are saved, B_active=H_A*H_B is tighter still. None enumerates all
coupling-valid states: the output coupling condition can only remove states.
The root and terminal require no additional vertices beyond these products.
The theory's O(n^3*m^3) counts potential CPS-2F configurations, not edges,
auxiliary locations or DP labels. It is displayed separately, never evaluated
with a hidden constant of one. CPS-3F's finite n^2*m^2 construction bound has
no auxiliary locations; auxiliary utilization 0/0 is not defined.

Discovered states are len(values), including the source. Dequeued states
count heap pops, including a terminal pop. completed_states_lower_bound counts
previously expanded nonterminal states; the last interrupted expansion may be
partial. Queue size is pending heap entries, not a progress percentage.
The product graph is implicit; all valid or reachable states are not enumerated
in advance. No claim that discovered states equal all valid states is made.

Candidate transition checks are counted only after excluding simultaneous waits
and testing destination coupling. A budget-triggering check is included even
when not relaxed. They are not reported as stored unique graph edges. A finite
bound on these checks, where both single graphs completed, is
E_A*E_B-H_A*H_B: each pair of outgoing single-graph entries defines one product
source and destination, the H_A*H_B simultaneous waits are excluded, and each
source is popped at most once because the queued set is never cleared.
Each single graph has unique (source,destination) entries by the construction's
index enumeration. Coupling and reachability further reduce checked candidates.
No inner-loop instrumentation was added to count stored product edges.
Single-graph edges including waits obey H*(H+1)/2 because every non-wait edge
strictly advances the dense topological index and occurs at most once.

DP labels X_v[r] are distinct from states: r counts advances in A's selected
index; at most n labels per state, each storing a best B count. Labels are not
instrumented. The objective reported here counts vertices including the first
one: max(r+1,z+1), not max(r,z). Unequal lengths and repeated coupling indices
are allowed. Feasibility validation alone does not prove optimality.

All percentages mean 100*observed/finite_bound, NOT search progress or a
prediction of remaining work. Interrupted counts are observations before
termination. Historical missing counters stay missing; no runs were repeated
to fill them. A missing percentage is explained by missing evidence, no graph
execution, or a zero/non-applicable denominator, never replaced by zero.
'''

def graph_tables(frame):
    records=[]
    for r in frame.to_dict('records'):
        if not r['method'].startswith('CPS'):continue
        d=ROOT/r['attempt_path'] if r['attempt_path'] else None
        raw=read(d/'raw_solver.json') if d and (d/'raw_solver.json').exists() else {}
        events=[json.loads(x) for x in (d/'phases.jsonl').read_text().splitlines()] if d and (d/'phases.jsonl').exists() else []
        obs={}
        for e in events:
            if 'solver_subphase' in e:obs.update(e)
        obs.update(raw.get('observation',{}));payload=raw.get('raw_result',{})
        n,m=int(r['nA']),int(r['nB']);joint2=r['method']=='CPS-2F'
        base={k:r[k] for k in ['configuration_id','plan_order','domain','pair','alpha','method','nA','nB','origin','attempt_path','status']}
        base.update(search_completed=raw.get('solver_status')=='optimal_configuration_graph',construction_completed=r.get('chain_graphs_completed')==2,
            count_scope='completed search' if raw.get('solver_status')=='optimal_configuration_graph' else 'observed before termination' if r['origin']!='unexecuted' else 'not executed')
        def add(quantity,count,bound,formula,assumptions,source='Thesis Observation 4 p.18 and section 5.1 p.26; manuscript section 4.4 pp.15-16'):
            available=count is not None and pd.notna(count)
            records.append({**base,'quantity':quantity,'observed':count if available else None,'finite_upper_bound':bound,
                'percent_of_bound':100*float(count)/bound if available and bound and bound>0 else None,
                'formula':formula,'assumptions':assumptions,'source':source,
                'asymptotic_states':'O(n^3 m^3)' if joint2 else 'O(n^2 m^2)',
                'missing_reason':'' if available else 'Counter unavailable or phase not reached; historical runs not rerun',
                'bound_missing_reason':'Missing completed single-graph counters' if bound is None else 'Zero/non-applicable bound' if bound==0 else ''})
        graphs=[]
        for side,k in [('A',n),('B',m)]:
            g=obs.get('single_graph_'+side,{})
            completed=r.get('chain_graphs_completed')
            aux=g.get('auxiliary',payload.get('auxiliary_'+side,0 if not joint2 and pd.notna(completed) and completed>=(1 if side=='A' else 2) else None))
            if not g and side=='B' and r.get('chain_graphs_completed')==2:
                aux=obs.get('last_completed_chain_auxiliary',aux)
            dense=g.get('dense_states',k*(k+aux) if aux is not None else None)
            active=g.get('active_states');edges=g.get('edges_including_waits')
            if edges is None and side=='B' and r.get('chain_graphs_completed')==2:edges=obs.get('last_completed_chain_edges')
            graphs.append(dict(aux=aux,dense=dense,active=active,edges=edges))
            add('auxiliary_'+side,aux,2*k*(k-1) if joint2 else 0,'2 k(k-1)' if joint2 else '0 (no auxiliary locations)', 'Interior curve-parameter intersections; CPS-3F percentage not applicable')
            add('single_dense_'+side,dense,k*(k+2*k*(k-1)) if joint2 else k*k,'k*(k+2k(k-1))' if joint2 else 'k^2','Dense indexing slots, not valid vertices')
            add('single_active_'+side,active,dense,'k*(k+a)','Anchored source/terminal included; single-graph construction complete')
            add('single_edges_with_waits_'+side,edges,active*(active+1)//2 if active is not None else None,'H*(H+1)/2','Unique forward transitions plus exactly H waits; completed graph','Finite derivation from chain_graph() enumeration and pruning')
        bound=n*m*(n+2*n*(n-1))*(m+2*m*(m-1)) if joint2 else n*n*m*m
        add('configuration_states_input_bound',r.get('states'),bound,'n*m*(n+2n(n-1))*(m+2m(m-1))' if joint2 else 'n^2*m^2','All index combinations; includes source/terminal; no hidden asymptotic constant')
        if joint2:
            ext=graphs[0]['dense']*graphs[1]['dense'] if all(g['dense'] is not None for g in graphs) else None
            add('configuration_states_extended_bound',r.get('states'),ext,'n*m*(n+aA)*(m+aB)','Actual extended-curve sizes; all index combinations')
        activeprod=graphs[0]['active']*graphs[1]['active'] if all(g['active'] is not None for g in graphs) else None
        add('configuration_states_active_bound',r.get('states'),activeprod,'H_A*H_B','Completed single graphs; includes coupling-invalid combinations')
        checks=graphs[0]['edges']*graphs[1]['edges']-activeprod if activeprod is not None and all(g['edges'] is not None for g in graphs) else None
        add('candidate_transition_checks',r.get('transitions'),checks,'E_A*E_B-H_A*H_B','No duplicate single-graph entries; each product source popped once; excludes simultaneous waits','Finite derivation from solve_cps() and chain_graph(); see graph_bounds.md')
        add('unique_configuration_edges',None,activeprod*(activeprod-1)//2 if activeprod is not None else None,'N*(N-1)/2, N=H_A*H_B','Product graph not materialized; no observed unique-edge count; strict forward topological order','Finite DAG counting derivation; not a transition-check measurement')
        add('dynamic_programming_labels',None,n*activeprod if activeprod is not None else None,'n*H_A*H_B','At most n distinct A-hop counts per state; labels not instrumented','Thesis section5.1; manuscript section4.4; solve_cps values dictionaries')
    table=pd.DataFrame(records);table.to_csv(OUT/'graph_bounds.csv',index=False)
    state=table[table.quantity=='configuration_states_input_bound'].copy()
    state[['configuration_id','domain','pair','alpha','method','status','observed','finite_upper_bound','percent_of_bound',
           'search_completed','construction_completed','count_scope','formula','source']].to_csv(OUT/'graph_state_summary.csv',index=False)
    (OUT/'graph_bounds.md').write_text(THEORY,encoding='utf-8');return table

def report():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import nbformat
    f=save();primary=(OUT/'execution_manifest_alpha05.json').exists();g=graph_tables(f);figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','text.color':'black','axes.labelcolor':'black',
        'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.bbox':'tight'})
    captions=[]
    def finish(fig,name,caption):
        fig.savefig(figdir/(name+'.png'),dpi=180);fig.savefig(figdir/(name+'.pdf'));plt.close(fig);captions.append((name,caption))
    # Denominators include every frozen configuration, including unexecuted ones.
    completion=f.groupby(['domain','method','termination_reason'],sort=False).size().rename('count').reset_index()
    completion['denominator']=completion.apply(lambda r:int(((f.domain==r.domain)&(f.method==r.method)).sum()),axis=1)
    completion.to_csv(OUT/'completion.csv',index=False)
    statuses=['validated']+sorted(set(f.termination_reason)-{'validated'})
    statuscolors={'validated':'#4daf4a','unexecuted':'#eeeeee','transition_budget_limit':'#e66101','state_budget_limit':'#984ea3',
        'single_chain_edge_budget_limit':'#377eb8','solver_timeout':'#b2182b','memory_limit':'#a6761d','startup_failure':'#555555'}
    fig,axes=plt.subplots(1,2,figsize=(12,4),sharey=True)
    for ax,domain in zip(axes,['protein','hurricane']):
        sub=f[f.domain==domain];bottom=np.zeros(4)
        for status in statuses:
            counts=np.array([int(((sub.method==m)&(sub.termination_reason==status)).sum()) for m in METHODS])
            ax.bar(np.arange(4),counts,bottom=bottom,label=status.replace('_',' '),color=statuscolors.get(status,'#666666'),edgecolor='white');bottom+=counts
        den=len(sub)//4
        ax.set_xticks(range(4),[m.replace('Independent ','Ind. ') for m in METHODS],rotation=25,ha='right');ax.set_title(domain.capitalize()+f' ({den} configurations/method)');ax.set_ylabel('Configurations');ax.set_ylim(0,den+.8)
        for x,m in enumerate(METHODS):ax.text(x,den+.15,f'{int(((sub.method==m)&sub.validated).sum())}/{den}',ha='center',fontsize=9)
    axes[1].legend(fontsize=8,bbox_to_anchor=(1.02,1));fig.tight_layout()
    finish(fig,'completion','Validated / intended configurations by domain and method. Colored stacks retain each distinct termination reason; white-grey denotes not executed.')
    cases=f[['pair','alpha','domain']].drop_duplicates().reset_index(drop=True);cases['x']=range(len(cases))
    plotted=f.merge(cases,on=['pair','alpha','domain'])
    fig,ax=plt.subplots(figsize=(15,5))
    for i,meth in enumerate(METHODS):
        s=plotted[(plotted.method==meth)&plotted.validated];x=s.x+(i-1.5)*.18
        ax.scatter(x,s.k,color=COLORS[meth],s=38,label=meth)
        if meth.startswith('Independent'):
            bad=~s.coupling_feasible.fillna(False).astype(bool);ax.scatter(x[bad],s.loc[bad,'k'],marker='x',s=95,color='black',linewidths=1.3)
    labels=[f'{r.pair}\nalpha={r.alpha:g}' for r in cases.itertuples()]
    if primary:
        ax.axvspan(3.5,9.5,color='#eeeeee',alpha=.4,zorder=-1)
        ax.text(.72,.62,'Shaded: newly selected pairs',ha='center',transform=ax.transAxes,color='#555555',fontsize=9)
    else:
        ax.axvspan(7.5,19.5,color='#eeeeee',alpha=.4,zorder=-1)
    handles,legendlabels=ax.get_legend_handles_labels()
    handles.append(Line2D([],[],marker='x',color='black',linestyle='',label='Independent coupling exceeds budget'))
    ax.set_xticks(cases.x,labels,rotation=90);ax.set_ylabel(r'$k=\max(|A\prime|,|B\prime|)$');ax.legend(handles=handles,ncol=3,loc='upper center',bbox_to_anchor=(.5,1.2),fontsize=9);ax.grid(axis='y',alpha=.2);fig.tight_layout()
    finish(fig,'solution_size',f'Validated outputs only on all {len(cases)} matched pair/alpha slots. Black crosses mark independent outputs exceeding delta3 (not baseline failures). Missing outputs are omitted, never zero; see completion table. Shading identifies newly selected pairs; alpha was selected using pilot evidence.')
    fig,ax=plt.subplots(figsize=(8,5))
    for r in f.itertuples():
        if pd.isna(r.solver_seconds):continue
        marker='o' if r.validated else 's' if str(r.raw_status).startswith('optimal') else 'X'
        ax.scatter(r.nA+r.nB,max(r.solver_seconds,1e-5),marker=marker,s=55,
            edgecolors=COLORS[r.method],facecolors='none' if r.origin=='reused' else COLORS[r.method],alpha=.8)
    ax.set_yscale('log');ax.set_xlabel('Total input vertices n + m');ax.set_ylabel('Solver seconds (log scale)');ax.grid(alpha=.2)
    handles=[Line2D([],[],marker='o',color=c,linestyle='',label=m) for m,c in COLORS.items()]
    handles += [Line2D([],[],marker='X',color='black',linestyle='',label='Solver resource exit'),Line2D([],[],marker='s',color='black',linestyle='',label='Solver returned; validation incomplete'),Line2D([],[],marker='o',markerfacecolor='none',color='black',linestyle='',label='Inherited observation')]
    ax.legend(handles=handles,fontsize=8,loc='upper left',bbox_to_anchor=(1.02,1.0));fig.tight_layout()
    finish(fig,'solver_time','Circles: validated solutions; X: solver resource exits (time to termination); square: solver returned but validation incomplete (time to raw return, not validation termination). Hollow markers: inherited observations. One observation/configuration; added boundary instrumentation precludes strict timing equivalence.')
    states=g[g.quantity=='configuration_states_input_bound'].dropna(subset=['observed'])
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for r in states.itertuples():
        marker='o' if r.search_completed else 'X';color=COLORS[r.method]
        axes[0].scatter(r.finite_upper_bound,max(r.observed,1),marker=marker,color=color,s=45)
        axes[1].scatter(r.nA+r.nB,r.percent_of_bound,marker=marker,color=color,s=45)
    axes[0].set(xscale='log',yscale='log',xlabel='Finite input-only state bound',ylabel='Discovered configuration states')
    axes[1].set(yscale='log',xlabel='Total input vertices n + m',ylabel='100 x discovered / finite bound (%)')
    if len(states):
        lo=max(1,min(states.observed.min(),states.finite_upper_bound.min()));hi=states.finite_upper_bound.max()
        axes[0].plot([lo,hi],[lo,hi],color='grey',ls='--',lw=1)
    handles=[Line2D([],[],marker='o',color=COLORS[m],linestyle='',label=m) for m in ['CPS-2F','CPS-3F']]+[Line2D([],[],marker='X',color='black',linestyle='',label='Observed before termination')]
    axes[0].legend(handles=handles,fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();finish(fig,'configuration_states','Input-only finite bound: CPS-2F n*m*(n+2n(n-1))*(m+2m(m-1)); CPS-3F n^2*m^2. Includes all index combinations and endpoints, not just valid states. X indicates interrupted observations. Dashed line is equality, not a fitted scaling law.')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for ax,quantity,title in zip(axes,['configuration_states_extended_bound','configuration_states_active_bound'],['Actual extended-curve bound','Active single-graph product bound']):
        sub=g[g.quantity==quantity].dropna(subset=['percent_of_bound'])
        for r in sub.itertuples():ax.scatter(r.nA+r.nB,r.percent_of_bound,color=COLORS[r.method],marker='o' if r.search_completed else 'X',s=45)
        ax.set(title=title,xlabel='Total input vertices n + m',ylabel='100 x discovered / finite bound (%)',yscale='log');ax.grid(alpha=.2)
        denominator=int(f.method.eq('CPS-2F').sum()) if quantity=='configuration_states_extended_bound' else int(f.method.str.startswith('CPS').sum())
        if len(sub):ax.text(.53,.97,f'{len(sub)} / {denominator} applicable slots with evidence',transform=ax.transAxes,ha='center',va='top',fontsize=8)
        else:
            ax.set_axis_off();ax.text(.5,.55,'Active-state counters unavailable\nfor inherited observations.\n\nNo new CPS graph was executed.\n0 / 40 slots with evidence.',transform=ax.transAxes,ha='center',va='center',fontsize=12)
    axes[0].legend(handles=handles,fontsize=8,loc='upper right',bbox_to_anchor=(.99,.85));fig.tight_layout()
    finish(fig,'state_utilization','Left: S_A*S_B from actual extended lengths (CPS-2F only; CPS-3F repeats its input bound and is omitted). Right: H_A*H_B from saved active counts. Missing historical counters stay missing. Percentages measure occupancy of bounds, never progress.')
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    for ax,prefix,title in zip(axes,['auxiliary_','single_active_','single_edges_with_waits_'],['Auxiliary matching locations','Active single-curve states','Single-curve edges, including waits']):
        sub=g[g.quantity.str.startswith(prefix)].dropna(subset=['observed'])
        for r in sub.itertuples():
            k=r.nA if r.quantity.endswith('_A') else r.nB
            ax.scatter(k,r.observed,color=COLORS[r.method],marker='o' if r.construction_completed else 'X',s=20,alpha=.65)
        ax.set(xlabel='Vertices in the corresponding input curve',ylabel=title);ax.set_yscale('symlog',linthresh=1);ax.grid(alpha=.2)
        if not len(sub):
            ax.set_axis_off();ax.text(.5,.5,'Active single-curve states\nnot recorded in inherited runs.\nNo new CPS graph executed.',transform=ax.transAxes,ha='center',va='center',fontsize=10)
    axes[0].legend(handles=handles[:2],fontsize=8);fig.tight_layout()
    finish(fig,'single_curve_sizes','Each point represents curve A or B for a CPS configuration; thresholds are not independent repetitions. No auxiliary locations in CPS-3F. Graph edges include self/wait entries. Only saved counters are shown; no implied linear growth or packedness claim.')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    edges=g[g.quantity.str.startswith('single_edges_with_waits_')].dropna(subset=['observed','finite_upper_bound'])
    for r in edges.itertuples():
        k=r.nA if r.quantity.endswith('_A') else r.nB
        axes[0].scatter(r.finite_upper_bound,r.observed,color=COLORS[r.method],marker='o' if r.search_completed else 'X',s=35,alpha=.8)
        axes[1].scatter(k,r.percent_of_bound,color=COLORS[r.method],marker='o' if r.search_completed else 'X',s=35,alpha=.8)
    axes[0].set(xscale='log',yscale='log',xlabel='Finite single-graph edge bound H(H+1)/2',ylabel='Observed edges, including H waits')
    axes[1].set(xlabel='Input curve vertices',ylabel='100 x observed edges / finite edge bound (%)')
    axes[0].legend(handles=handles[:3],fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();finish(fig,'single_curve_edge_bounds','Completed single-curve graphs only. H is the active vertex count; unique forward edges plus one wait at each active vertex obey H(H+1)/2. Circles: completed product search; X: subsequently interrupted product search. These are single-graph edges, not product-transition checks. Historical missing H counts remain missing.')
    fig,ax=plt.subplots(figsize=(11,2.8));ax.axis('off')
    texts=['Input curves\nOriginal vertices A, B','Single-curve graphs\nOriginal output choices +\nauxiliary matching positions','Configuration search\nContinuous/discrete fidelity\n+ discrete output coupling','Simplified curves\nOriginal vertices only\nUnequal lengths allowed']
    for x,txt in zip([.1,.36,.64,.9],texts):
        ax.text(x,.6,txt,ha='center',va='center',transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.7',facecolor='#f4f4f4',edgecolor='#555555'),fontsize=9)
    for x,y in [(.2,.245),(.475,.52),(.755,.795)]:ax.annotate('',xy=(y,.6),xytext=(x,.6),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=1.5))
    ax.text(.5,.08,'Auxiliary locations encode matching positions; they are never selectable output vertices. CPS-3F adds none.',ha='center',transform=ax.transAxes,fontsize=10)
    finish(fig,'algorithm_diagram','Anchored global-fidelity configuration construction and coupled output. Independent baselines solve the single-curve problems without enforcing output coupling.')
    # Requested fixed illustration, never selected based on new outcomes.
    example=f[(f.pair=='AL031854__AL011991')&(f.alpha==.5)&f.method.isin(['Independent continuous','CPS-2F'])]
    fig,axes=plt.subplots(1,2,figsize=(11,5));exrecords=[]
    for ax,meth in zip(axes,['Independent continuous','CPS-2F']):
        r=example[example.method==meth].iloc[0];raw=read(ROOT/r.attempt_path/'raw_solver.json')['raw_result'];data=np.load(ROOT/r.input_file);A,B=data['A'],data['B'];val=read(ROOT/r.attempt_path/'validation.json')
        a,b=A[raw['indices_A']],B[raw['indices_B']]
        for P,Q,col,name in [(A,a,'#377eb8','A'),(B,b,'#e66101','B')]:
            ax.plot(P[:,0],P[:,1],color=col,alpha=.25,lw=1,label=name+' input');ax.plot(Q[:,0],Q[:,1],'-o',color=col,lw=2,ms=4,label=name+' simplified')
        for i,j in val['coupling']:ax.plot([a[i,0],b[j,0]],[a[i,1],b[j,1]],color='grey',ls=':',alpha=.3,lw=.7)
        ax.set_title(f'{meth}: k={int(r.k)}\ncoupling={r.coupling_distance:.3f} km; budget={r.delta3:.3f} km',fontsize=10)
        ax.set_aspect('equal',adjustable='datalim');ax.set_xlabel('Basin x (km)');ax.set_ylabel('Basin y (km)');ax.legend(fontsize=8)
        exrecords.append(dict(method=meth,k=int(r.k),coupling_distance=float(r.coupling_distance),delta3=float(r.delta3),delta1=float(r.delta1),delta2=float(r.delta2),attempt_path=r.attempt_path))
    fig.tight_layout();finish(fig,'atlantic_coupling_example','Pre-specified illustration AL031854 / AL011991, alpha=0.5, full supplied 25-point tracks. Faint lines: inputs; bold lines: simplified curves; dotted lines: saved discrete coupling. Equal k does not imply the same coupling feasibility. Coordinates do not imply temporal synchronization.')
    atomic(OUT/'atlantic_example.json',exrecords)
    fig,ax=plt.subplots(figsize=(8,4))
    for i,meth in enumerate(METHODS):
        sub=f[(f.method==meth)&f.peak_rss_mib.notna()]
        for r in sub.itertuples():ax.scatter(i,r.peak_rss_mib,color=COLORS[meth],marker='o' if r.validated else 'X',alpha=.65)
    ax.axhline(512,color='black',ls='--',lw=1);ax.set_xticks(range(4),METHODS,rotation=15);ax.set_ylabel('Sampled process-tree peak RSS (MiB)');fig.tight_layout()
    finish(fig,'memory','Sampled whole-worker-tree RSS, including preparation and warm-up; not exact algorithm memory and not guaranteed peak capture. Dashed line: 512 MiB cap. Circles validated, X unsuccessful.')
    comparisons=[]
    for r in f[f.method.str.startswith('CPS')].to_dict('records'):
        bm='Independent continuous' if r['method']=='CPS-2F' else 'Independent discrete'
        baseline=f[(f.pair==r['pair'])&(f.alpha==r['alpha'])&(f.method==bm)].iloc[0]
        matched=bool(r['validated'] and baseline.validated)
        comparisons.append(dict(pair=r['pair'],domain=r['domain'],alpha=r['alpha'],method=r['method'],baseline=bm,
            delta1=r['delta1'],delta2=r['delta2'],delta3=r['delta3'],joint_status=r['status'],baseline_status=baseline.status,
            matched_validated=matched,k_joint=r['k'] if r['validated'] else None,k_baseline=baseline.k if baseline.validated else None,
            vertex_cost=r['k']-baseline.k if matched else None,baseline_coupling_feasible=baseline.coupling_feasible,
            joint_attempt=r['attempt_path'],baseline_attempt=baseline.attempt_path))
    comp=pd.DataFrame(comparisons);comp.to_csv(OUT/'matched_comparisons.csv',index=False)
    (OUT/'figure_captions.md').write_text('\n\n'.join(f'## {n}\n\n{c}' for n,c in captions),encoding='utf-8')
    counts=read(OUT/'counts.json');summary=[]
    runtime=[]
    for method in METHODS:
        s=f[(f.method==method)&f.validated]
        runtime.append(dict(method=method,validated_count=len(s),solver_median_seconds=s.solver_seconds.median(),
            solver_min_seconds=s.solver_seconds.min(),solver_max_seconds=s.solver_seconds.max(),
            sampled_peak_min_mib=s.peak_rss_mib.min(),sampled_peak_max_mib=s.peak_rss_mib.max(),
            scope='Validated single observations across distinct configurations, including disclosed reused observations; not timing repetitions'))
    pd.DataFrame(runtime).to_csv(OUT/'runtime_summary.csv',index=False)
    for method in METHODS:
        s=f[f.method==method];summary.append(f'| {method} | {int(s.validated.sum())}/{len(s)} | {int(s.resource_limited.sum())} | {int(s.origin.eq("unexecuted").sum())} |')
    quality=[]
    for method in ['CPS-2F','CPS-3F']:
        s=comp[(comp.method==method)&comp.matched_validated]
        quality.append(f'- {method}: {len(s)}/{len(cases)} validated matched comparisons; '+(f'vertex cost range {s.vertex_cost.min():g} to {s.vertex_cost.max():g}; {int(s.vertex_cost.eq(0).sum())} equal-size outcomes.' if len(s) else 'no paired quality claim available.'))
    table='| Method | Validated / planned | Resource exits | Unexecuted |\n|---|---:|---:|---:|\n'+'\n'.join(summary)
    report_text=f'''# Expanded bounded CPS comparison

## Scope and completion
Ten fixed pairs, two alpha values, four methods: 80 configurations.
Inherited observations: {counts['inherited']}; new attempts: {counts['new_executed']};
validated: {counts['validated']}; resource-limited: {counts['resource_limited']}; unexecuted: {counts['unexecuted']}.
No automatic retries, budget increases, new thresholds or timing repetitions.
Inherited observations retain their original attempt paths and are not new measurements.

{table}

## Inputs and protocol
[Input manifest](input_manifest.csv), [exact configurations](configuration_manifest.json),
[selection rule](selection_rule.md), [candidate ranking evidence](candidate_selection.json).
Three new pairs per domain were selected before solving, using input geometry and
reuse penalties only. Protein sizes are all 22/22 after supplied w=16 sampling of
already prepared 325-point curves; these are not full-resolution backbones.
The common 1o7j.a anchor is unavoidable in the supplied aligned pair list. Structure
reuse means pairs are not fully independent. Hurricanes preserve complete supplied
tracks and basin coordinates, including unequal lengths. No temporal synchronization
is inferred from coordinates. Gaps/provenance remain those of the supplied datasets.

alpha is 0.5 or 1; delta1=alpha*mean edge length of A, delta2 analogously for B;
delta3 is the unrounded discrete Frechet distance of saved inputs. Identical thresholds,
input arrays, tolerances and anchored endpoints apply to all methods in a configuration.
All methods minimize max(kA,kB). Independent continuous/discrete methods impose their
respective GLOBAL fidelity constraints separately and do NOT enforce delta3. CPS-2F
has continuous fidelity and discrete output coupling; CPS-3F has discrete fidelity
and discrete output coupling. Original vertices alone are output candidates; unequal
output lengths and repeated coupling indices are allowed. Certificates are disabled.

Shared external solver allowance 60s and sampled tree RSS 512MiB. CPS soft time 55s,
100000 states, 20000000 transition checks and 800000 main single-graph edges.
Independent discrete and CPS-3F independent preprocessing retain their 2000000-edge
defaults. Independent continuous has no corresponding graph-work cap. Internal soft
time is checked after expanded states, not continuously; external monitoring remains
the backstop. Separate budgets: startup60, loading15, warmup120, checkpoint5,
validation30, optional30, reporting10 seconds. Warmup is unchanged and output-local
caches are fresh per new slot. These practical bounds do not guarantee completion.

## Observed compression and coupling
{chr(10).join(quality)}
Only paired validated outputs support these comparisons; failures stay in the
denominator. Full retained fractions, selected sizes and coupling values are in
[results.csv](results.csv); [matched comparisons](matched_comparisons.csv) uses the
baseline with the same fidelity metric. Comparing CPS-2F with CPS-3F also changes
the fidelity constraint; size differences are not pure implementation speedups.
The pre-specified Atlantic alpha=0.5 illustration has independent continuous and
CPS-2F k=5: coupling 725.571465 km versus 426.184218 km with delta3=426.184218 km.
Thus joint feasibility can change even when optimum output size does not increase.
This illustration is not a general zero-cost claim or evidence of topology preservation.

## Computational cost and graph size
[Graph-bound table](graph_bounds.csv), [derivation and source references](graph_bounds.md).
[Compact state comparison](graph_state_summary.csv); [source files and hashes](source_references.json).
Observed discovered states are compared with explicit finite products of index sets,
not a big-O expression with constant one. Active states, dense slots, auxiliary
locations, edges including waits and candidate checks are distinct. Missing historical
counters stay missing. Interrupted counts are observed before termination; occupancy
percentages do not estimate search completion. DP labels and all coupling-valid states
were not separately counted. No packedness claim or general complexity law is inferred.
See figures/configuration_states and figures/state_utilization for conservative
input-only and tighter data-dependent bounds respectively.

## Validation and limitations
Raw outputs are atomically saved before external validation. CPS finish() already
performs internal fidelity/coupling checks; this work remains in solver time.
External validation checks original indices, order, endpoints, counts, applicable
fidelity, and discrete coupling; baseline coupling violations are informational.
Independent continuous also checks each returned simplification internally; this
check is included in solver time. Warmup compiles C-contiguous continuous-decision
inputs, whereas the frozen protein arrays are Fortran-contiguous. The stopped
slot's cache index contains only the C/C continuous-decision signature: validation
may compile an additional layout specialization. The existing warmup procedure
was preserved, not silently extended. Validation time is therefore not guaranteed
to exclude compilation. This and concurrent host activity are plausible contributors,
not a proven causal attribution for the timeout.
Validation establishes feasibility, not optimality; optimality is the solver's
claim supported by existing exhaustive small-instance checks, not independently
proved for each large input. Resource exits are not infeasibility. Original inputs
are feasible joint outputs because delta3 equals their discrete Frechet distance.
Runtime observations are single attempts with separate phase clocks; resource exits
are times to termination, never successful solution runtimes. RSS includes the worker
tree and startup/warmup, is sampled, and is not exact algorithm memory.

## Provenance and changes
Algorithms, worker, supervisor and warmup are unchanged. common.py adds only scalar
metadata at main single-graph completion; it does not instrument inner transition loops.
The added boundary observation may add small unquantified timing overhead relative to
inherited observations. expanded_comparison.py supplies selection, orchestration and
presentation using the existing runner/worker/validation. Hashes are in provenance.json
and each attempt configuration. instrumentation_checks.json records two small graph
regressions. This is not a repeated timing benchmark; additional thresholds or repeated
runs of a pair do not create independent input pairs. Manuscript experiments are not
used as verified evidence. No further experiments are authorized or launched here.

## Files
[Complete results](results.csv), [completion](completion.csv), [captions](figure_captions.md),
figures/*.pdf and *.png, raw per-attempt checkpoints/validation/phase/resource logs in slots/.
'''
    if (OUT/'pause.json').exists():report_text+='''\n## Execution pause and limits on conclusions

The expansion is **not complete**. One new independent-continuous solver returned
and saved its raw output (protein 1o7j.a / 1qd1.b, alpha=0.5), but external validation
exceeded 30 seconds. The raw objective is unverified and excluded from all quality
comparisons. The next readiness check recorded sustained PyCharm activity:
603.8, 621.4, 243.7, 153.2, 157.3, 148.1, 135.3, 90.9, 42.0, 37.5 percent of one CPU
core. The runner stopped before launching slot 34. Host activity at that check does
not prove the cause of the earlier validation timeout. No process was terminated
other than the attempt worker tree at its validation budget; user kernels remain.
No new CPS graph was run, so graph figures currently describe inherited evidence
only, and tighter active-graph bounds lack observed counters. The 29 verified
outputs are all inherited. No expanded-cohort compression, timing, or graph-scaling
conclusion is justified yet. Resume only the 47 untouched slots after host readiness;
do not retry the completed resource-limited observation automatically.
'''+json.dumps(read(OUT/'pause.json'),indent=2)+'\n'
    if primary:
        recovery=read(ROOT/read(OUT/'validation_recovery/index.json')['attempt_directory']/'validation.json')
        complete=int(f.get('raw_status',pd.Series(dtype=str)).fillna('').str.startswith('optimal').sum())
        termination=f.groupby('termination_reason').size().to_dict()
        graph_findings=[]
        for method in ['CPS-2F','CPS-3F']:
            s=g[(g.method==method)&(g.quantity=='configuration_states_active_bound')].dropna(subset=['percent_of_bound'])
            if len(s):graph_findings.append(f'- {method}: {len(s)}/10 configurations have the tighter active-state product bound; observed occupancy ranges from {s.percent_of_bound.min():.2f}% to {s.percent_of_bound.max():.2f}%, with {int((~s.search_completed).sum())} interrupted search(es) included explicitly.')
        stopped=f[f.resource_limited]
        stop_details='\n'.join(f"- {r.pair}, {r.method}: {r.status}; {int(r.transitions):,} transition checks, {int(r.states):,} discovered states, {r.solver_seconds:.3f} s to termination, {r.peak_rss_mib:.2f} MiB sampled peak RSS; no returned objective." for r in stopped.itertuples())
        report_text=f'''# Primary expanded comparison: alpha=0.5

## Coverage and provenance
Ten frozen pairs x four methods = **40 primary configurations**. alpha=0.5 was
chosen using earlier pilot evidence; this is an exploratory, evidence-informed
comparison, not a prespecified independent confirmatory study.

{table}

- Inherited observations: {counts['inherited']} (original solver/validation paths retained).
- New expansion slots executed: {counts['new_executed']}/24, including the previously executed slot 33.
- Solvers returned an optimal-status output: {complete}/40.
- Successfully validated: {counts['validated']}/40.
- Resource-limited observations without validated output: {counts['resource_limited']}.
- Unexecuted: {counts['unexecuted']}; termination categories: {json.dumps(termination)}.

All alpha=1 artifacts are preserved separately in alpha1_pilot_evidence.csv and
all_alpha_results.csv. No pending alpha=1 configuration was launched. Original
configuration_manifest.json is unchanged; execution_manifest_alpha05.json is the
revised execution plan. presentation_before_alpha05 preserves the prior presentation.
results.csv is the primary 40-row table, including failures and missing outputs.

## Validation recovery and targeted fix
Slot 33 returned in 0.0478647 s and saved 16 A indices and 20 B indices. The original
validation exceeded 30 s before optional measurements. Its frozen protein arrays
are Fortran-contiguous; warmup compiled the C/C layout of continuous_decision, and
the saved cache lacked the F/C specialization required by the original validator.
This exposed compilation inside validation. Contention was also observed later;
the original log does not identify the exact time spent compiling versus scheduled
CPU work, so a unique historical cause cannot be proved.

The targeted fix copies validation inputs to C-contiguous layout without changing
coordinate values, dtype, fidelity predicate, tolerances, endpoints, ordering or
membership checks. Five focused checks compare layouts and reject invalid fidelity,
endpoints, ordering and membership. Independent coupling is descriptive, not required.
Achieved continuous-distance searches remain in optional measurements outside validation.
No validation budget was increased. Validation-only replay passed in
{recovery.get('validation_seconds',float('nan')):.6f} s; kA={recovery.get('kA')},
kB={recovery.get('kB')}, k={recovery.get('k')}. No optimization was rerun for this output.
The original timeout remains unchanged, and validation_recovery contains the separate
replay, source checkpoint hash, worker code provenance, phase and resource logs.
The row retains original solver/warmup/total-process measurements; recovery timing
is a separate measurement, not a second solver repetition.

Slot34 subsequently saved an optimal-status independent-discrete output but hit
WinError5 while atomically replacing progress.json immediately before validation.
The orchestrator was paused after the current worker finished; no worker or unrelated
process was killed for this repair. A bounded 0.7-second retry of the I/O rename
(not a solver retry) handles transient Windows sharing denials within unchanged
phase budgets. The exact competing reader is unknown. A focused injected-denial
test passed. Slot34 validation-only recovery passed in approximately 0.0030 seconds;
its original solver output and error remain preserved. Outcome aggregation now
exposes worker software errors even when a raw optimal-status checkpoint exists.
See orchestration_repair.json, atomic_replacement_checks.json and validation_recovery_slot34/.

## Inputs, methods and resource policy
All ten input pairs and exact thresholds were frozen before expansion. Every protein
pair shares 1o7j and has 22 vertices per curve after existing w=16 sampling of supplied
prepared backbones. These pairs establish neither protein input-size scaling nor
independent biological replication. Newly selected hurricane tracks retain every
supplied point (including 25/82, 62/81 and 30/29 lengths), without alignment or cropping.
Shared basin coordinates are preserved; coordinates alone do not imply synchronized time.
Pair provenance and input hashes are in input_manifest.json. Original pilot pairs and
newly selected pairs are identified in results.csv and the size figure.

delta1=0.5*mean edge length(A); delta2=0.5*mean edge length(B);
delta3=unrounded discrete Frechet distance(A,B). Identical inputs, endpoints and
thresholds apply across methods. Objective: min max(kA,kB). Independent continuous
and discrete baselines enforce their respective global fidelity only. CPS-2F enforces
continuous fidelity plus discrete output coupling; CPS-3F enforces discrete fidelity
plus discrete coupling. All outputs use original vertices; unequal lengths and repeated
coupling indices are allowed. certificate=False disables both shortcuts. No greedy
processing, vertex replacement, threshold tuning or reference-guided stopping was added.

External solver60s, internal CPS55s, sampled tree RSS512MiB; CPS100000 states,
20000000 candidate checks and 800000 main single-graph edges. Independent discrete
and CPS-3F independent preprocessing retain the 2000000-edge default. Independent
continuous has no equivalent graph-work cap. These are shared external limits with
disclosed method-specific internal limits, not equal computational work guarantees.
Separate unchanged phase budgets: startup60, input15, warmup120, checkpoint5,
validation30, optional measurements30, reporting10 seconds. Internal soft time is checked
after state expansion; external monitoring is the backstop. No automatic retries.

## Observed size and coupling comparisons
{chr(10).join(quality)}

Only matched validated outputs enter vertex-cost comparisons. retained_A, retained_B
and retained_total record retained fractions. CPS-2F is compared against independent
continuous; CPS-3F against independent discrete. Direct CPS-2F/CPS-3F comparisons
change fidelity constraints as well as the algorithm. Missing outputs are not zero.
In the fixed AL031854/AL011991 illustration both independent continuous and CPS-2F
have k=5, but coupling is 725.571465 versus 426.184218 km for a 426.184218 km budget.
This demonstrates equal size with different coupling feasibility, not universal
zero vertex cost or topological preservation.

## Graph bounds and computational cost
{chr(10).join(graph_findings)}
{stop_details}

graph_bounds.md cites thesis Observation4 (printed p18/PDF26) and section5.1
(printed p26/PDF34), and manuscript sections3.1/4.1/4.4 (PDF7-8/11-12/15-16).
For CPS-2F, auxiliary bound U(n)=2n(n-1); finite potential-state bound
n*m*(n+U(n))*(m+U(m)). CPS-3F has no auxiliary points and bound n^2*m^2.
Actual extended-size bound is n*m*(n+aA)*(m+aB); completed active single-graph
counts give H_A*H_B. Neither product enumerates coupling-valid states. Source and
terminal are included. Auxiliary matching locations are not selectable output vertices.
Big-O expressions are not assigned a hidden constant of one.

graph_bounds.csv separates discovered states, dense single-graph slots, active states,
single-graph edges including waits, auxiliary locations and candidate transition checks.
Single-graph edge bound H(H+1)/2 counts the same objects (unique forward edges plus
H waits). Candidate-check bound E_A*E_B-H_A*H_B is derived separately; the implicit
product graph is not materialized, and candidate checks are not relabelled as unique
stored edges. Unsupported product-edge and DP-label observations remain missing.
Percentages compare matching counting quantities and are never completion estimates.
Interrupted counts are observed before termination. Historical missing counters are
not filled by rerunning solvers. No packedness or general scaling law is inferred.

Solver time includes existing internal fidelity/coupling checks. External validation
has its own timer; recovery is separately identified. Successful times, solver resource
exits and returned-but-unvalidated outputs use distinct markers. Inherited observations
are hollow. RSS is sampled whole-process-tree memory, including warmup, not exact
algorithm memory. Boundary graph counters and five-second supervisor status updates
add small unquantified overhead, so these single attempts are not a timing benchmark.
Feasibility validation alone does not establish optimality; solver optimality claims
remain separate, supported by prior small-instance exhaustive checks.

## Live progress and files
progress.json is atomically refreshed approximately every five seconds outside solver
loops; progress_history.jsonl retains updates. New slots finished /24 counts recorded
terminal outcomes, not just successes. Reused observations are separate. Graph counters
carry timestamps; budget-used percentages never mean search completion. The saved-only
presentation notebook displays progress and results without starting optimization.

[Results](results.csv), [completion](completion.csv), [matched comparisons](matched_comparisons.csv),
[graph counts and formulas](graph_bounds.csv), [compact state table](graph_state_summary.csv),
[captions](figure_captions.md), figures/*.pdf and *.png, and original per-attempt checkpoints.
Algorithms and scientific parameters are unchanged. common.py has the layout-only
validation fix; worker.py supports validation replay; supervisor.py has an optional
five-second observation callback; expanded_comparison.py manages the revised manifest
and reporting. No duplicate infrastructure or new versioned code directory was created.

The batch was resumed after an orchestration pause. Historical progress entries named
total_batch_elapsed_seconds measured their current execution segment. The final progress
record separates that segment from total wall time since manifest freeze (including
pauses and recovery); timing_definitions.json records this clarification. No solver
timers or historical progress records were rewritten.
'''
        if counts['unexecuted']:
            report_text+='\n## Pause\n\n'+json.dumps(read(OUT/'progress.json'),indent=2)+'\n'
    (OUT/'report.md').write_text(report_text,encoding='utf-8')
    (OUT/'status.md').write_text(f"# Expansion status\n\n{json.dumps(counts,indent=2)}\n\nSee report.md and progress.json for the current primary cohort status and validation recovery.\n",encoding='utf-8')
    notebook=ROOT/'notebooks_paper/expanded_comparison_saved_results.ipynb'
    nb=nbformat.v4.new_notebook(cells=[nbformat.v4.new_markdown_cell('# Expanded CPS comparison: saved results\n\nResearch question: how do fidelity/coupling constraints, bounded completion and graph sizes compare?\n\nThis notebook reads frozen results only. It never invokes optimization.'),
        nbformat.v4.new_code_cell("from pathlib import Path\nimport pandas as pd\nfrom IPython.display import display, Markdown, Image\nROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p/'output/paper_v1/expanded_comparison/results.csv').exists())\nSAVED = ROOT/'output/paper_v1/expanded_comparison'"),
        nbformat.v4.new_markdown_cell('## Inputs and definitions\n\nThe central parameter above selects the saved experiment directory. Scientific parameters are frozen in configuration_manifest.json. No notebook threshold overrides.'),
        nbformat.v4.new_code_cell("display(pd.read_csv(SAVED/'input_manifest.csv')[['domain','pair','nA','nB','subsampling']])\ndisplay(Markdown((SAVED/'report.md').read_text()))"),
        nbformat.v4.new_code_cell("display(pd.read_csv(SAVED/'completion.csv'))\nfor name in ['completion','solution_size','solver_time','memory','configuration_states','state_utilization','single_curve_sizes','single_curve_edge_bounds','algorithm_diagram','atlantic_coupling_example']:\n    display(Image(filename=str(SAVED/'figures'/f'{name}.png')))"),
        nbformat.v4.new_markdown_cell('## Supported findings and remaining gaps\n\nSee the saved report and captions. Resource-limited and unexecuted configurations remain explicit; feasibility checks do not independently prove optimality. Full tables: results.csv, matched_comparisons.csv and graph_bounds.csv.')])
    if primary:
        nb.cells.insert(2,nbformat.v4.new_code_cell("import json\np = json.loads((SAVED/'progress.json').read_text())\ndisplay(pd.Series({k:p.get(k) for k in ['utc','new_slots_finished','new_slots_total','new_slots_finished_percent','reused_observations','validated','resource_limited','validation_incomplete','failed','pending','phase','pair','method','pause_reason']}))\n# Re-execute this read-only cell to refresh; no optimization is launched."))
        nb.cells[0].source='# Primary CPS comparison: alpha=0.5\n\nTen frozen pairs, four methods; alpha selected using pilot evidence. This notebook reads saved artifacts only and never launches optimization.'
    nbformat.write(nb,notebook)
    print('REPORT',counts,flush=True)

TIMING=OUT/'timing_study'

def timing_prepare():
    import platform,importlib.metadata as meta,psutil
    runner.cohort_active_check();TIMING.mkdir(exist_ok=True)
    if (TIMING/'schedule.json').exists():
        print('Existing timing schedule retained; no new schedule or extra repetitions.');return
    base=read(OUT/'execution_manifest_alpha05.json')['rows'];assert len(base)==40
    pairs=list(dict.fromkeys(r['pair'] for r in base));schedule=[]
    for repetition in [1,2,3]:
        methods=METHODS[repetition-1:]+METHODS[:repetition-1]
        for pair in pairs:
            for method in methods:
                source=next(r for r in base if r['pair']==pair and r['method']==method)
                order=len(schedule)+1;identifier=source['configuration_id']+f'|repetition{repetition}'
                c={k:v for k,v in source['configuration'].items() if k not in ['code_version','attempt_id','test_mode','validation_replay_raw','validation_replay_sha256']}
                c.update(configuration_id=identifier,diagnostic=f'timing{order:03d}',certificate=False,graph_limits=LIMITS)
                r={k:v for k,v in source.items() if k not in ['configuration','reuse_path','original_origin']}
                r.update(plan_order=order,configuration_id=identifier,base_configuration_id=source['configuration_id'],repetition=repetition,
                    block_id=f'r{repetition}:{pair}',configuration=c,reuse_path=None)
                for p,h in c['algorithms'].items():assert sha(ROOT/p)==h
                assert sha(ROOT/c['input_file'])==c['input_sha256'] and c['alpha']==.5
                schedule.append(r)
    atomic(TIMING/'schedule.json',dict(rows=schedule,**stamp(),rounds=3,planned_repetitions=120,
        ordering='Fixed pair order from primary manifest; method order rotated left by round-1; one worker at a time. No outcome-based changes.',
        cache_policy='Fresh output-local cache for each repetition; same worker warmup and thread settings as completed cohort.',
        readiness_policy='Existing 20-second process-aware readiness check before each four-method pair/round block and after resume; active worker check before every attempt. Readiness is outside measured phases.',
        source_manifest_sha256=sha(OUT/'execution_manifest_alpha05.json'),original_results_sha256=sha(OUT/'results.csv'),
        environment=dict(executable=sys.executable,executable_sha256=sha(sys.executable),python=sys.version,platform=platform.platform(),
            logical_cpus=psutil.cpu_count(),memory_total_mib=psutil.virtual_memory().total/1024**2,
            threads={'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'},packages={n:meta.version(n) for n in ['numpy','numba','llvmlite','psutil','pandas','matplotlib']}),
        execution_code={n:sha(Path(__file__).with_name(n)) for n in ['common.py','worker.py','supervisor.py','controlled_followup.py','expanded_comparison.py']}))
    pd.DataFrame([{k:v for k,v in r.items() if k!='configuration'} for r in schedule]).to_csv(TIMING/'schedule.csv',index=False)
    for r in schedule:
        slot=TIMING/'slots'/f"{r['plan_order']:02d}";slot.mkdir(parents=True,exist_ok=True)
        atomic(slot/'slot.json',dict(configuration_id=r['configuration_id'],plan_order=r['plan_order'],reuse_path=None))
        atomic(slot/'configuration.json',r['configuration'])
    timing_save();print('Frozen timing schedule: 120 fresh repetitions, no original timings reused.',flush=True)

def timing_save():
    original_cohort=runner.COHORT;runner.COHORT=TIMING
    schedule=read(TIMING/'schedule.json')['rows'];records=[]
    try:
        for r in schedule:
            x=runner.cohort_result(r)
            x.update(base_configuration_id=r['base_configuration_id'],repetition=r['repetition'],block_id=r['block_id'])
            x['termination_reason']=x['status'];x['optional_measurement_status']=None
            if x['attempt_path']:
                d=ROOT/x['attempt_path'];sup=read(d/'supervisor.json') if (d/'supervisor.json').exists() else {}
                if sup.get('limiting_phase') and not x['validated']:
                    x['resource_limited']=True;x['termination_reason']=sup['status']
                    if x['status'] not in ['validation_failed','validation_incomplete']:x['status']=sup['status']
                if (d/'measurements.json').exists():x['optional_measurement_status']=read(d/'measurements.json')['status']
                x['configuration_sha256']=sha(d/'configuration.json')
            for k in ['kA','kB','k']:
                x['raw_'+k]=x[k]
                if not x['validated']:x[k]=None
            records.append(x)
    finally:runner.COHORT=original_cohort
    f=pd.DataFrame(records);primary=pd.read_csv(OUT/'results.csv').set_index('configuration_id')
    f['objective_consistency']='not_validated'
    for cid,sub in f.groupby('base_configuration_id'):
        successful=sub[sub.validated];values=set(successful.k)
        reference=primary.loc[cid]
        if bool(reference.validated):values.add(float(reference.k))
        if len(values)>1:f.loc[successful.index,'objective_consistency']='unexpected_objective_difference'
        else:f.loc[successful.index,'objective_consistency']='agrees_with_primary' if reference.validated else 'consistent_among_available_validated_repetitions'
    f.to_csv(TIMING/'results.csv',index=False)
    atomic(TIMING/'counts.json',dict(planned=120,finished=int((f.origin.eq('new')&f.status.ne('attempt_incomplete')).sum()),validated=int(f.validated.sum()),
        resource_limited=int((f.resource_limited&~f.status.eq('validation_incomplete')).sum()),validation_incomplete=int(f.status.eq('validation_incomplete').sum()),
        errors=int(f.status.isin(['software_error','validation_failed','attempt_incomplete']).sum()),pending=int(f.origin.eq('unexecuted').sum())))
    return f

def timing_execute():
    from supervisor import attempt
    runner.cohort_active_check();timing_prepare();schedule=read(TIMING/'schedule.json');last_ready_block=None
    for n,h in schedule['execution_code'].items():assert sha(Path(__file__).with_name(n))==h,('Timing dependencies changed',n)
    assert sys.executable==schedule['environment']['executable']
    with LiveProgress(frame_loader=timing_save,output=TIMING,total_new=120,batch_started_unix=schedule['unix']) as progress:
        for r in schedule['rows']:
            slot=TIMING/'slots'/f"{r['plan_order']:02d}";existing=list((slot/'attempts').glob('*/configuration.json'))
            if existing:
                assert len(existing)==1
                result=progress.frame[progress.frame.configuration_id.eq(r['configuration_id'])].iloc[0]
                if result.status in ['software_error','validation_failed','attempt_incomplete'] or result.objective_consistency=='unexpected_objective_difference':
                    progress.finish('paused','Existing systemic error/objective difference requires investigation; no duplicate execution');return
                continue
            if (slot/'launch.json').exists():
                progress.finish('paused','Launch marker without finalized attempt; no automatic retry');return
            runner.cohort_active_check();progress.current(r,'host_readiness')
            if last_ready_block!=r['block_id']:
                if not runner.cohort_readiness(slot):
                    progress.finish('paused','Sustained competing host activity; see readiness record');return
                last_ready_block=r['block_id']
            config=read(slot/'configuration.json');assert config==r['configuration']
            atomic(slot/'launch.json',dict(**stamp(),configuration_sha256=sha(slot/'configuration.json'),automatic_retries=False))
            progress.current(r,'startup');d=attempt(config,output_directory=slot,monitor=progress.observe)
            atomic(slot/'attempt_index.json',dict(attempt_directory=str(d.relative_to(ROOT))))
            progress.refresh();result=progress.frame[progress.frame.configuration_id.eq(r['configuration_id'])].iloc[0]
            print('TIMING END',r['plan_order'],r['repetition'],r['pair'],r['method'],result.status,result.k,flush=True)
            if result.status in ['software_error','validation_failed','attempt_incomplete'] or str(result.get('raw_status')).startswith('infeasible') or result.objective_consistency=='unexpected_objective_difference':
                atomic(TIMING/'pause.json',dict(**stamp(),configuration_id=r['configuration_id'],status=result.status,objective_consistency=result.objective_consistency))
                progress.finish('paused','Systemic software/feasibility/objective consistency investigation required');return
        progress.refresh();progress.finish('completed');atomic(TIMING/'completed.json',dict(**stamp(),counts=read(TIMING/'counts.json')))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','execute','report','primary_prepare','recover','primary_execute','timing_prepare','timing_execute']);parser.add_argument('--slot',type=int,default=33);args=parser.parse_args()
    if args.action=='recover':recover(args.slot)
    else:globals()[args.action]()
