"""Frozen coverage extension, using the existing worker, monitor and validator.

No optimizer implementation is defined here. No timing repetitions are launched.
"""
from pathlib import Path
import sys,json,time,argparse,shutil
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'notebooks_paper/graph_only_v3')]
import numpy as np
import pandas as pd
import expanded_comparison as existing
import controlled_followup as runner
from common import atomic,sha,stamp,BUDGETS
from supervisor import attempt
OUT=ROOT/'output/paper_v1/coverage_extension'
OLD=ROOT/'output/paper_v1/expanded_comparison'
METHODS=existing.METHODS
TARGETS=['4eca.b','4eca.c','3ntx.a','1wls.a','2eq5.a','2zsk.a','1zq1.a','3jq0.a','2fep.a']

RULE='''# Coverage extension: rule fixed before optimization

Alpha=0.5 only. Preserve the original ten pairs and their observations. Add the
nine requested paper protein identities; preparation eligibility is determined
from sources, not solver outcomes. Reuse the existing aligned w=16 arrays for
4eca.b and 4eca.c, including their documented omissions. Do not reconstruct them.
For new unequal-length inputs reuse small_pair_data.structural_fit (US-align,
rigid sequential structural superposition) on entire observed chains, then reuse
decimation_algorithms.decimate independently with w=16, retaining both endpoints.
No equal-length cropping, padding, guessed index correspondences or trimming to
the paper's length is allowed. Internal breaks in a newly prepared chain block
the pair: the existing residue preparation splits those chains, and using a
longest fragment or bridging gaps would introduce a different input policy.
Existing prepared gapped geometries remain explicitly identified as such.

Hurricanes: consider only unordered-unique pairs in the existing pairs.csv,
with finite 2D tracks, >=2 vertices, matching declared lengths, and matching
basin classifications (AL=Atlantic; EP/CP=Pacific). Retain every supplied row and the common basin frame.
Exclude the five existing pairs. Prefer unused storms. To cover ordinary sizes,
use the existing cohort's range of maximum track lengths when it contains enough
eligible pairs; split that range's candidate size ranks into five quantile bands.
Choose one pair from each band, in increasing band order. Within a band first
minimize reused storms, then maximize minimum standardized feature distance to
the original and already selected pairs, then use pair ID as deterministic tie
break. Features are log(max length), log(length ratio), mean chord/path ratio,
and centroid separation divided by mean path length. Existing feature() is reused.
No fidelity results, compression, runtime, or solver success enter selection.

Thresholds: delta1=0.5 mean_edge(A), delta2=0.5 mean_edge(B),
delta3=discrete_Frechet(A,B), full precision. All adopted budgets, thread settings,
fresh output-local caches, warm-up, validation and certificate=False are retained.
Run each eligible new slot once, sequentially. A resource exit is an observation.
Blocked inputs are neither observations nor infeasible instances. Readiness is
checked before each new configuration, outside measured solver phases.
'''

def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def manifest():return read(OUT/'configuration_manifest.json')

def prepare():
    runner.cohort_active_check()
    assert read(OLD/'timing_study/counts.json')['pending']==0,'Finish old timing schedule first'
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'configuration_manifest.json').exists():print('Frozen manifest retained.');return
    (OUT/'selection_rule.md').write_text(RULE,encoding='utf-8')
    for name in ['inputs','protein_audit','slots','figures','provenance']:(OUT/name).mkdir(exist_ok=True)
    import pymupdf
    paper=ROOT/'papers/chain_pair_simplification.pdf';doc=pymupdf.open(paper);pages=[]
    for i,page in enumerate(doc):
        t=page.get_text()
        if any(x in t for x in ['Table 1:','Table 2:','Table 3:']):pages.append(dict(pdf_page=i+1,text=t))
    alltext='\n'.join(x['text'] for x in pages)
    assert all(x in alltext for x in ['3ntx','1wls','2eq5','2zsk','1zq1','3jq0','2fep','1d9q'])
    atomic(OUT/'provenance/paper_tables.json',dict(path=str(paper.relative_to(ROOT)),sha256=sha(paper),pages=pages))
    (OUT/'provenance/paper_tables.txt').write_text(alltext,encoding='utf-8');doc.close()
    from protein_data import prepare as audit_proteins,paper_experiments
    from small_pair_data import structural_fit
    from decimation_algorithms import decimate
    from curve_algorithms import discrete_frechet
    chains,audit,breaks,missing=audit_proteins(ROOT)
    for f,name in [(audit,'identifiers'),(breaks,'chain_breaks'),(missing,'missing_residues'),(paper_experiments(),'paper_parameters_reference_only')]:f.to_csv(OUT/'protein_audit'/f'{name}.csv',index=False)
    for key,c in chains.items():c['mapping'].to_csv(OUT/'protein_audit'/f'{key}.residues.csv',index=False)
    curves_path=ROOT/'data/decimated protein backbones/decimated_curves.csv';curves=pd.read_csv(curves_path)
    pair_path=ROOT/'data/decimated protein backbones/decimated_pairs.csv'
    exe=ROOT/'tools/usalign/USalign.exe';assert exe.exists()
    def supplied(key):return curves[(curves.chain==key)&curves.frame.eq('aligned')&curves.w.eq(16)].sort_values('t')[['x','y','z']].to_numpy(float)
    oldinputs=read(OLD/'input_manifest.json')['inputs'];inputs=[];arrays={};eligibility=[]
    amap=chains['1o7j.a']['mapping'];Afull=chains['1o7j.a']['xyz'];assert len(chains['1o7j.a']['fragments'])==1
    for b in TARGETS:
        pair='1o7j.a__'+b;c=chains[b];row=audit.set_index('verified').loc[b]
        r=dict(domain='protein',pair=pair,individuals=['1o7j',b.split('.')[0]],nA_full=len(Afull),nB_full=len(c['xyz']),paper_reported_B=int(row.reported_length),
            B_fragments=int(row.fragments),A_author_chain='A',B_author_chain=row.author_chain,model=int(row.model),units='angstrom',subsampling='w=16 with last endpoint appended if necessary',
            existing_prepared=b in ['4eca.b','4eca.c'],eligible=True,blocked_reason=None)
        if r['existing_prepared']:
            A,B=supplied('1o7j.a'),supplied(b);assert len(A)==len(B)==22
            Adep=decimate(Afull,16);Bdep=decimate(c['xyz'],16)
            # Rigid preparation must preserve within-chain distances and endpoint selection.
            errors=[]
            for key,X,Y in [('1o7j.a',A,Adep),(b,B,Bdep)]:
                selected=curves[curves.chain.eq(key)&curves.frame.eq('aligned')&curves.w.eq(16)].sort_values('t')
                assert np.array_equal(selected.source_row.to_numpy(),decimate(np.arange(325),16))
                error=float(np.max(np.abs(np.linalg.norm(X[:,None]-X[None,:],axis=2)-np.linalg.norm(Y[:,None]-Y[None,:],axis=2))))
                # Source CSV stores aligned coordinates to 0.001 A; two-point
                # distance error is bounded by 2*sqrt(3)*0.0005, not solver EPS.
                assert error<=np.sqrt(3)*.001+1e-9
                errors.append(error)
            r.update(alignment='existing supplied index-Kabsch alignment; arrays reused exactly',gap_policy='historical supplied geometric chain has omitted residue 12; retained as supplied, not claimed biologically contiguous',source_coordinate_rounding_A=.001,source_rigid_distance_discrepancies=errors,source_files=[str(curves_path.relative_to(ROOT)),str(pair_path.relative_to(ROOT))])
        elif len(c['fragments'])>1:
            r.update(eligible=False,blocked_reason='Internal chain break(s); no existing uncropped whole-chain gap-handling policy. Existing preparation splits fragments; cropping/bridging not introduced.',alignment='not performed',gap_policy='blocked; no new gap bridging',nA=len(decimate(Afull,16)),nB=len(decimate(c['xyz'],16)))
            inputs.append(r);eligibility.append(r);continue
        else:
            folder=OUT/'inputs'/pair;folder.mkdir(exist_ok=True)
            Bfit,fit=structural_fit(Afull,c['xyz'],amap,c['mapping'],folder/'alignment',exe)
            ia=decimate(np.arange(len(Afull)),16);ib=decimate(np.arange(len(Bfit)),16)
            assert ia[0]==ib[0]==0 and ia[-1]==len(Afull)-1 and ib[-1]==len(Bfit)-1
            A,B=Afull[ia],Bfit[ib]
            assert len(fit['alignment_pairs'])>0 and np.allclose(np.linalg.norm(np.diff(Bfit,axis=0),axis=1),np.linalg.norm(np.diff(c['xyz'],axis=0),axis=1),atol=1e-5)
            np.savez(folder/'preparation.npz',A_full=Afull,B_deposited=c['xyz'],B_aligned=Bfit,indices_A=ia,indices_B=ib)
            r.update(alignment='existing US-align rigid sequential structural_fit applied to entire observed unequal-length chains before w=16',gap_policy='no internal breaks in selected observed chain; terminal missing residues documented',alignment_record=str((folder/'alignment/alignment.json').relative_to(ROOT)),source_files=[f'data/pdb_cache/{b.split(".")[0]}.cif','data/raw protein backbones/1o7j.a.raw'])
        r.update(nA=len(A),nB=len(B));arrays[pair]=(A,B);inputs.append(r);eligibility.append(r)
    pd.DataFrame(eligibility).to_csv(OUT/'protein_audit/pair_eligibility.csv',index=False)
    # Hurricane candidate construction uses coordinates only, never outcome columns.
    hp_path=ROOT/'data/hurricanes/pairs.csv';hp=pd.read_csv(hp_path);pool=[];seen=set()
    for sr in hp.sort_values(['curve_A','curve_B']).to_dict('records'):
        a,b=sr['curve_A'],sr['curve_B'];key=tuple(sorted([a,b]))
        if key in seen:continue
        seen.add(key);A,B=[np.loadtxt(ROOT/'data/hurricanes'/f'{z}.txt') for z in [a,b]]
        assert A.ndim==B.ndim==2 and A.shape[1]==B.shape[1]==2 and np.isfinite(A).all() and np.isfinite(B).all()
        assert len(A)==sr['n_A'] and len(B)==sr['n_B'] and min(len(A),len(B))>=2
        basin_for={'AL':'atlantic','EP':'pacific','CP':'pacific'}
        assert basin_for[a[:2]]==basin_for[b[:2]]==sr['basin']
        pool.append(dict(domain='hurricane',pair=a+'__'+b,nA=len(A),nB=len(B),individuals=[a,b],features=existing.feature(A,B),basin=sr['basin'],A=A,B=B))
    retained=[p for p in pool if p['pair'] in {r['pair'] for r in oldinputs if r['domain']=='hurricane'}];assert len(retained)==5
    X=np.asarray([p['features'] for p in pool]);sd=X.std(0);sd[sd<1e-12]=1
    for p,z in zip(pool,(X-X.mean(0))/sd):p['standardized_features']=z.tolist()
    low=min(max(p['nA'],p['nB']) for p in retained);high=max(max(p['nA'],p['nB']) for p in retained)
    candidates=[p for p in pool if p['pair'] not in {q['pair'] for q in retained} and low<=max(p['nA'],p['nB'])<=high]
    if len(candidates)<5:candidates=[p for p in pool if p['pair'] not in {x['pair'] for x in retained}]
    sizes=np.array([max(p['nA'],p['nB']) for p in candidates]);cut=np.quantile(sizes,[.2,.4,.6,.8])
    for p in candidates:p['size_band']=int(np.searchsorted(cut,max(p['nA'],p['nB']),side='left'))
    chosen=list(retained);used={z for p in chosen for z in p['individuals']};new=[]
    for band in range(5):
        available=[p for p in candidates if p['pair'] not in {x['pair'] for x in chosen} and p['size_band']==band]
        if not available:available=[p for p in candidates if p['pair'] not in {x['pair'] for x in chosen}]
        def score(p):return (sum(z in used for z in p['individuals']),-min(float(np.linalg.norm(np.array(p['standardized_features'])-q['standardized_features'])) for q in chosen),p['pair'])
        pick=min(available,key=score);pick['selection_score']=score(pick);chosen.append(pick);new.append(pick);used.update(pick['individuals'])
    atomic(OUT/'hurricane_selection.json',dict(size_range=[low,high],band_cutpoints=cut.tolist(),candidates=[{k:v for k,v in p.items() if k not in ['A','B']}|dict(selected=p['pair'] in {x['pair'] for x in new}) for p in candidates],rule_sha256=sha(OUT/'selection_rule.md')))
    for p in new:
        A,B=p['A'],p['B'];arrays[p['pair']]=(A,B)
        inputs.append({k:v for k,v in p.items() if k not in ['A','B']}|dict(eligible=True,blocked_reason=None,units='projected km',subsampling='none; every supplied track vertex',alignment='unchanged common basin frame; no independent alignment',source_files=[f'data/hurricanes/{x}.txt' for x in p['individuals']]))
    oldrows=read(OLD/'execution_manifest_alpha05.json')['rows'];template=oldrows[0]['configuration']
    alg=template['algorithms'];rows=[]
    for r in inputs:
        if r['eligible']:
            A,B=arrays[r['pair']];path=OUT/'inputs'/f"{r['domain']}_{r['pair']}.npz";np.savez(path,A=np.ascontiguousarray(A),B=np.ascontiguousarray(B))
            r.update(input_file=str(path.relative_to(ROOT)),input_sha256=sha(path),s_A=float(np.linalg.norm(np.diff(A,axis=0),axis=1).mean()),s_B=float(np.linalg.norm(np.diff(B,axis=0),axis=1).mean()),delta3=float(discrete_frechet(A,B)))
            r.update(delta1=.5*r['s_A'],delta2=.5*r['s_B'])
            r['source_hashes']={p:sha(ROOT/p) for p in r['source_files']}
        for method in METHODS:
            order=len(rows)+1;cid=f"{r['domain']}:{r['pair']}:a0.5|{method}"
            c=None
            if r['eligible']:
                c={k:v for k,v in template.items() if k not in ['code_version','attempt_id','validation_replay_raw','validation_replay_sha256']}
                c.update({k:r[k] for k in ['domain','pair','nA','nB','units','input_file','input_sha256','s_A','s_B','delta1','delta2','delta3']})
                c.update(alpha=.5,method=method,configuration_id=cid,case=cid.split('|')[0],diagnostic=f'coverage{order:02d}',certificate=False,graph_limits=existing.LIMITS,phase_budgets=BUDGETS,memory_mib=512)
            row={k:r.get(k) for k in ['domain','pair','nA','nB','units','input_file','input_sha256','delta1','delta2','delta3','eligible','blocked_reason']}
            row.update(plan_order=order,configuration_id=cid,alpha=.5,method=method,configuration=c,reuse_path=None)
            rows.append(row)
    assert len(inputs)==14 and len(rows)==56
    for name,h in alg.items():assert sha(ROOT/name)==h
    atomic(OUT/'input_manifest.json',dict(inputs=inputs,**stamp()))
    pd.DataFrame(inputs).to_csv(OUT/'input_manifest.csv',index=False)
    atomic(OUT/'configuration_manifest.json',dict(rows=rows,algorithms=alg,source_old_manifest=sha(OLD/'execution_manifest_alpha05.json'),
        old_timing_completed=read(OLD/'timing_study/completed.json'),old_pending_slots=[117,118,119,120],selection_rule_sha256=sha(OUT/'selection_rule.md'),
        execution_dependencies={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'protein_data.py',ROOT/'small_pair_data.py',ROOT/'decimation_algorithms.py',exe]+list((ROOT/'notebooks_paper/graph_only_v3').glob('*.py'))},**stamp()))
    pd.DataFrame([{k:v for k,v in r.items() if k!='configuration'}|dict(configuration_json=json.dumps(r['configuration'])) for r in rows]).to_csv(OUT/'run_manifest.csv',index=False)
    for r in rows:
        slot=OUT/'slots'/f"{r['plan_order']:02d}";slot.mkdir(exist_ok=True)
        atomic(slot/'slot.json',dict(configuration_id=r['configuration_id'],plan_order=r['plan_order'],reuse_path=None))
        if r['eligible']:atomic(slot/'configuration.json',r['configuration'])
    save();print('FROZEN',len(rows),'new slots;',sum(r['eligible'] for r in rows),'eligible;',sum(not r['eligible'] for r in rows),'blocked',flush=True)
    print(pd.DataFrame(inputs)[['domain','pair','nA','nB','eligible']].to_string(index=False),flush=True)

def save():
    runner.COHORT=OUT;records=[]
    for r in manifest()['rows']:
        if not r['eligible']:
            x={k:v for k,v in r.items() if k!='configuration'}
            x.update(origin='blocked',status='preparation_blocked',termination_reason='preparation_blocked',validated=False,resource_limited=False,attempt_path=None,kA=None,kB=None,k=None,solver_seconds=None)
        else:
            x=runner.cohort_result(r);x['termination_reason']=x['status']
            for key in ['kA','kB','k']:
                x['raw_'+key]=x[key]
                if not x['validated']:x[key]=None
        x['observation_cohort']=OUT.name+' single observation';records.append(x)
    f=pd.DataFrame(records);f.to_csv(OUT/'results.csv',index=False)
    counts=dict(planned=len(f),eligible=int(sum(r['eligible'] for r in manifest()['rows'])),finished=int(f.origin.eq('new').sum()),validated=int(f.validated.sum()),
        resource_limited=int(f.resource_limited.sum()),validation_incomplete=int(f.status.eq('validation_incomplete').sum()),blocked=int(f.origin.eq('blocked').sum()),pending=int(f.origin.eq('unexecuted').sum()))
    atomic(OUT/'counts.json',counts);return f

class Progress(existing.LiveProgress):
    def __init__(self):
        self.host_wait_total=0.;self.host_wait_start=None;self.next_readiness=None
        super().__init__(frame_loader=save,output=OUT,total_new=len(manifest()['rows']),batch_started_unix=manifest()['unix'])
    def publish(self):
        with self.lock:
            waiting=0. if self.host_wait_start is None else time.monotonic()-self.host_wait_start
            self.data.update(continuous_host_wait_seconds=waiting,cumulative_host_wait_seconds=self.host_wait_total+waiting,
                continuous_host_wait_budget_seconds=900,cumulative_host_wait_budget_seconds=3600,next_readiness_check_unix=self.next_readiness,
                execution_state='waiting_for_host' if self.host_wait_start is not None else ('completed' if str(self.data.get('phase','')).startswith('completed') else 'paused' if self.data.get('phase')=='paused' else 'running'))
            launched_pending=sum((OUT/'slots'/f"{int(r.plan_order):02d}"/'launch.json').exists() for r in self.frame[self.frame.origin.eq('unexecuted')].itertuples())
            self.data.update(pending=int(self.frame.origin.eq('unexecuted').sum())-launched_pending,in_progress=launched_pending)
            self.data.update(old_pending_repetitions_completed=4,old_pending_repetitions_total=4,blocked=int(self.frame.origin.eq('blocked').sum()),
                new_configurations_planned=len(manifest()['rows']),new_configurations_eligible=int(sum(r['eligible'] for r in manifest()['rows'])))
            super().publish()
            print('OLD pending repetitions 4/4 completed; preparation-blocked new slots:',self.data['blocked'],flush=True)
            print('HOST WAIT continuous=%.1fs/900s cumulative=%.1fs/3600s next_check=%s'%(waiting,self.host_wait_total+waiting,self.next_readiness),flush=True)

def await_host(slot,progress):
    """Unchanged 20-second readiness samples, starts about 30 seconds apart.

    Waiting is outside every worker phase; no optimizer retry occurs here.
    Reserve the complete readiness sample within each waiting allowance.
    """
    while True:
        now=time.monotonic()
        used=0. if progress.host_wait_start is None else now-progress.host_wait_start
        if min(900-used,3600-progress.host_wait_total-used)<30:
            progress.host_wait_total+=used;progress.host_wait_start=None;progress.next_readiness=None
            progress.finish('paused','Bounded host waiting exhausted (15 min continuous / 60 min session); pending slot never launched')
            return False
        started=time.monotonic();progress.next_readiness=time.time()
        ready=runner.cohort_readiness(slot)
        if ready:
            if progress.host_wait_start is not None:progress.host_wait_total+=time.monotonic()-progress.host_wait_start
            progress.host_wait_start=None;progress.next_readiness=None
            return True
        if progress.host_wait_start is None:progress.host_wait_start=started
        with progress.lock:
            progress.data.update(phase='waiting_for_host',phase_started=progress.host_wait_start,phase_budget_seconds=900,
                pause_reason='Readiness failed; measurements saved in current slot; automatic bounded recheck')
        delay=max(0.,30-(time.monotonic()-started));progress.next_readiness=time.time()+delay
        progress.publish()
        # Short interruptible sleeps; the existing progress thread publishes every 5s.
        until=time.monotonic()+delay
        while time.monotonic()<until:time.sleep(min(1.,until-time.monotonic()))

def execute():
    runner.cohort_active_check();m=manifest()
    if (OUT/'completed.json').exists():
        print('Already completed; no execution:',read(OUT/'completed.json'),flush=True);return
    amendment=read(OUT/'provenance/bounded_wait_amendment.json') if (OUT/'provenance/bounded_wait_amendment.json').exists() else {}
    for p,h in m['execution_dependencies'].items():
        actual=sha(ROOT/p)
        assert actual==h or (p==amendment.get('path') and h==amendment.get('frozen_sha256') and actual==amendment.get('new_sha256')),('Changed frozen dependency',p)
    with Progress() as progress:
        for r in m['rows']:
            if not r['eligible']:continue
            slot=OUT/'slots'/f"{r['plan_order']:02d}";prior=list((slot/'attempts').glob('*/configuration.json'))
            if prior:
                assert len(prior)==1
                outcome=progress.frame.set_index('configuration_id').loc[r['configuration_id']]
                if outcome.status in ['software_error','validation_failed','attempt_incomplete']:
                    progress.finish('paused','Existing attempt requires investigation; no optimizer repeat');return
                continue
            if (slot/'launch.json').exists():progress.finish('paused','Launch marker requires investigation; no duplicate');return
            runner.cohort_active_check();progress.current(r,'host_readiness')
            if not await_host(slot,progress):return
            c=read(slot/'configuration.json');assert c==r['configuration'] and sha(ROOT/c['input_file'])==c['input_sha256']
            atomic(slot/'launch.json',dict(**stamp(),configuration_sha256=sha(slot/'configuration.json'),automatic_retries=False))
            progress.current(r,'startup');d=attempt(c,output_directory=slot,monitor=progress.observe)
            atomic(slot/'attempt_index.json',dict(attempt_directory=str(d.relative_to(ROOT))))
            progress.refresh();outcome=progress.frame.set_index('configuration_id').loc[r['configuration_id']]
            print('COVERAGE END',r['plan_order'],r['pair'],r['method'],outcome.status,outcome.k,flush=True)
            if outcome.status in ['software_error','validation_failed','attempt_incomplete'] or str(outcome.get('raw_status')).startswith('infeasible'):
                progress.finish('paused','Systemic software/feasibility investigation required');return
        progress.refresh();progress.finish('completed_with_preparation_gaps' if any(not r['eligible'] for r in m['rows']) else 'completed')
        atomic(OUT/'completed.json',dict(**stamp(),counts=read(OUT/'counts.json')))

def prepare_hurricane3():
    """Same coordinate-only selection rule; three size bands for three new pairs."""
    from curve_algorithms import discrete_frechet
    runner.cohort_active_check()
    prior=ROOT/'output/paper_v1/coverage_extension'
    assert OUT.name=='hurricane_extension_3'
    done=read(prior/'counts.json');assert done['pending']==0 and done['finished']==done['eligible']==48
    if (OUT/'configuration_manifest.json').exists():print('Frozen manifest already exists; preserved');return
    for sub in ['inputs','slots','provenance']:(OUT/sub).mkdir(parents=True,exist_ok=True)
    rule=('Existing coordinate-only selection procedure, with three quantile bands for three requested pairs. '
          'Use unordered-unique eligible pairs from pairs.csv, full finite 2D tracks with matching source lengths, common basin, positive mean edge length. '
          'Exclude all ten existing hurricane pairs. Use the current min/max maximum-track-length range. '
          'Within each size band minimize reused storm identities, then maximize the minimum standardized feature distance to existing and newly selected pairs, then break ties by pair ID. '
          'Features: log(max length), log(length ratio), mean chord/path ratio, centroid separation / mean path length. '
          'No solver outcomes or timing enter selection. No new preprocessing or additional subsampling. '
          'Alpha=0.5; all thresholds, phase budgets, graph budgets, thread settings, warmup and cache policy unchanged.\n')
    (OUT/'selection_rule.md').write_text(rule,encoding='utf-8')
    existing_inputs=read(OLD/'input_manifest.json')['inputs']+read(prior/'input_manifest.json')['inputs']
    existing_pairs={r['pair'] for r in existing_inputs if r['domain']=='hurricane'};assert len(existing_pairs)==10
    pool=[];seen=set();source=ROOT/'data/hurricanes/pairs.csv'
    for sr in pd.read_csv(source).sort_values(['curve_A','curve_B']).to_dict('records'):
        a,b=sr['curve_A'],sr['curve_B'];key=tuple(sorted([a,b]))
        if key in seen:continue
        seen.add(key);A,B=[np.loadtxt(ROOT/'data/hurricanes'/f'{z}.txt') for z in [a,b]]
        assert A.ndim==B.ndim==2 and A.shape[1]==B.shape[1]==2 and np.isfinite(A).all() and np.isfinite(B).all()
        assert len(A)==sr['n_A'] and len(B)==sr['n_B'] and min(len(A),len(B))>=2
        assert {'AL':'atlantic','EP':'pacific','CP':'pacific'}[a[:2]]=={'AL':'atlantic','EP':'pacific','CP':'pacific'}[b[:2]]==sr['basin']
        assert all(np.linalg.norm(np.diff(P,axis=0),axis=1).mean()>0 for P in [A,B])
        pool.append(dict(domain='hurricane',pair=a+'__'+b,nA=len(A),nB=len(B),individuals=[a,b],basin=sr['basin'],features=existing.feature(A,B),A=A,B=B))
    retained=[p for p in pool if p['pair'] in existing_pairs];assert len(retained)==10
    X=np.asarray([p['features'] for p in pool]);sd=X.std(0);sd[sd<1e-12]=1
    for p,z in zip(pool,(X-X.mean(0))/sd):p['standardized_features']=z.tolist()
    low=min(max(p['nA'],p['nB']) for p in retained);high=max(max(p['nA'],p['nB']) for p in retained)
    candidates=[p for p in pool if p['pair'] not in existing_pairs and low<=max(p['nA'],p['nB'])<=high];assert len(candidates)>=3
    cuts=np.quantile([max(p['nA'],p['nB']) for p in candidates],[1/3,2/3])
    for p in candidates:p['size_band']=int(np.searchsorted(cuts,max(p['nA'],p['nB']),side='left'))
    chosen=list(retained);used={z for p in chosen for z in p['individuals']};new=[]
    for band in range(3):
        available=[p for p in candidates if p['pair'] not in {q['pair'] for q in chosen} and p['size_band']==band]
        if not available:available=[p for p in candidates if p['pair'] not in {q['pair'] for q in chosen}]
        def score(p):return (sum(z in used for z in p['individuals']),-min(float(np.linalg.norm(np.array(p['standardized_features'])-q['standardized_features'])) for q in chosen),p['pair'])
        pick=min(available,key=score);pick['selection_score']=score(pick);chosen.append(pick);new.append(pick);used.update(pick['individuals'])
    atomic(OUT/'hurricane_selection.json',dict(size_range=[low,high],band_cutpoints=cuts.tolist(),existing_pairs=sorted(existing_pairs),feature_mean=X.mean(0).tolist(),feature_std=sd.tolist(),candidates=[{k:v for k,v in p.items() if k not in ['A','B']}|dict(selected=p in new) for p in candidates],rule_sha256=sha(OUT/'selection_rule.md')))
    template=next(r['configuration'] for r in read(prior/'configuration_manifest.json')['rows'] if r['eligible']);inputs=[];rows=[]
    for p in new:
        A,B=p['A'],p['B'];path=OUT/'inputs'/f"hurricane_{p['pair']}.npz";np.savez(path,A=np.ascontiguousarray(A),B=np.ascontiguousarray(B))
        r={k:v for k,v in p.items() if k not in ['A','B']}
        r.update(eligible=True,blocked_reason=None,units='projected km',subsampling='none; every supplied track vertex',alignment='unchanged common basin frame; no independent alignment',source_files=[f'data/hurricanes/{z}.txt' for z in p['individuals']],input_file=str(path.relative_to(ROOT)),input_sha256=sha(path),s_A=float(np.linalg.norm(np.diff(A,axis=0),axis=1).mean()),s_B=float(np.linalg.norm(np.diff(B,axis=0),axis=1).mean()),delta3=float(discrete_frechet(A,B)))
        r.update(delta1=.5*r['s_A'],delta2=.5*r['s_B']);r['source_hashes']={z:sha(ROOT/z) for z in r['source_files']};inputs.append(r)
        for method in METHODS:
            order=len(rows)+1;cid=f"hurricane:{r['pair']}:a0.5|{method}";c=json.loads(json.dumps(template))
            c.update({k:r[k] for k in ['domain','pair','nA','nB','units','input_file','input_sha256','s_A','s_B','delta1','delta2','delta3']})
            c.update(method=method,case=cid.split('|')[0],configuration_id=cid,diagnostic=f'hurricane3_{order:02d}')
            assert c['certificate'] is False and c['alpha']==.5 and c['graph_limits']==existing.LIMITS and c['phase_budgets']==BUDGETS
            rows.append({k:r.get(k) for k in ['domain','pair','nA','nB','units','input_file','input_sha256','delta1','delta2','delta3','eligible','blocked_reason']}|dict(plan_order=order,configuration_id=cid,alpha=.5,method=method,configuration=c,reuse_path=None))
    atomic(OUT/'input_manifest.json',dict(inputs=inputs,**stamp()));pd.DataFrame(inputs).to_csv(OUT/'input_manifest.csv',index=False)
    atomic(OUT/'configuration_manifest.json',dict(rows=rows,algorithms=template['algorithms'],prior_completed=read(prior/'completed.json'),source_manifest_sha256=sha(prior/'configuration_manifest.json'),selection_rule_sha256=sha(OUT/'selection_rule.md'),execution_dependencies={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__)]+list((ROOT/'notebooks_paper/graph_only_v3').glob('*.py'))},orchestration_change='Manifest-selectable output and schedule length only; worker/validator/algorithms unchanged',**stamp()))
    pd.DataFrame([{k:v for k,v in r.items() if k!='configuration'}|dict(configuration_json=json.dumps(r['configuration'])) for r in rows]).to_csv(OUT/'run_manifest.csv',index=False)
    for r in rows:
        slot=OUT/'slots'/f"{r['plan_order']:02d}";slot.mkdir(exist_ok=True);atomic(slot/'slot.json',dict(configuration_id=r['configuration_id'],plan_order=r['plan_order'],reuse_path=None));atomic(slot/'configuration.json',r['configuration'])
    save();print('FROZEN BEFORE SOLVING',[(p['pair'],p['nA'],p['nB'],p['selection_score']) for p in inputs],flush=True)

def prepare_final20():
    """Seeded stratification from input metadata only; no outcome-dependent selection."""
    from curve_algorithms import discrete_frechet
    runner.cohort_active_check();assert OUT.name=='hurricane_final_20'
    if (OUT/'configuration_manifest.json').exists():print('Frozen manifest exists; preserved');return
    for sub in ['inputs','slots','provenance']:(OUT/sub).mkdir(parents=True,exist_ok=True)
    previous=[OLD,ROOT/'output/paper_v1/coverage_extension',ROOT/'output/paper_v1/hurricane_extension_3']
    oldinputs=sum([read(p/'input_manifest.json')['inputs'] for p in previous],[])
    tested={tuple(sorted(r['pair'].split('__'))) for r in oldinputs if r['domain']=='hurricane'}
    assert len(tested)==13
    used_prior={z for pair in tested for z in pair};used=set(used_prior)
    rule='''# Final hurricane cohort: selection frozen before optimization
Seed: 20260916, NumPy default_rng (PCG64), candidates sorted by identifiers.
Sampling frame: unordered pairs in the supplied pairs.csv; not all possible storm combinations.
Eligibility: finite two-dimensional full tracks, at least two vertices, declared sizes matching actual sizes, compatible AL/EP/CP basin frame, and positive mean edge lengths. No upper size cap. Previously tested unordered pairs and duplicate source rows are excluded with reasons.
Quartiles: numpy.quantile(max(nA,nB), [.25,.5,.75]); ties enter the lower interval via searchsorted(side=left). Each quartile is subdivided by ratio <=1.5 versus >1.5.
Allocation: round-robin one slot at a time across lexicographically sorted nonempty strata until 20 slots are allocated. Strata at capacity are skipped, redistributing deficits deterministically among remaining nonempty strata; empty strata are retained with zero allocation. No selected pair is replaced after solving.
Selection: round-robin across allocated strata. Within a stratum first minimize the number of storm identities already used in previous or current cohorts, then choose uniformly using the seeded generator among tied candidates in sorted pair order. Record every reuse. This preference-constrained stratified sample is not a simple random sample; unweighted percentages describe the selected cohort, not the archive.
Method order rotates by pair: the four existing method names shifted by pair index modulo four. One attempt each, no repeats. Alpha=.5 was chosen using earlier pilot evidence. All scientific and resource settings remain unchanged.
'''
    (OUT/'selection_rule.md').write_text(rule,encoding='utf-8');pool=[];audit=[];seen=set()
    source=ROOT/'data/hurricanes/pairs.csv'
    for sr in pd.read_csv(source).sort_values(['curve_A','curve_B']).to_dict('records'):
        a,b=sr['curve_A'],sr['curve_B'];key=tuple(sorted([a,b]));reason=None;A=B=None
        if key in seen:reason='duplicate unordered source pair'
        elif key in tested:reason='previously tested unordered pair'
        seen.add(key)
        try:
            A,B=[np.loadtxt(ROOT/'data/hurricanes'/f'{z}.txt') for z in [a,b]]
            if not all(P.ndim==2 and P.shape[1]==2 and len(P)>=2 and np.isfinite(P).all() for P in [A,B]):reason=reason or 'invalid finite 2D track / fewer than two vertices'
            elif len(A)!=sr['n_A'] or len(B)!=sr['n_B']:reason=reason or 'declared size mismatch'
            elif {'AL':'atlantic','EP':'pacific','CP':'pacific'}.get(a[:2])!={'AL':'atlantic','EP':'pacific','CP':'pacific'}.get(b[:2]) or {'AL':'atlantic','EP':'pacific','CP':'pacific'}.get(a[:2])!=sr['basin']:reason=reason or 'incompatible basin frame'
            elif min(np.linalg.norm(np.diff(P,axis=0),axis=1).mean() for P in [A,B])<=0:reason=reason or 'zero mean edge length'
        except (OSError,ValueError) as error:reason=reason or 'unreadable track: '+str(error)
        record=dict(pair=a+'__'+b,individuals=[a,b],domain='hurricane',basin=sr['basin'],nA=len(A) if A is not None else None,nB=len(B) if B is not None else None,eligible=reason is None,exclusion_reason=reason)
        audit.append(record.copy())
        if reason is None:pool.append(record|dict(A=A,B=B,max_length=max(len(A),len(B)),length_ratio=max(len(A),len(B))/min(len(A),len(B))))
    assert len(pool)>=20,'Insufficient eligible pairs; do not substitute duplicates'
    cuts=np.quantile([p['max_length'] for p in pool],[.25,.5,.75]);strata={f'Q{q}_{ratio}':[] for q in range(1,5) for ratio in ['balanced','unequal']}
    for p in pool:
        p['size_quartile']=int(np.searchsorted(cuts,p['max_length'],side='left'))+1;p['ratio_group']='balanced' if p['length_ratio']<=1.5 else 'unequal';p['stratum']=f"Q{p['size_quartile']}_{p['ratio_group']}";strata[p['stratum']].append(p)
    allocation={s:0 for s in strata};nonempty=[s for s in strata if strata[s]]
    while sum(allocation.values())<20:
        for s in nonempty:
            if sum(allocation.values())==20:break
            if allocation[s]<len(strata[s]):allocation[s]+=1
    rng=np.random.default_rng(20260916);chosen=[];remaining={s:list(v) for s,v in strata.items()};taken={s:0 for s in strata}
    while len(chosen)<20:
        for s in nonempty:
            if taken[s]>=allocation[s]:continue
            candidates=remaining[s];reuse=min(sum(z in used for z in p['individuals']) for p in candidates)
            tied=[p for p in candidates if sum(z in used for z in p['individuals'])==reuse]
            p=tied[int(rng.integers(len(tied)))];remaining[s]=[q for q in candidates if q['pair']!=p['pair']]
            p['reused_prior_storms']=[z for z in p['individuals'] if z in used_prior];p['reused_current_storms']=[z for z in p['individuals'] if z in used and z not in used_prior];p['selection_order']=len(chosen)+1;p['tie_pool_size']=len(tied)
            chosen.append(p);used.update(p['individuals']);taken[s]+=1
    selected={p['pair'] for p in chosen}
    atomic(OUT/'selection_pool.json',dict(seed=20260916,rng='NumPy default_rng PCG64',source_sha256=sha(source),quartile_cutpoints=cuts.tolist(),eligible_pairs=len(pool),track_size_range=[min(min(p['nA'],p['nB']) for p in pool),max(p['max_length'] for p in pool)],allocation=allocation,stratum_sizes={s:len(v) for s,v in strata.items()},prior_pairs=[list(p) for p in sorted(tested)],pool=[{k:v for k,v in p.items() if k not in ['A','B']}|dict(selected=p['pair'] in selected) for p in pool],**stamp()))
    pd.DataFrame(audit).to_csv(OUT/'inclusion_exclusion.csv',index=False)
    template=read(previous[-1]/'configuration_manifest.json')['rows'][0]['configuration'];inputs=[];rows=[]
    for p in chosen:
        A,B=p['A'],p['B'];path=OUT/'inputs'/f"hurricane_{p['pair']}.npz";np.savez(path,A=np.ascontiguousarray(A),B=np.ascontiguousarray(B));r={k:v for k,v in p.items() if k not in ['A','B']}
        r.update(blocked_reason=None,units='projected km',subsampling='none; every supplied track vertex',alignment='unchanged shared basin frame',source_files=[f'data/hurricanes/{z}.txt' for z in p['individuals']],input_file=str(path.relative_to(ROOT)),input_sha256=sha(path),s_A=float(np.linalg.norm(np.diff(A,axis=0),axis=1).mean()),s_B=float(np.linalg.norm(np.diff(B,axis=0),axis=1).mean()),delta3=float(discrete_frechet(A,B)))
        r.update(delta1=.5*r['s_A'],delta2=.5*r['s_B']);r['source_hashes']={z:sha(ROOT/z) for z in r['source_files']};inputs.append(r)
        shift=(r['selection_order']-1)%4;order=METHODS[shift:]+METHODS[:shift]
        for method in order:
            cid=f"hurricane:{r['pair']}:a0.5|{method}";c=json.loads(json.dumps(template));c.update({k:r[k] for k in ['domain','pair','nA','nB','units','input_file','input_sha256','s_A','s_B','delta1','delta2','delta3']});c.update(method=method,case=cid.split('|')[0],configuration_id=cid,diagnostic=f'final20_{len(rows)+1:02d}')
            assert c['certificate'] is False and c['alpha']==.5 and c['graph_limits']==existing.LIMITS and c['phase_budgets']==BUDGETS
            rows.append({k:r.get(k) for k in ['domain','pair','nA','nB','units','input_file','input_sha256','delta1','delta2','delta3','eligible','blocked_reason','stratum','size_quartile','ratio_group','selection_order']}|dict(plan_order=len(rows)+1,configuration_id=cid,alpha=.5,method=method,configuration=c,reuse_path=None))
    atomic(OUT/'input_manifest.json',dict(inputs=inputs,**stamp()));pd.DataFrame(inputs).to_csv(OUT/'input_manifest.csv',index=False)
    atomic(OUT/'configuration_manifest.json',dict(rows=rows,algorithms=template['algorithms'],seed=20260916,selection_pool_sha256=sha(OUT/'selection_pool.json'),selection_rule_sha256=sha(OUT/'selection_rule.md'),execution_dependencies={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__)]+list((ROOT/'notebooks_paper/graph_only_v3').glob('*.py'))},**stamp()))
    pd.DataFrame([{k:v for k,v in r.items() if k!='configuration'}|dict(configuration_json=json.dumps(r['configuration'])) for r in rows]).to_csv(OUT/'run_manifest.csv',index=False)
    for r in rows:
        slot=OUT/'slots'/f"{r['plan_order']:02d}";slot.mkdir(exist_ok=True);atomic(slot/'slot.json',dict(configuration_id=r['configuration_id'],plan_order=r['plan_order'],reuse_path=None));atomic(slot/'configuration.json',r['configuration'])
    save();print('FROZEN 80 slots',[(p['pair'],p['nA'],p['nB'],p['stratum']) for p in inputs],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','prepare_hurricane3','prepare_final20','execute','save']);p.add_argument('--cohort',choices=['coverage_extension','hurricane_extension_3','hurricane_final_20'],default='coverage_extension');args=p.parse_args()
    OUT=ROOT/'output/paper_v1'/args.cohort;globals()[args.action]()
