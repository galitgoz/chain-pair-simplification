"""Saved-results-only Stage 1 audit and report. No optimization."""
from paper_core import ROOT,OUT,METHODS,sha,COLORS
import pandas as pd,numpy as np,json,time,nbformat,platform
import matplotlib.pyplot as plt
frames=[]
for domain in ['protein','hurricane']:
    p=OUT/f'{domain}_results.csv';d=pd.read_csv(p)
    for idx,row in d[d.status.eq('software_error')].iterrows():
        stem=row['case'].replace(':','_')+'_'+row['method'].replace(' ','_')
        log=OUT/'logs'/f'{stem}.log'
        if row.measured_execution_count==0 and log.exists() and 'ResourceLimit: wall-time budget exceeded' in log.read_text():
            d.loc[idx,'status']='startup_resource_limit'
            d.loc[idx,'failure_reason']='Warm-up hit the 8-second soft limit before the pilot solve; see '+str(log.relative_to(OUT))
    if domain=='hurricane':d['resource_monitoring_valid']=True
    d['peak_memory_scope']='sum of RSS over launcher and descendants, including imports/warmup; sampled' if domain=='hurricane' else 'unavailable: original launcher-only monitor invalid'
    for col in ['kA','kB','k','compression_removed_pct','continuous_A_lower','continuous_A_upper','continuous_B_lower','continuous_B_upper','coupling_distance','coupling_budget_ratio','vertex_cost_vs_matching_baseline','solve_seconds']:
        if col not in d:d[col]=np.nan
    d.to_csv(p,index=False);frames.append(d)
allrows=pd.concat(frames,ignore_index=True)
assert len(allrows)==32 and not allrows.duplicated(['case','method']).any()
assert allrows.groupby('case')[['delta1','delta2','delta3']].nunique().eq(1).all().all()
assert allrows.groupby('case').method.nunique().eq(4).all()
for r in allrows[allrows.status.eq('optimal')].itertuples():
    assert r.k==max(r.kA,r.kB) and r.constraints_valid
    if r.method.startswith('CPS'):assert r.coupling_distance<=r.delta3+1e-9
    if r.route in ['independent_certificate','alternate_certificate']:assert r.k==r.lower_bound
allrows.to_csv(OUT/'pilot_results.csv',index=False)
summary=allrows[['domain','pair','alpha','method','k','status']].copy()
summary['entry']=[str(int(k)) if s=='optimal' else 'missing: '+s for k,s in zip(summary.k,summary.status)]
wide=summary.pivot(index=['domain','pair','alpha'],columns='method',values='entry').reindex(columns=METHODS)
wide.to_csv(OUT/'pilot_summary.csv')
validation=pd.read_csv(OUT/'validation.csv')
fig,ax=plt.subplots(figsize=(7,2.6),facecolor='white')
counts=validation.groupby('affected_method').outcome.apply(lambda x:x.eq('PASS').sum())
ax.barh(counts.index,counts.values,color=[COLORS.get(m,'#555555') for m in counts.index])
for y,v in enumerate(counts.values):ax.text(v+.3,y,str(v),va='center')
ax.set_xlim(0,max(counts.values)*1.15);ax.set_xlabel('Completed checks with PASS');ax.set_facecolor('white')
fig.tight_layout();fig.savefig(OUT/'figures/correctness/check_coverage.png',dpi=160);plt.close(fig)
nbpath=ROOT/'notebooks_paper/01_correctness.ipynb';nb=nbformat.read(nbpath,as_version=4)
if not any('![Finite validation coverage]' in c.source for c in nb.cells):
    nb.cells.append(nbformat.v4.new_markdown_cell('![Finite validation coverage](../output/paper_v1/figures/correctness/check_coverage.png)\n\nCounts describe check coverage, not a statistical correctness probability.'))
nbformat.write(nb,nbpath)
changed=[]
for filename in ['existing_inventory.csv','protected_results_snapshot.csv']:
    for r in pd.read_csv(OUT/filename).itertuples():
        p=ROOT/r.path
        if not p.exists() or sha(p)!=r.sha256:changed.append(r.path)
(OUT/'preservation_audit.json').write_text(json.dumps(dict(changed=sorted(set(changed)),inventory_entries=219,protected_result_files=1054),indent=2))
assert not changed,changed
saved=[]
for p in sorted((ROOT/'notebooks_paper').glob('*.py')):saved.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
pd.DataFrame(saved).to_csv(OUT/'stage1_modules.csv',index=False)
optimal=allrows.status.eq('optimal');missing=int((~optimal).sum());peak=allrows.peak_process_rss_mib.max()
prior=pd.read_csv(OUT/'prior_experiments.csv')
mask=prior.experiment.eq('Stage 1 pilot')
prior.loc[mask,'coverage']='completed' if missing==0 else 'partial'
prior.loc[mask,'optimal']=int(optimal.sum());prior.loc[mask,'nonoptimal']=missing
prior.loc[mask,'note']='All 32 first attempts recorded; see Stage 1 status for resource failures'
prior.to_csv(OUT/'prior_experiments.csv',index=False)
total=allrows.process_wall_seconds.sum();solvertime=allrows.solve_seconds.sum()
status=f'''# Stage 1 — completed execution, partial pilot coverage

## Completed
- New notebooks 00–03 and shared adapters under notebooks_paper. Existing algorithm implementations reused unchanged.
- Frozen inputs: two protein pairs at existing w=16 (22 vertices per curve), two hurricane pairs (25 per curve); no further sampling or alignment.
- 8 pair/alpha configurations, 4 methods each, 32 attempted configurations. alpha=0.5,1; separate mean-edge fidelity thresholds; exact unrounded input discrete coupling budget.
- Correctness: {validation.outcome.eq('PASS').sum()}/{len(validation)} checks PASS; {validation.outcome.eq('FAIL').sum()} failures; {validation.outcome.eq('INCOMPLETE').sum()} incomplete. Exhaustive small-input checks force certificates off. Finite checks do not prove the draft theorem or eliminate shared distance-oracle errors.
- Pilot optimal outputs: {optimal.sum()}/32; missing: {missing}/32. See [pilot summary](pilot_summary.csv), [full results](pilot_results.csv), [parameters](pilot_parameters.csv), and [validation](validation.csv).
- Preservation audit: no changes to the 219 inventoried old files or 1,054 old output artifacts. Existing Jupyter kernel was not attached to or stopped.

## Failures and limitations
- Protein attempts used a launcher-only Windows memory monitor. Memory values are invalid and excluded; no solution was returned. Original attempts remain in logs/protein_initial_monitoring_raw.csv. They were not repeated. Two warm-ups (one protein, one hurricane) hit the soft resource cap before the pilot solver started and are classified startup_resource_limit after reviewing their tracebacks; these are not algorithmic correctness failures.
- Corrected hurricane monitoring sums RSS over the launcher and descendants and terminates the worker tree on a cap. Polling limits can overshoot while the OS delays scheduling. The pre-existing kernel is excluded.
- {allrows.status.value_counts().to_dict()}
- Resource exits do not establish mathematical infeasibility. No failed result is used in a quality comparison. Certificate routes are distinct from graph execution.
- Both protein pairs reuse 1o7j.a. Existing long-edge flags are reported as supplied, not newly verified missing-residue diagnoses. The dataset README's older index-wise coupling convention is not used: this pilot computes discrete Fréchet coupling with repeated indices.
- No full-resolution protein experiment, Contact F1/topology claim, phase timing, exact algorithm peak-memory measurement, or population-level conclusion.

## Runtime and resources
- One attempt per configuration; no timing repetitions or cached repetitions.
- Sum of worker wall times: {total:.3f} s. Recorded completed solver/soft-limit times: {solvertime:.3f} s (excludes unavailable killed-worker solve times). Notebook timings are in logs/notebook_execution.jsonl.
- Maximum valid sampled worker-tree RSS: {peak:.2f} MiB. Protein memory unavailable. Limits: 512 MiB, 20 s after warm-up, 90 s startup; joint graph soft limit 8 s, 30,000 states, 500,000 transitions, 200,000 chain edges.
- First attempts used below-normal priority; corrected remaining workers use ordinary priority. Thus these diagnostics are unsuitable for method runtime rankings.

## Planned, not executed
- 04_graph_cost.ipynb and 05_paper_summary.ipynb are not created or run in Stage 1. No paper conclusions drafted.
- Stage 2 proposal: first resolve the {missing} missing configurations on the same frozen inputs, with corrected monitoring and a proposed 60 s execution cap: at most {missing*60/60:.0f} worker-minutes plus startup. Do not enlarge the cohort until coverage is complete. This is a proposal only and requires a further request.
- Future runtime study: 10 measured repetitions after warm-up, increased only after observing variability; separate graph phases and certificate conditions. Different alpha settings are not repetitions.
'''
(OUT/'status.md').write_text(status,encoding='utf-8')
print(wide.to_string());print(status)
