"""Report only saved Stage 1 graph results; no optimization."""
from graph_core import ROOT,OUT,METHODS,sha
from pathlib import Path
import ast,json,pandas as pd,numpy as np
HERE=Path(__file__).parent
frames=[pd.read_csv(OUT/f'{d}_results.csv') for d in ['protein','hurricane']]
r=pd.concat(frames,ignore_index=True)
assert len(r)==32 and not r.duplicated(['case','method']).any()
assert not r.certificate_enabled.any()
assert r.groupby('case')[['delta1','delta2','delta3']].nunique().eq(1).all().all()
joint=r[r.method.str.startswith('CPS')];success=joint[joint.status.eq('optimal')]
assert success.raw_status.eq('optimal_configuration_graph').all()
assert success.route.eq('configuration_graph').all() and success.configuration_graph_entered.all()
assert success.chain_graphs_completed.eq(2).all()
assert success.states_discovered.ge(1).all()
assert (success.coupling_distance<=success.delta3+1e-9).all()
assert not r.route.astype(str).str.contains('certificate').any()
r.to_csv(OUT/'pilot_results.csv',index=False)
summary=r[['pair','alpha','method','status','k']].copy()
summary['entry']=[str(int(k)) if s=='optimal' else 'missing: '+s for k,s in zip(summary.k,summary.status)]
wide=summary.pivot(index=['pair','alpha'],columns='method',values='entry').reindex(columns=METHODS)
wide.to_csv(OUT/'pilot_summary.csv')
calls=[]
for p in HERE.glob('*.py'):
    for node in ast.walk(ast.parse(p.read_text(encoding='utf-8'))):
        if not isinstance(node,ast.Call):continue
        name=node.func.id if isinstance(node.func,ast.Name) else node.func.attr if isinstance(node.func,ast.Attribute) else ''
        if name not in ['solve','solve_cps']:continue
        cert=[k.value for k in node.keywords if k.arg=='certificate']
        explicit=len(cert)==1 and isinstance(cert[0],ast.Constant) and cert[0].value is False
        assert explicit,(p,node.lineno)
        calls.append(dict(file=str(p.relative_to(ROOT)),line=node.lineno,target=name,explicit_certificate_false=explicit))
pd.DataFrame(calls).to_csv(OUT/'callsite_audit.csv',index=False)
changed=[]
for old in pd.read_csv(OUT/'preservation_snapshot.csv').itertuples():
    p=ROOT/old.path
    if not p.exists() or sha(p)!=old.sha256:changed.append(old.path)
(OUT/'preservation_audit.json').write_text(json.dumps(dict(changed=changed,files_checked=len(pd.read_csv(OUT/'preservation_snapshot.csv'))),indent=2))
assert not changed,changed
v=pd.read_csv(OUT/'validation.csv');assert len(v)==100 and v.outcome.eq('PASS').all()
oldcheck=pd.read_csv(OUT/'logs/initial_observer_validation.csv')
initialfail=int(oldcheck.outcome.ne('PASS').sum())
missing=int(r.status.ne('optimal').sum());total=r.process_wall_seconds.sum();peak=r.peak_process_rss_mib.max()
joint.groupby(['status','route']).size().rename('rows').to_csv(OUT/'graph_route_summary.csv')
text=f'''# Graph-only Stage 1 — final report

## Protocol correction
All CPS calls, including warm-up and correctness checks, explicitly pass certificate=False. The adapter rejects True. The independent early return and alternative-subsequence search are disabled; no shortcut functions were deleted. Independent continuous and discrete baselines are separate runs. Existing lower-bound preprocessing remains inside solve_cps and is included in total time, never used to bypass graph execution.

The same frozen cohort and numerical thresholds were reused: two protein pairs at existing w=16, two hurricane pairs, alpha=0.5 and 1, four methods, 32 attempts. No larger or later-stage experiment was started.

## Verification
- 100/100 final small-input checks PASS, including exhaustive enumeration, independently feasible examples, and the 3-vertex example eligible for alternative certification.
- The guarded alternative-subsequence function was called zero times.
- An initial observer revision failed {initialfail} assertions because it counted baseline single-curve graphs as joint construction. The observer was corrected before the pilot. The initial check log is preserved; this was an instrumentation error, not an optimization mismatch.
- The observer now counts only direct chain-graph calls from solve_cps and marks product-graph entry at its first heap pop. Successful joint rows additionally require raw status optimal_configuration_graph. Baseline graph instrumentation is not a performance metric.
- {len(success)}/{len(joint)} CPS attempts returned verified optimal graph solutions. All {len(success)} successful CPS rows actually entered the product graph and completed both chain graphs. No certificate success is included.
- {r.status.eq('optimal').sum()}/32 total verified solutions, {missing}/32 missing. Status counts: {r.status.value_counts().to_dict()}.
- Joint route outcomes: {joint.groupby(['status','route']).size().to_dict()}.

## Preservation and exclusion
Six previous successful CPS certificate rows are listed in excluded_certificate_results.csv. All 32 previous pilot rows are tracked in prior_result_compatibility.csv and none is reused as a new timing measurement. Original notebooks, results, and root algorithm modules are unchanged; preservation audit passed for {len(pd.read_csv(OUT/'preservation_snapshot.csv'))} files. Read-only process inspection found the pre-existing Jupyter kernel, which was neither attached to nor stopped.

## Call sites and thesis text
See callsite_audit.csv for every executable solve/solve_cps call in the revised modules. The new call path replaces the previous certificate-default pilot adapter: graph_core.solve -> solve_cps(certificate=False), pilot_worker measured call, graph_core warm-up, graph_checks tests, and shortcut_probe all use explicit False. Prior experiment modules are preserved as historical software.

The methods/protocol insert is notebooks_paper/graph_only_v1/thesis_methods_protocol.md. No editable thesis draft was found, so reference PDFs were preserved. The alternative-subsequence search is not a studied thesis method, contribution, or local replacement algorithm. No paper conclusions are drafted.

## Resources and limitations
One attempt per configuration, diagnostic timings only. Worker wall-time sum {total:.3f} s; maximum sampled worker-tree RSS {peak:.2f} MiB. Limits: 20 s execution watchdog after warm-up, 90 s startup, 512 MiB RSS; graph caps 8 s soft time, 30,000 states, 500,000 product transitions and 200,000 chain edges. Watchdog polling can overshoot under OS scheduling delays. Warm-up, preprocessing, graph construction, and graph resource exits are distinguished by route/phase; missing outputs are not infeasible. No per-phase timing benchmark is claimed.

Protein curves are already subsampled, share 1o7j.a, and retain supplied gap bridges; 1d9q.d has the existing documented count discrepancy (325 supplied versus paper length 297). Hurricane inputs retain the common basin frame and lack timestamps. These remain pipeline examples, not evidence of topological preservation or broad performance conclusions.

## Outputs and next step
New notebooks 00–03 are under notebooks_paper/graph_only_v1. Tables, figures, logs, exclusions and validation are under output/paper_v1/graph_only_v1. See pilot_summary.csv, pilot_results.csv, graph_route_summary.csv and validation.csv.

Stage 1 stops here. Any retry of the {missing} missing configurations, change of resource caps, or larger experiment requires a further request. Notebooks 04/05 and later stages remain unexecuted.
'''
(OUT/'status.md').write_text(text,encoding='utf-8')
print(wide.to_string());print(text)
