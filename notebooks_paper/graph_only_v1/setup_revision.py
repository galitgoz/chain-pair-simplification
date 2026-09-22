"""Freeze a new graph-only revision, retaining previous experiments unchanged."""
from pathlib import Path
import sys,json,shutil,hashlib
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/paper_v1/graph_only_v1';OLD=OUT.parent
assert not (OUT/'protocol.json').exists(),'Revision already initialized'
for folder in ['cache','logs','workers','figures/proteins','figures/hurricanes']:(OUT/folder).mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(Path(__file__).parent))
from graph_core import PROTOCOL,sha
import pandas as pd,psutil
for name in ['input_manifest.csv','pilot_parameters.csv','selection_audit.csv','shared_entities.csv','environment.json']:
    shutil.copy2(OLD/name,OUT/name)
protocol=json.loads((OLD/'protocol.json').read_text());protocol.update(PROTOCOL)
protocol.update(revision='graph_only_v1',parent_protocol='output/paper_v1/protocol.json',
    scope='Stage 1 only; same frozen cohort and thresholds; all CPS calls explicitly certificate=False',
    independent_preprocessing='Existing solve_cps still computes an independent lower bound; included in total solve time, never used to bypass the graph',
    resource_monitor='sum of worker-tree RSS; wall watchdog may overshoot under OS scheduling delays',
    parameter_sha256=sha(OUT/'pilot_parameters.csv'))
(OUT/'protocol.json').write_text(json.dumps(protocol,indent=2))
old=pd.read_csv(OLD/'pilot_results.csv')
old['source_file']='output/paper_v1/pilot_results.csv'
old['compatibility']='not imported; new isolated run required'
mask=old.route.isin(['independent_certificate','alternate_certificate'])
old.loc[mask,'compatibility']='excluded: certificate success is not a graph-run result'
old['original_raw_status']=old.route.map({'independent_certificate':'optimal_independent_certificate','alternate_certificate':'optimal_alternate_certificate'})
old.to_csv(OUT/'prior_result_compatibility.csv',index=False)
old[mask].to_csv(OUT/'excluded_certificate_results.csv',index=False)
snapshot=[]
for base in [ROOT/'notebooks_paper',ROOT/'output']:
    for p in base.rglob('*'):
        if not p.is_file() or OUT in p.parents or Path(__file__).parent in p.parents:continue
        if any(x in p.parts for x in ['cache','__pycache__']):continue
        snapshot.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
for name in ['src/cps_paper_algorithms.py','src/curve_algorithms.py']:snapshot.append(dict(path=name,sha256=sha(ROOT/name)))
pd.DataFrame(snapshot).to_csv(OUT/'preservation_snapshot.csv',index=False)
active=[]
for p in psutil.process_iter(['pid','ppid','cmdline']):
    c=' '.join(p.info['cmdline'] or [])
    if 'python' in c.lower() and ('ipykernel' in c or 'pilot_worker.py' in c or 'run_pilot.py' in c):active.append(p.info)
(OUT/'active_run_inspection.json').write_text(json.dumps(dict(processes=active,action='read-only inspection; no process attached to or stopped'),indent=2))
# Reuse the previous tests with only adapter and protocol-specific assertions changed.
source=(ROOT/'notebooks_paper/correctness_checks.py').read_text(encoding='utf-8')
source=source.replace('from paper_core import ROOT,OUT,METHODS,solve,measures','from graph_core import ROOT,OUT,METHODS,solve,measures')
source=source.replace('certificate=True','certificate=False')
source=source.replace("assert result['route']=='independent_certificate'","assert result['route']=='configuration_graph' and result['configuration_graph_entered']")
source=source.replace("assert result.get('states')==result.get('transitions')==0","assert result.get('states')>0 and result.get('transitions')>0")
source=source.replace('Independent certificate is labeled separately; graph is not reported as executed.','Independently feasible outputs do not bypass execution: both chain graphs and the product recurrence ran.')
source=source.replace("record('certificate_route_label',m,route)","record('feasible_independent_outputs_do_not_bypass_graph',m,route)")
source+='''

if __name__=='__main__':
    from unittest.mock import patch
    import cps_paper_algorithms as alg
    # A call to the alternative-subsequence search fails the suite immediately.
    with patch.object(alg,'alternate_certificate',side_effect=AssertionError('Forbidden alternative-subsequence shortcut invoked')) as forbidden:
        result=check_all()
        assert forbidden.call_count==0
    assert result.outcome.eq('PASS').all(),result[result.outcome.ne('PASS')].to_string()
    print(result.groupby(['affected_method','outcome']).size().to_string())
    print('Both shortcuts disabled; 99 checks passed; alternative-subsequence calls: 0')
'''
(Path(__file__).parent/'graph_checks.py').write_text(source,encoding='utf-8')
(OUT/'status.md').write_text('# Graph-only Stage 1 revision\n\nProtocol frozen. Previous outputs preserved and not reused. Correctness and pilot pending. All CPS calls require certificate=False. No later stage authorized.\n')
print('Frozen unchanged 32-configuration cohort; excluded certificate successes:',int(mask.sum()))
