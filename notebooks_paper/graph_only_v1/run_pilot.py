"""Only the frozen 32-configuration Stage 1 cohort. No cached timing repetitions."""
from graph_core import ROOT,OUT,METHODS,PROTOCOL,sha
from pathlib import Path
import sys,json,time,subprocess,psutil
import pandas as pd

def run(domain):
    target=OUT/f'{domain}_results.csv'
    if target.exists():raise RuntimeError('Results already exist; read them without rerunning optimization.')
    gate=pd.read_csv(OUT/'validation.csv');assert gate.outcome.eq('PASS').all(),'Correctness gate failed'
    params=pd.read_csv(OUT/'pilot_parameters.csv');protocol=json.loads((OUT/'protocol.json').read_text())
    assert sha(OUT/'pilot_parameters.csv')==protocol['parameter_sha256']
    for name,digest in protocol['algorithms'].items():assert sha(ROOT/name)==digest
    for r in pd.read_csv(OUT/'input_manifest.csv').itertuples():assert sha(ROOT/r.input_file)==r.input_sha256
    rows=[];paths={}
    for p in params[params.domain.eq(domain)].to_dict('records'):
        for method in METHODS:
            stem=p['case'].replace(':','_')+'_'+method.replace(' ','_')
            task=OUT/'workers'/f'{stem}.json';ready=task.with_suffix('.ready');result=task.with_suffix('.result.json');progress=task.with_suffix('.progress.json')
            assert not task.exists(),'No repeated attempts in this run'
            task.write_text(json.dumps(dict(param=p,method=method,ready=str(ready),result_path=str(result),progress=str(progress))))
            row=dict(**p,method=method,certificate_enabled=False,status='not_run',route='startup',phase='startup',
                     configuration_graph_entered=False,chain_graph_calls=0,chain_graphs_completed=0,repetition=1,cached=False,
                     hard_seconds=PROTOCOL['hard_seconds'],memory_limit_mib=PROTOCOL['memory_mib'],failure_reason='')
            start=time.perf_counter();peak=0;began=None;limit=None;observed={}
            with (OUT/'logs'/f'{stem}.log').open('w') as log:
                process=subprocess.Popen([sys.executable,str(Path(__file__).with_name('pilot_worker.py')),str(task)],cwd=ROOT,
                                         stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
                while process.poll() is None:
                    try:
                        tree=psutil.Process(process.pid).children(recursive=True)+[psutil.Process(process.pid)]
                        rss=0
                        for child in tree:
                            try:rss+=child.memory_info().rss
                            except psutil.NoSuchProcess:pass
                        peak=max(peak,rss)
                    except psutil.NoSuchProcess:pass
                    try:
                        if ready.exists() and began is None:began=json.loads(ready.read_text())['started_unix']
                        if progress.exists():observed=json.loads(progress.read_text())
                    except (ValueError,OSError):pass  # a concurrent short write is retried
                    if peak>PROTOCOL['memory_mib']*1024**2:limit='memory_limit'
                    elif began is not None and time.time()-began>PROTOCOL['hard_seconds']:limit='time_limit'
                    elif began is None and time.perf_counter()-start>PROTOCOL['startup_seconds']:limit='startup_time_limit'
                    if limit:
                        for child in reversed(tree):
                            try:child.kill()
                            except psutil.NoSuchProcess:pass
                        process.wait();break
                    time.sleep(.05)
            row.update(process_wall_seconds=time.perf_counter()-start,peak_process_rss_mib=peak/1024**2,measured_execution_count=int(began is not None),**observed)
            if limit:row.update(status=limit,raw_status=limit,route=row['phase'],failure_reason='Worker-tree watchdog resource exit; not infeasibility')
            elif result.exists():row.update(json.loads(result.read_text()))
            else:row.update(status='software_error',raw_status='software_error',route=row['phase'],failure_reason=f'Worker exit {process.returncode}; see log')
            if row['status']=='optimal' and method.startswith('CPS'):
                assert row['raw_status']=='optimal_configuration_graph' and row['configuration_graph_entered']
            paths[p['case']+'|'+method]={k:row.pop(k) for k in ['indices_A','indices_B','coupling'] if k in row}
            rows.append(row);pd.DataFrame(rows).to_csv(target,index=False)
            (OUT/f'{domain}_paths.json').write_text(json.dumps(paths,indent=2))
            print(p['pair'],p['alpha'],method,row['status'],row['route'],flush=True)
    frame=pd.DataFrame(rows)
    for i,r in frame.iterrows():
        baseline='Independent continuous' if r.method=='CPS-2F' else 'Independent discrete' if r.method=='CPS-3F' else r.method
        b=frame[(frame.case==r.case)&(frame.method==baseline)].iloc[0]
        frame.loc[i,'baseline_method']=baseline;frame.loc[i,'baseline_missing']=b.status!='optimal'
        if r.status==b.status=='optimal':frame.loc[i,'vertex_cost_vs_matching_baseline']=r.k-b.k
    frame.to_csv(target,index=False)
    return frame

if __name__=='__main__':
    assert sys.argv[1:] in [['protein'],['hurricane'],['protein','hurricane']]
    for domain in sys.argv[1:]:run(domain)
