"""Validate saved reference witnesses; no optimization or certificate search."""
from common import ROOT,OUT,validate,atomic,sha
import pandas as pd,numpy as np,json,time
configs=json.loads((OUT/'diagnostics.json').read_text());records=[]
def witness(source,domain,config,method,status):
    table=pd.read_csv(source/'pilot_results.csv');row=table[table.case.eq(config['case']) & table.method.eq(method)].iloc[0]
    assert row.status=='optimal'
    for key in ['delta1','delta2','delta3']:assert row[key]==config[key]
    paths=json.loads((source/f'{domain}_paths.json').read_text())[config['case']+'|'+method]
    raw=dict(**paths,status=status,k=int(row.k),kA=int(row.kA),kB=int(row.kB))
    data=np.load(ROOT/config['input_file']);start=time.perf_counter()
    validation=validate(data['A'],data['B'],raw,{**config,'method':method})
    evidence=dict(source_table=str((source/'pilot_results.csv').relative_to(ROOT)),source_sha256=sha(source/'pilot_results.csv'),source_paths_sha256=sha(source/f'{domain}_paths.json'),source_method=method,raw_result=raw,validation=validation,validation_seconds=time.perf_counter()-start)
    return evidence
for c in configs:
    if c['diagnostic']=='A':
        certificate=witness(OUT.parent,'hurricane',c,'CPS-2F','optimal_alternate_certificate')
        baseline=witness(OUT.parent/'graph_only_v1','hurricane',c,'Independent continuous','optimal_independent_baseline')
        assert certificate['raw_result']['k']==baseline['raw_result']['k']==5
        record=dict(diagnostic='A',reference_k=5,evidence_kind='saved certificate witness and independently optimal continuous lower bound; NOT a graph result',certificate=certificate,baseline=baseline,reference_verified=True)
    elif c['diagnostic']=='B':
        baseline=witness(OUT.parent/'graph_only_v1','protein',c,'Independent continuous','optimal_independent_baseline')
        assert baseline['raw_result']['k']==5 and baseline['validation']['coupling_feasible']
        record=dict(diagnostic='B',reference_k=5,evidence_kind='saved independently optimal continuous pair also satisfies delta3; NOT a new graph run',baseline=baseline,reference_verified=True)
    else:record=dict(diagnostic='C',reference_k=None,reference_verified=False,evidence_kind='No verified identical-configuration reference in the frozen-protocol prior runs')
    records.append(record);atomic(OUT/'references'/f"{c['diagnostic']}.json",record)
atomic(OUT/'references.json',records)
print('Reference A: verified k=5 witness and lower bound; B: verified k=5 independent feasible optimum; C: unavailable.')
