"""Check expanded scope, cache provenance, and agreement between both notebooks."""

# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path[:0] = [str(_Path(__file__).resolve().parents[2] / 'src'), str(_Path(__file__).resolve().parents[2])]

from pathlib import Path
import json,ast,shutil
import nbformat
import pandas as pd
import numpy as np
root=Path(__file__).resolve().parents[2]
base=root/'output/hurricane_cps_expanded';aux=root/'output/hurricane_auxiliary_expanded'
r=pd.read_csv(base/'comparison.csv');old=pd.read_csv(root/'output/hurricane_cps/comparison.csv')
manifest=json.loads((base/'manifest.json').read_text())
assert manifest['selected_pairs']==32 and manifest['pairs_per_basin']==16
assert manifest['budgets']['seconds']==120 and manifest['budgets']['max_transitions']==20_000_000
assert set(old.pair)<=set(r.pair)
assert len(r)==192 and r.pair.nunique()==32
assert (base/'comparison.csv').read_bytes()==(aux/'comparison.csv').read_bytes()
fresh=r[~r.reused_verified_result.astype(bool)]
assert fresh.soft_time_limit_seconds.eq(120).all() and fresh.transition_limit.eq(20_000_000).all()
prior=old[old.status.str.startswith('optimal')].set_index(['case','method'])
current=r.set_index(['case','method'])
for key,row in prior.iterrows():
    assert current.loc[key,'reused_verified_result']
    assert current.loc[key,'k']==row.k
    assert np.isclose(current.loc[key,'dDF_simplified'],row.dDF_simplified)
previous_limits=old[old.status.eq('resource_limit')][['case','method','pair','alpha']]
retry=previous_limits.merge(r,on=['case','method','pair','alpha'])
assert (~retry.reused_verified_result.astype(bool)).all()
retry[['pair','alpha','method','status','k','seconds','reason']].to_csv(base/'retried_cases.csv',index=False)
shutil.copy2(base/'retried_cases.csv',aux/'retried_cases.csv')
for name in ['notebooks/legacy/1.analyze_CPS.ipynb','notebooks/legacy/0.analyze_hurricans.ipynb']:
    path=root/name;nb=nbformat.read(path,4)
    config=next(c.source for c in nb.cells if c.cell_type=='code' and 'PAIRS_PER_BASIN=' in c.source)
    assert 'PAIRS_PER_BASIN=16' in config and 'seconds=120.' in config and 'max_transitions=20_000_000' in config
    for c in nb.cells:
        c.source=c.source.replace('verified output/hurricane_cps results','verified output/hurricane_cps_expanded results')
        if c.cell_type=='code' and c.source.startswith('columns='):
            c.source=c.source.replace("'seconds','status'","'seconds','reused_verified_result','status'")
            columns=ast.literal_eval(ast.parse(c.source).body[0].value)
            table_output=next(o for o in c.outputs if 'text/html' in o.get('data',{}))
            table_output.data['text/html']=r[columns].to_html(index=False,float_format=lambda x:f'{x:.3f}')
        if c.cell_type=='markdown' and c.source.startswith('## 6.'):
            c.source+='\n\n**Timing provenance:** `reused_verified_result=True` identifies cached exact solutions. Their seconds are from the original solve; they were not timed again. New and previously unfinished cases use the increased limits. See `retried_cases.csv` for the three original NA cases.'
    if name=='notebooks/legacy/0.analyze_hurricans.ipynb':
        # Copy refreshed shared output links to this notebook's own artifact folder.
        nb=nbformat.reads(nbformat.writes(nb).replace('output/hurricane_cps_expanded/','output/hurricane_auxiliary_expanded/'),4)
    nbformat.validate(nb);nbformat.write(nb,path)
am=aux/'auxiliary_manifest.json'
a=json.loads(am.read_text());a['source']='verified output/hurricane_cps_expanded results; same inputs and settings';am.write_text(json.dumps(a,indent=2))
report=dict(selected_pairs=32,available_pairs=2222,method_results=len(r),
    reused_optimal=int(r.reused_verified_result.sum()),fresh_runs=int((~r.reused_verified_result).sum()),
    statuses=r.status.value_counts().to_dict(),previous_limit_retries=retry[['pair','alpha','method','status']].to_dict('records'))
(base/'expansion_validation.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
