from pathlib import Path
from cps_notebook_data import load_inputs
from cps_notebook_run import run_comparison
data=load_inputs(Path('.'),ws=(16,))
print(data['inputs'][['B','nA','nB','dropped_decimated_vertices_per_chain']].to_string(index=False),flush=True)
frame,paths=run_comparison(data,'output/cps_papers',budgets=dict(max_states=300000,max_transitions=15000000,seconds=45.))
print(frame.groupby(['method','status']).size().to_string(),flush=True)
