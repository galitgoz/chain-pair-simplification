from pathlib import Path
from protein_data import prepare
from curve_algorithms import shortest_independent,continuous_distance,discrete_frechet,local_shortcut_path
import numpy as np
import pandas as pd
import json,time

def run(chains,thresholds=(2.,4.,8.),directory='output/independent'):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    rows=[]; paths={}
    for key,chain in chains.items():
        for fragment in chain['fragments']:
            P=fragment['xyz']
            for delta in thresholds:
                t=time.perf_counter(); indices=shortest_independent(P,delta); elapsed=time.perf_counter()-t
                lo,hi=continuous_distance(P,P[indices])
                record=dict(chain=key,fragment=fragment['id'],delta_A=delta,n=len(P),k=len(indices),
                    retained_fraction=len(indices)/len(P),continuous_F_lower_A=lo,continuous_F_upper_A=hi,
                    discrete_F_A=float(discrete_frechet(P,P[indices])),runtime_s=elapsed,
                    local_k=len(local_shortcut_path(P,delta)))
                rows.append(record)
                label=f"{fragment['id']}@{delta:g}"
                paths[label]=dict(fragment_indices=indices.tolist(),chain_indices=fragment['chain_indices'][indices].tolist())
                print(f"{label}: {len(P)} -> {len(indices)} ({elapsed:.2f}s)",flush=True)
                pd.DataFrame(rows).to_csv(directory/'results.csv',index=False)
                (directory/'paths.json').write_text(json.dumps(paths,indent=2))
    return pd.DataFrame(rows),paths

if __name__=='__main__':
    chains,*_=prepare('.')
    run(chains)
