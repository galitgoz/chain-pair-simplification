"""Preserve failed attempts; do not repeat pilot optimizations."""
from paper_core import OUT,ROOT
import pandas as pd,nbformat,json
p=OUT/'protein_results.csv';d=pd.read_csv(p)
d.to_csv(OUT/'logs/protein_initial_monitoring_raw.csv',index=False)
d['resource_monitoring_valid']=False
d['peak_process_rss_mib']=float('nan')
d['failure_reason']=d.failure_reason.fillna('')+'; launcher-only monitor: child memory unmeasured. No solution recovered; not repeated.'
d.to_csv(p,index=False)
nb=nbformat.read(ROOT/'notebooks_paper/02_proteins.ipynb',as_version=4)
for c in nb.cells:
    if c.cell_type=='code' and 'results = run_pilot(DOMAIN)' in c.source:
        c.source=c.source.replace('results = run_pilot(DOMAIN)',"results = pd.read_csv(OUT / 'protein_results.csv')  # preserved first attempts, no repeat")
nb.cells.insert(5,nbformat.v4.new_markdown_cell('**Execution limitation:** The 16 protein configurations were attempted once. The original Windows monitor measured only the Python launcher, so its memory values are invalid and have been removed. No solution was returned before the resource exits. The attempts and original logs are preserved; no extra timing repetition was launched. The worker monitor was corrected to account for the full child-process tree before the hurricane pilot. Protein objectives remain explicitly missing.'))
nbformat.write(nb,ROOT/'notebooks_paper/02_proteins.ipynb')
print(d.status.value_counts().to_dict())
