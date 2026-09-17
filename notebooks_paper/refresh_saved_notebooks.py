"""Prepare saved-results rendering; this script never runs optimization."""
from paper_core import ROOT
import nbformat
for name in ['02_proteins.ipynb','03_hurricanes.ipynb']:
    p=ROOT/'notebooks_paper'/name;nb=nbformat.read(p,as_version=4)
    for c in nb.cells:
        if c.cell_type=='code' and 'from paper_core import run_pilot' in c.source:
            c.source="""# Read this run's saved results; this is not another timing repetition.
results = pd.read_csv(OUT / f'{DOMAIN}_results.csv')
display(compact_results(results))
display(results.groupby(['status','route']).size().rename('configurations').to_frame())"""
    nb.cells.append(nbformat.v4.new_markdown_cell('**Run status:** Read [status.md](../output/paper_v1/status.md) for verified coverage and limitations. Small horizontal offsets in paired plots make equal method values visible; offsets do not change any measured value. All original attempts and logs are retained.'))
    nbformat.write(nb,p)
print('Domain notebooks configured to read saved results only.')
