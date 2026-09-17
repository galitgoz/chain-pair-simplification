from paper_core import ROOT,OUT
import nbformat,pandas as pd
p=ROOT/'notebooks_paper/00_data_and_protocol.ipynb';nb=nbformat.read(p,as_version=4)
for c in nb.cells:
    if c.cell_type=='code' and "manifest = pd.read_csv" in c.source:
        c.source=c.source.replace("display(manifest[['pair','subsampling','gaps_A','gaps_B','alignment']])", "brief = manifest[['pair','subsampling','gaps_A','gaps_B','alignment']].copy()\nbrief['subsampling'] = brief['subsampling'].str.split(';').str[0]\ndisplay(brief)\nselection = pd.read_csv(OUT / 'selection_audit.csv')\ndisplay(selection.groupby('domain').agg(eligible_pairs=('pair','size'), selected_pairs=('selected','sum')))")
nb.cells.append(nbformat.v4.new_markdown_cell('**Upstream protein coverage caveat:** The existing identifier audit records 325 supplied Cα coordinates for 1d9q.d versus a reported paper length of 297, and three contiguous fragments before representing the stored chain as a single curve. This pilot retains the supplied gap-bridging edges and applies no new fragmentation. That audit is provenance metadata, not a new structural verification in Stage 1. “Entire supplied curve” does not mean a complete biological polymer.'))
nbformat.write(nb,p)
status=OUT/'status.md'
status.write_text(status.read_text(encoding='utf-8')+'\n## Additional source caveat\nThe pre-existing identifier audit records 1d9q.d as 325 supplied coordinates versus paper length 297, with three contiguous fragments. The supplied single curve retains gap-bridging edges. No new fragmentation or biological completeness claim is made.\n',encoding='utf-8')
print('Documented upstream protein count and gaps; no input coordinates or thresholds changed.')
