
# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path[:0] = [str(_Path(__file__).resolve().parents[2] / 'src'), str(_Path(__file__).resolve().parents[2])]

from pathlib import Path
import pymupdf
out=Path('output/cps_papers/source_checks');out.mkdir(parents=True,exist_ok=True)
for name,pages,prefix in [('M_Sc__Thesis___Galit_Gozoltzani.pdf',[31,32],'thesis'),('global curve simplification.pdf',[12],'gcs'),('chain_pair_simplification.pdf',[8,17,18],'fan')]:
    doc=pymupdf.open(Path('papers')/name)
    for p in pages:doc[p].get_pixmap(matrix=pymupdf.Matrix(1.25,1.25)).save(str(out/f'{prefix}_{p+1}.png'))
print('Rendered source algorithm/table pages')
