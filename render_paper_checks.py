from pathlib import Path
import pymupdf
out=Path('output/paper_checks'); out.mkdir(parents=True,exist_ok=True)
for filename,pages,prefix in [('papers/global curve simplification.pdf',[12],'gcs'),('papers/chain_pair_simplification.pdf',[17,18],'cps')]:
    with pymupdf.open(filename) as doc:
        for page in pages:
            doc[page].get_pixmap(matrix=pymupdf.Matrix(1.4,1.4)).save(str(out/f'{prefix}_{page+1}.png'))
