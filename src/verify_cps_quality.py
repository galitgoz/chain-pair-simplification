import numpy as np
from cps_quality import contact_metrics
from cps_paper_algorithms import extended_curve

P=np.array([[0,0,0],[2,0,0],[2,2,0],[0,.1,0]],float)
seq=np.array([0,5,10,15])
m=contact_metrics(P,np.arange(4),seq,.2,5)
assert m['contact_TP']==1 and m['contact_FP']==m['contact_FN']==0
assert m['contact_F1']==m['contact_precision']==m['contact_recall']==1
m=contact_metrics(P,[0,3],seq,.2,5)
assert (m['contact_TP'],m['contact_FP'],m['contact_FN'])==(1,5,0)
assert np.isclose(m['contact_F1'],2/7)
assert np.isnan(contact_metrics(P,np.arange(4),seq,0,5)['contact_F1'])
line=np.array([[0,0,0],[2,0,0],[4,0,0]],float)
X,t=extended_curve(line,1.)
assert len(X)-len(line)==2 and np.allclose(t,[0,.5,1,1.5,2])
assert len(extended_curve(line,0)[0])==len(line)
print('Contact identity, false contacts, undefined F1, and auxiliary deduplication checks passed.')
