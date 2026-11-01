from paths import OUTPUT, ROOT, RESEARCH
from pathlib import Path
import json
import numpy as np

root=OUTPUT;td=np.load(root/'tetra-v1/data.npz');tl=np.load(root/'tetra-long/data.npz');ad=np.load(root/'anderson-v2/data.npz');al=np.load(root/'anderson-long/data.npz')
facts={'tetra':[],'anderson':[]}
for i in [421,565,705,117]:
    facts['tetra'].append({'id':i,'seed_center':td['centers'][i,3].tolist(),'seed_axis':td['axes'][i,3].tolist(),'sigma':float(td['sigma'][i]),'gamma':float(td['gamma'][i]),'start':td['starts'][i].tolist(),'absorbing_radius_threshold':float(12*td['sigma'][i]*np.exp(-.5)/td['gamma'][i])})
for i in [54,449,2268]:
    p=ad['parameters'][i];z=ad['starts'][i,0];before=float(z)
    if p[3]==2:
        for _ in range(10000):z=p[1]+p[0]*np.tanh(z)
    facts['anderson'].append({'id':i,'parameters_k_c_lambda_kind':p.tolist(),'cold_history':ad['starts'][i].tolist(),'plain_iteration_start':before,'plain_iteration_after_10000':float(z) if p[3]==2 else None})
facts['tetra_survey']={'searched':1024,'scout_positive':int(((td['windows'].min(1)>.003)&(td['speed']>1e-4)).sum()),'long_tested':len(tl['ids']),'long_positive':int(((tl['windows'].min(1)>.003)&(tl['speed']>1e-4)).sum()),'long_failed_ids':tl['ids'][~((tl['windows'].min(1)>.003)&(tl['speed']>1e-4))].tolist()}
facts['anderson_survey']={'contraction_search':4096,'extended_search':4096,'scout_positive':len(ad['positive_ids']),'long_positive':int(((al['windows'].min(1)>.005)&(al['status']==1)).sum()),'long_fixed':int((al['status']==0).sum()),'long_nonpositive_nonfixed':int(((al['status']==1)&~(al['windows'].min(1)>.005)).sum()),'long_failed_ids':al['ids'][~((al['windows'].min(1)>.005)&(al['status']==1))].tolist()}
(root/'facts.json').write_text(json.dumps(facts,indent=2));print(json.dumps(facts,indent=2))
