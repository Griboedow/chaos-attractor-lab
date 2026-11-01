"""Compare complete old and accelerated 128-segment histories."""
from pathlib import Path
import hashlib
import json
import numpy as np

root=Path(__file__).resolve().parent
names=['restart128-dt006','accel128-dt006']
rows={}
for filename in ('trajectory.npz','checkpoint.npz'):
    with np.load(root/names[0]/filename) as a, np.load(root/names[1]/filename) as b:
        assert a.files==b.files
        errors={key:float(np.max(abs(a[key]-b[key]))) for key in a.files}
        assert max(errors.values())<1e-8,errors
        rows[filename]=errors
out=dict(runs=names,maximum_absolute_differences=rows,
    input_sha256={f'{name}/{file}':hashlib.sha256((root/name/file).read_bytes()).hexdigest()
                  for name in names for file in ('trajectory.npz','checkpoint.npz')},
    scope='Identical mesh, time step and endpoint-initialized IVP. Complete sampled histories and final BDF checkpoints compared; not a continuum refinement test.')
(root/'assembly-history-check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out))
