import hashlib,json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
source=root/'results/motor-postcritical/six-cluster-cycle-rho-3p268/cycle.npz'
output=root/'results/motor-saddle-proof/multiple-shooting-seed.json'
if output.exists():raise FileExistsError(output)
with np.load(source) as data:
    indices=[0,125,250,375]
    record=dict(source=source.relative_to(root).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        indices=indices,states=data['cluster_states'][indices].tolist(),times=data['time'][indices].tolist(),
        period=float(data['time'][-1]),parameters=data['parameters'].tolist(),
        limitation='Ordinary floating-point centers only; no proof depends on their accuracy.')
output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
output.with_suffix('.source.py').write_bytes(Path(__file__).read_bytes())
print(record['period'],record['times'])
