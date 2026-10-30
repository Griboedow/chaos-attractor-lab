import hashlib,json,math
from pathlib import Path
import numpy as np
from scipy.linalg import schur
root=Path(__file__).resolve().parents[1]
source=root/'results/motor-saddle-proof/attempt-02.json'
output=root/'results/motor-saddle-proof/spectral-basis.json'
if output.exists():raise FileExistsError(output)
d=json.loads(source.read_text())
m=np.array([[math.ldexp(float(int(v['midpoint'][0])),v['midpoint'][1]) for v in row] for row in d['spectrum']['quotient_full_monodromy']])
t,q=schur(m,output='real')
record=dict(source=source.relative_to(root).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    basis=[[[str(a),str(b)] for a,b in (float(v).as_integer_ratio() for v in row)] for row in q],
    method='SciPy real Schur vectors of the enclosed midpoint, used only as an exact dyadic candidate similarity basis.',
    orthogonality_error=float(np.max(abs(q.T@q-np.eye(17)))),
    similarity_error=float(np.max(abs(q@t@q.T-m))),
    limitation='No approximate eigenvalue or orthogonality assertion is trusted; interval inversion must certify the proposed exact basis.')
output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
output.with_suffix('.source.py').write_bytes(Path(__file__).read_bytes())
print(record['orthogonality_error'],record['similarity_error'])
