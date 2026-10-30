import hashlib,json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'research'))
from flint import arb
from validate_motor_saddle_period import count_windows
target=root/'results/motor-saddle-proof/zero-window-controls.json'
if target.exists():raise FileExistsError(target)
cases=[('owned initial zero',0,1,arb(1),True,False,True,1,0),
       ('excluded final zero',-1,0,arb(1),False,True,True,0,0),
       ('interior upward crossing',-1,1,arb(1),False,False,True,1,0),
       ('interior downward crossing',1,-1,arb(-1),False,False,True,0,1),
       ('overestimated tube with positive endpoint signs',1,1,arb(1),False,False,True,0,0),
       ('unresolved derivative sign',-1,1,arb(0,2),False,False,False,0,0),
       ('unowned initial endpoint',0,1,arb(1),False,False,False,0,0)]
checks=[]
for name,left,right,derivative,initial,final,passed,up,down in cases:
    actual=count_windows(0,[[arb(0,2)]],[[derivative]],arb(left),arb(right),initial,final)
    checks.append(dict(name=name,passed=(actual['passed'],actual['upward'],actual['downward'])==(passed,up,down),result=actual))
record=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,
    period_source_sha256=hashlib.sha256((root/'research/validate_motor_saddle_period.py').read_bytes()).hexdigest(),
    limitation='Boundary and rejection controls only; the actual period theorem also requires validated tubes, exact endpoint identities and the certified orbit.')
target.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
target.with_suffix('.source.py').write_bytes(Path(__file__).read_bytes())
print(record['status'],len(checks))
if record['status']!='passed':raise SystemExit(2)
