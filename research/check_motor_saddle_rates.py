"""Known-cubic controls for the interval saddle-rate postprocessor."""
import argparse
import hashlib
import json
from pathlib import Path
from flint import arb, ctx
import certify_motor_saddle_rates as rates


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output')
    ctx.prec=256
    checks=[]
    cases=[('stable_real','.8',['-.005','-.05','1'],True),
        ('unstable_real','1.125',['-.005','-.05','1'],True),
        ('stable_complex','.8',['.0041','-.1','1'],True),
        ('extra_large_root','.8',['-.01','-.15','1'],False),
        ('repeated_large_root','.8',['-.04','-.75','1'],False)]
    for label,a,q,want in cases:
        a=arb(a);q=list(map(arb,q))
        p=[-a*q[0],q[0]-a*q[1],q[1]-a,arb(1)]
        lo='1.12' if label=='unstable_real' else '.79'
        hi='1.13' if label=='unstable_real' else '.81'
        try:
            root,result=rates.isolate(p,lo,hi)
            passed=result['passed'];contains=bool(root.contains(a))
        except ArithmeticError:
            passed=False;contains=None
        checks.append(dict(name=label,expected_pass=want,actual_pass=passed,
            exact_root_contained=contains,passed=passed==want and (contains if want else True)))
    sources=[Path(__file__),Path(rates.__file__)]
    record=dict(status='passed' if all(x['passed'] for x in checks) else 'failed',checks=checks,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        scope='Known exact cubic controls. The original midpoint-only routine correctly declined coarse intervals when the root equalled the midpoint; the refinement probes strict quarter endpoints without changing acceptance gates.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for p in sources:
        snapshot=args.output.with_name(args.output.stem+'.'+p.name)
        if snapshot.exists():raise FileExistsError('Snapshot exists')
        snapshot.write_bytes(p.read_bytes())
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(record['status'],len(checks),'checks')
    if record['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
