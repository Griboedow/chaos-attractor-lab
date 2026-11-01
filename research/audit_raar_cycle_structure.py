"""Interval consequences of the isolated base seven-cycle certificate."""
import argparse,hashlib,json
from pathlib import Path
import mpmath as mp


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);mp.iv.dps=200;iv=mp.iv;r=json.loads(args.source.read_text())
    def parse(s):return iv.mpf([v.strip() for v in s.strip('[]').split(',')])
    a1,a2,a3=map(parse,r['characteristic_coefficients'])
    def q(x):return ((x+a1)*x+a2)*x+a3
    root_brackets=[['-2.571','-2.569'],['0.008','0.009'],['0.784','0.785']];rows=[]
    for left,right in root_brackets:
        values=[q(iv.mpf(v)) for v in [left,right]];signs=[]
        for v in values:
            assert v.a>0 or v.b<0;signs.append(1 if v.a>0 else -1)
        assert signs[0]*signs[1]<0;rows.append(dict(bracket=[left,right],endpoint_values=[str(v) for v in values],signs=signs))
    determinant=-a3;assert determinant.b<0 and q(iv.mpf(1)).a>0
    source=Path(__file__);(args.out/source.name).write_bytes(source.read_bytes());(args.out/'input.json').write_bytes(args.source.read_bytes())
    report=dict(status='passed',period=7,interval_digits=200,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),input_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),multiplier_brackets=rows,determinant=str(determinant),characteristic_at_plus_one=str(q(iv.mpf(1))),unstable_real_multipliers=1,stable_real_multipliers=2,scope='The cubic has three distinct real roots, one below -1 and two in (0,1). Its determinant is negative. Thus the seventh-iterate germ is locally invertible and orientation reversing. No C1 conjugacy near this cycle to a fixed-time map of a C1 ODE on R3 is possible: ODE solution-map Jacobians have positive determinant, preserved by periodic-orbit conjugacy. This is an application of standard invariants, not a new attractor-class or historical-priority proof.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
