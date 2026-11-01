"""Second map realization for the full critical-crossing face calculation.

Scalar circle coordinates, nested reflections and automatic differentiation
replace the vector projector/Jacobian model. The enclosure, subdivision driver
and mpmath interval library remain shared, explicitly recorded below.
"""
import argparse
import json
from pathlib import Path
import mpmath as mp
import certify_raar_critical_crossing as driver
from certify_raar_critical_doubleton import image
from check_raar_coverings_ad import Model as ScalarModel, Dual


class Model:
    def __init__(self,parameters):
        self.scalar=ScalarModel(parameters)

    def step(self,x,jacobian=True):
        if not jacobian:
            y,axis=self.scalar.field(x)
            return y,None,axis
        variables=[Dual(x[i],[mp.iv.mpf(int(i==j)) for j in range(3)]) for i in range(3)]
        values,axis=self.scalar.field(variables)
        return [v.v for v in values],[v.d for v in values],axis


def main(args):
    driver.image=image;driver.Model=Model
    sources=['check_raar_critical_crossing_ad.py','check_raar_coverings_ad.py','certify_raar_critical_doubleton.py']
    try:
        driver.main(args)
    finally:
        if args.out.exists():
            for name in sources:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    path=args.out/'report.json';report=json.loads(path.read_text())
    report.update(second_implementation='Scalar circle coordinates, nested reflectors and forward automatic derivatives. Shared doubleton enclosure, partition driver and mpmath interval library.',second_sources={name:driver.sha(args.out/name) for name in sources})
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'second_implementation':report['second_implementation'],'seconds':report['seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input',type=Path,default=Path('results/raar-symbolic-dynamics/runs/raar-covering-refined01/input.json'))
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--window',type=int,default=0)
    p.add_argument('--steps',type=int,default=42)
    p.add_argument('--orientation',type=int,choices=[-1,1],default=1)
    p.add_argument('--face',type=int,choices=[-1,0,1],default=0)
    p.add_argument('--digits',type=int,default=100)
    p.add_argument('--max-cells',type=int,default=20000)
    main(p.parse_args())
