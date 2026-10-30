"""Replay instability brackets from the exact serialized full return matrices."""
import argparse,datetime,hashlib,json
from pathlib import Path
from flint import arb,arb_mat,ctx
from validate_motor_cycle import read_ball,ball_record

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose an unused output')
    ctx.prec=768
    source=json.loads(args.input.read_text())
    results=[]
    for name,block in source['singleton_blocks'].items():
        matrix=arb_mat([[read_ball(v) for v in row] for row in block['full_return']])
        polynomial=matrix.charpoly()
        left,right=polynomial(arb(-32)),polynomial(arb(-1))
        passed=bool(left<0 and right>0)
        results.append(dict(cluster=name,passed=passed,bracket=[-32,-1],
            polynomial_at_left=ball_record(left),polynomial_at_right=ball_record(right)))
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status='passed' if all(v['passed'] for v in results) else 'failed',bits=768,
        input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=results,
        scope='Recompute characteristic polynomials from exact saved interval matrices at higher precision. This checks serialization and strict real-root brackets; it does not reintegrate the variational equations.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    args.output.with_suffix('.source.py').write_bytes(Path(__file__).read_bytes())
    print(result['status'],[(v['cluster'],v['bracket']) for v in results])
    assert result['status']=='passed'

if __name__=='__main__':main()
