"""Numerical seed and exact dyadic similarity basis; neither is proof."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import schur
import motor_postcritical_cycles as cycles


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new seed directory')
    args.output.mkdir(parents=True)
    root=Path(__file__).resolve().parents[1]
    source=root/'results/motor-activity-continuation/probe-destinations/case-16/periodic/cycle.npz'
    data=np.load(source);initial=data['quotient_initial_state'];period=float(data['time'][-1]);p=data['parameters']
    weights=np.array([2,2,1,2,1])/8;permutation=[3,0,4,1,2]
    flow=cycles.integrate(initial,period/6,p,weights,augmented=True,rtol=5e-13,max_step=.01,dense=True)
    pm=np.zeros((14,14))
    for row,col in enumerate(permutation):
        pm[row,col]=1;pm[5+row,5+col]=1
        if row<4:
            if col<4:pm[10+row,10+col]=1
            else:pm[10+row,10:]=-weights[:-1]/weights[-1]
    matrix=np.linalg.solve(pm,flow.y[14:,-1].reshape(14,14))
    triangular,basis=schur(matrix,output='real')
    seed=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_path=source.relative_to(root).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        exact_rho='129/55',exact_gain='1+(19/20)/sqrt(32001/500000)',segments=3,
        states=[flow.sol(i*period/18)[:14].tolist() for i in range(3)],segment_duration=period/18,
        relative_period_guess=period/6,full_period_guess=period,permutation=permutation,
        relative_eigenvalues=[dict(real=float(z.real),imag=float(z.imag),modulus=float(abs(z))) for z in np.linalg.eigvals(matrix)],
        basis=[[list(float(v).as_integer_ratio()) for v in row] for row in basis],
        qualification='Floating-point proposals only. The interval driver must validate every inclusion and spectrum.',
        source_code_sha256={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [Path(__file__),Path(cycles.__file__)]})
    (args.output/'seed.json').write_text(json.dumps(seed,indent=2)+'\n',encoding='utf-8')
    for x in [Path(__file__),Path(cycles.__file__)]: (args.output/x.name).write_bytes(x.read_bytes())
    print('segment',period/18,'relative eigenvalue moduli',sorted(abs(np.linalg.eigvals(matrix)),reverse=True))


if __name__=='__main__':main()
