"""Refine a numerical fold of periodic orbits using a bordered shooting system.

This is floating-point evidence, not an interval-validated bifurcation proof.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from continue_motor_braid_cycles import shoot, FREE, B_SCALE
from motor_rows import rhs


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=args.source/'cycles.npz';data=np.load(source)
    protocol=json.loads((args.source/'protocol.json').read_text())
    original=np.load(protocol['source']);rho,r=original['parameters'][:2]
    index=int(np.argmin(data['coordinates'][:,11]));y0=data['coordinates'][index]
    g,j,*_=shoot(y0,rho,r);v0=np.linalg.svd(j[:,:11])[2][-1]
    controls=dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        initial_branch_index=index,unknowns='10 initial coordinates, period,100*b;11D right nullvector',
        equations='G=0; G_(initial,period) v=0; dot(v,v)=1',
        scope='Numerical fold evidence only, not an interval certificate or proof of global dynamics.')
    (out/'protocol.json').write_text(json.dumps(controls,indent=2)+'\n')
    for name in ['refine_motor_cycle_fold.py','continue_motor_braid_cycles.py','check_motor_asymmetry.py','motor_rows.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    count=0;began=time.perf_counter()
    def equation(z):
        nonlocal count
        g,j,*_=shoot(z[:12],rho,r);v=z[12:];count+=1
        return np.r_[g,j[:,:11]@v,np.dot(v,v)-1]
    sol=root(equation,np.r_[y0,v0],tol=1e-9)
    residual=equation(sol.x);y=sol.x[:12];v=sol.x[12:]
    g,j,m,q,period,p=shoot(y,rho,r);u,s,vh=np.linalg.svd(j[:,:11]);left=u[:,-1]
    quadratic=[]
    for eps in [1e-3,3e-4,1e-4]:
        delta=np.r_[v,0.]*eps
        value=(shoot(y+delta,rho,r)[0]-2*g+shoot(y-delta,rho,r)[0])/eps**2
        quadratic.append(dict(epsilon=eps,left_projection=float(left@value)))
    independent=solve_ivp(lambda t,z:rhs(t,z,p,4),(0,period),q,method='Radau',rtol=2e-12,atol=2e-14,max_step=.005)
    assert independent.success
    report=dict(success=bool(sol.success),message=sol.message,evaluations=count,wall_seconds=time.perf_counter()-began,
        parameters=p.tolist(),period=float(period),maximum_residual=float(np.max(abs(residual))),
        shooting_singular_values=s.tolist(),b_transversality=float(left@j[:,11]*B_SCALE),
        quadratic_projections=quadratic,section_velocity=float(rhs(0,q,p,4)[8]),
        independent_Radau_closure=float(np.max(abs(independent.y[:,-1]-q))),
        multipliers=[[float(z.real),float(z.imag)] for z in np.linalg.eigvals(m)],
        limitation='Finite differences of floating-point shooting derivatives; no rigorous error enclosure or normal-form theorem.')
    np.savez_compressed(out/'fold.npz',coordinates=y,nullvector=v,left_nullvector=left,seed=q,period=period,parameters=p,monodromy=m,derivative=j)
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    assert report['maximum_residual']<1e-8 and report['independent_Radau_closure']<1e-8


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
