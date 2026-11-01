"""Resolve cancellation in the original five-point ice-cloud Jacobian control."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ice_cloud import coefficients,make_normalized,field_jac,physical_field,to_physical


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    original=json.loads((args.source/'progress.json').read_text());assert len(original)==3
    for row in original:assert hashlib.sha256((args.source/row['file']).read_bytes()).hexdigest()==row['sha256']
    rng=np.random.default_rng(2026110111);rows=[]
    for w in [.01,.1,1.]:
        c=coefficients(updraft=w);p,scale,omega=make_normalized(c)
        for index in range(16):
            z=np.r_[rng.uniform(-2,2,2),rng.uniform(.85,1.6)]
            f,j=field_jac(z,p)
            def independent(v):
                x=to_physical(v,scale);return physical_field(x,c)/(np.r_[x[:2],1.]*omega)
            exact=np.column_stack([independent(z+1e-25j*np.eye(3)[i]).imag/1e-25 for i in range(3)])
            eps=1e-5;eye=np.eye(3)*eps
            finite=np.column_stack([(-field_jac(z+2*h,p)[0]+8*field_jac(z+h,p)[0]-8*field_jac(z-h,p)[0]+field_jac(z-2*h,p)[0])/(12*eps) for h in eye])
            ae=float(np.max(abs(f-independent(z))/(1+abs(f))));ce=float(np.max(abs(j-exact)/(1+abs(j))))
            fe=float(np.max(abs(j-finite)/(1+abs(j))))
            assert ae<1e-10 and ce<1e-10
            rows.append(dict(updraft=w,index=index,state=z.tolist(),field_error=ae,complex_step_scaled_error=ce,
                five_point_scaled_error=fe,analytic=j.tolist(),independent_complex_step=exact.tolist(),five_point=finite.tolist()))
    result=dict(passed=True,source=str(args.source),source_progress_sha256=hashlib.sha256((args.source/'progress.json').read_bytes()).hexdigest(),
        maximum_field_error=max(r['field_error'] for r in rows),maximum_complex_step_error=max(r['complex_step_scaled_error'] for r in rows),
        maximum_five_point_error=max(r['five_point_scaled_error'] for r in rows),rows=rows,
        explanation='Original long source integrations completed and independent trajectories agreed, but the final componentwise five-point Jacobian assertion failed. At far-from-equilibrium random states, exponential nucleation creates terms up to1e23 while some partial derivatives are order1; subtraction destroys these derivatives. Independent dimensional complex-step derivatives resolve every tested entry. Original failed packet remains unchanged apart from an added failure record.',
        scope='All test states stay away from s=1; complex-step verification applies to the selected smooth side of the evaporation law. No derivative is asserted across its kink.')
    (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in ['check_ice_cloud_algebra.py','ice_cloud.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print({k:v for k,v in result.items() if k not in ['rows','explanation','scope']})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
