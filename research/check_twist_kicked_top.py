"""Check the explicit affine embedding of the 1988 pumped kicked-top map.

The source Eq. (3) is transcribed from the author's public full-text preview.
No source PDF download, interval proof or exhaustive equivalence search is used.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from twist_map import step, step_tangent


S = np.array([[0.,0.,1.],[0.,-1.,0.],[1.,0.,0.]])
E = np.array([1.,0.,0.])


def original(z,b,a):
    x,y,w=z; c,s=np.cos(x),np.sin(x)
    return np.array([b*(y*c-w*s),-b*x,a*(1-b)+b*(y*s+w*c)])


def source_to_twist(z,b,a):
    return (S@z/a-(1-b)*E)/b


def twist_to_source(z,b,a):
    return a*S.T@(b*z+(1-b)*E)


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(2026110120)
    pairs=[(.1,35.),(.5,-5.),(.5,15.),(.7,32.),(.9,128.)]
    states=rng.uniform(-2,2,(1024,3));records=[]
    for b,a in pairs:
        p=np.array([b,0.,-np.pi/2,0.,0.,0.,a,0.,0.,0.])
        errors=[];jacerrors=[];deterrors=[]
        for z in states:
            source=twist_to_source(z,b,a)
            expected=source_to_twist(original(source,b,a),b,a)
            actual,j=step_tangent(z,np.eye(3),p)
            js=np.column_stack([np.imag(original(source.astype(complex)+1e-25j*np.eye(3)[i],b,a))/1e-25 for i in range(3)])
            errors.append(np.max(abs(actual-expected))/(1+np.max(abs(expected))))
            jacerrors.append(np.max(abs(j-S@js@S.T))/(1+np.max(abs(j))))
            deterrors.append(abs(np.linalg.det(js)-b**3))
        assert max(errors)<1e-12 and max(jacerrors)<1e-11 and max(deterrors)<1e-10
        # Source stability thresholds at the trivial fixed point (0,0,A).
        # Eigenvalue b and roots of lambda^2+b*A*lambda+b^2.
        fixed=np.array([0.,0.,a]);res=np.linalg.norm(original(fixed,b,a)-fixed)
        js=np.column_stack([np.imag(original(fixed.astype(complex)+1e-25j*np.eye(3)[i],b,a))/1e-25 for i in range(3)])
        poly=np.poly(js);expectedpoly=np.polymul([1.,-b],[1.,b*a,b*b])
        assert np.max(abs(poly-expectedpoly))<1e-10 and res<1e-14
        records.append(dict(b=b,A=a,parameters=p.tolist(),maximum_scaled_state_error=max(errors),
            maximum_scaled_jacobian_error=max(jacerrors),maximum_determinant_error=max(deterrors),
            trivial_fixed_point_residual=float(res),characteristic_polynomial_error=float(np.max(abs(poly-expectedpoly)))))
    np.savez_compressed(out/'controls.npz',states=states,pairs=pairs,rotation=S)
    payload=dict(source='Peplowski and Stefanski (1988), Physics Letters A132,408-413,Eq.(3),b_perp=b_parallel=b.',
        doi='10.1016/0375-9601(88)90503-8',url='https://www.academia.edu/113918730/Kicked_top_with_pumping_and_dissipation_Classical_results',
        access='Public author-linked full-text OCR read; typeset PDF download was restricted and not bypassed.',
        original_map=['b*(Y*cos(X)-Z*sin(X))','-b*X','A*(1-b)+b*(Y*sin(X)+Z*cos(X))'],
        coordinate_change='C(X,Y,Z)=((Z,-Y,X)/A-(1-b)*(1,0,0))/b; A!=0,0<b<1.',
        identity='C composed with G equals F composed with C, with p=(b,0,-pi/2,0,0,0,A,0,0,0).',
        seed=2026110120,states_per_pair=len(states),checks=records,
        scope='Explicit algebraic family inclusion plus floating-point controls. The selected1709/3254tuples are not this special slice; no conjugacy classification of their attractors is established.')
    (out/'analysis.json').write_text(json.dumps(payload,indent=2)+'\n')
    for name in ['check_twist_kicked_top.py','twist_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(dict(checks=len(pairs)*len(states),maximum_state_error=max(r['maximum_scaled_state_error'] for r in records),maximum_jacobian_error=max(r['maximum_scaled_jacobian_error'] for r in records))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
