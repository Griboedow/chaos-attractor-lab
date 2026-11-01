"""Independent symbolic Jacobian identity and interval circle-link certificate."""
import argparse, hashlib, json
from pathlib import Path
import sympy as sp
import mpmath as mp


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    beta,a,b,c,s,lam=sp.symbols('beta a b c s lam', real=True)
    u=sp.Matrix([1,0,0]);v=sp.Matrix([c,s,0])
    J=beta*sp.eye(3)-beta*a*u*u.T+(1-2*beta)*b*v*v.T+2*beta*a*b*c*u*v.T
    reduce=lambda p:sp.rem(sp.expand(p),s*s+c*c-1,s).expand()
    T=2*beta-beta*a+(1-2*beta)*b+2*beta*a*b*c*c
    D=beta*((1-a)*(beta+(1-2*beta)*b)+a*b*c*c)
    assert reduce(J.det()-beta*D)==0
    assert reduce((lam*sp.eye(3)-J).det()-(lam-beta)*(lam*lam-T*lam+D))==0
    iv=mp.iv;iv.dps=80
    def exact(x):
        n,d=float(x).as_integer_ratio();return iv.mpf(n)/d
    rows=[]
    for case in json.loads(args.input.read_text())['cases']:
        be,r,theta,cx,cy,cz=map(exact,case['params']);st=iv.sin(theta);ct=iv.cos(theta)
        sn=-cz/(r*st);crossings=[];link=0
        if sn.a>1 or sn.b < -1:
            reason='B does not cross the plane of A'
        else:
            assert sn.a > -1 and sn.b < 1
            csabs=iv.sqrt(1-sn*sn)
            for sign in [-1,1]:
                cs=sign*csabs;xx=cx+r*cs;yy=cy+r*ct*sn
                radial=xx*xx+yy*yy;vertical=r*st*cs
                assert radial.b<1 or radial.a>1
                assert vertical.a>0 or vertical.b<0
                inside=bool(radial.b<1)
                crossing_sign=1 if vertical.a>0 else -1
                if inside:link+=crossing_sign
                crossings.append(dict(cosine_sign=sign,radius_squared=str(radial),vertical_derivative=str(vertical),inside_A_disk=inside,oriented_sign=crossing_sign))
            reason='Signed transverse intersections with the oriented unit disk'
        rows.append(dict(id=case['id'],linking_number=link,orientation='A is counterclockwise in xy; B follows increasing gamma; A disk normal +z',reason=reason,crossings=crossings))
    source=Path(__file__);(args.out/source.name).write_bytes(source.read_bytes());(args.out/'input.json').write_bytes(args.input.read_bytes())
    result=dict(status='passed',sympy=sp.__version__,mpmath=mp.__version__,interval_digits=80,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),identities=['det J = beta D','characteristic polynomial = (lambda-beta)(lambda^2-T lambda+D)'],circle_links=rows,scope='Algebra and geometric linking of constraints only; no classification of attractor topology or conjugacy.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
