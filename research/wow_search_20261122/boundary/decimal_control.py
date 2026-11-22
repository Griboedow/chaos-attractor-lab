"""Independent standard-library Decimal control using Newton's original z - p/p' update."""
from decimal import Decimal as Q,localcontext
import json
from pathlib import Path
from selected import D,K
from search import OUTPUT

def add(u,v):return u[0]+v[0],u[1]+v[1]
def sub(u,v):return u[0]-v[0],u[1]-v[1]
def mul(u,v):return u[0]*v[0]-u[1]*v[1],u[0]*v[1]+u[1]*v[0]
def div(u,v):
    q=v[0]**2+v[1]**2
    return (u[0]*v[0]+u[1]*v[1])/q,(u[1]*v[0]-u[0]*v[1])/q
def norm(u):return u[0]**2+u[1]**2
def decode(value):return tuple(Q(int(m))*(Q(2)**int(e)) for m,e in value)
def step(z,a):
    x=(Q(1),Q(0))
    for i in range(D-1):x=mul(x,z)
    derivative=sub((D*x[0],D*x[1]),a)
    return sub(z,div(add(sub(mul(x,z),mul(a,z)),a),derivative))
def main():
    q=json.loads((OUTPUT/'certificate.json').read_text(encoding='utf-8'))
    with localcontext() as ctx:
        ctx.prec=160;origin=(Q(0),Q(0));a=decode(q['parameter_certificate']['center_dyadic']);r=[decode(m['center_dyadic']) for m in q['root_certificates']]
        outcomes=[]
        for witness in q['witnesses']:
            z=decode(witness['center_dyadic'])
            for it in range(451):
                rs=[i for i,root in enumerate(r) if norm(sub(z,root))<Q('1e-14')]
                if rs:label=rs[0];break
                if norm(z)<Q('1e-14'):label=D+((-it)%K);break
                z=step(z,a)
            else:raise RuntimeError('Did not reach any target')
            assert label==witness['certified_outcome'] and it==witness['capture_iteration']
            outcomes.append(dict(label=label,capture_iteration=it))
    (OUTPUT/'independent-decimal.json').write_text(json.dumps(dict(precision_decimal_digits=160,method="Newton numerator and derivative evaluated separately with Decimal arithmetic; z - p/p'.",records=outcomes),indent=2)+'\n')
if __name__=='__main__':main()
