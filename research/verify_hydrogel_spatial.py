"""Check the retained spatial hydrogel packet, optionally replaying algebra."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PACKET=ROOT/'results/hydrogel-spatial-controls'


def numerical():
    import numpy as np
    from scipy.optimize import brentq
    import hydrogel_spatial as m
    rng=np.random.default_rng(260409487);n=8
    beta,zeta,xi,chi,eta=m.parameters()
    q=np.array([m.expq(v) for v in rng.normal(scale=0.1,size=(n,3))])
    k=rng.normal(scale=0.5,size=(n,2))
    grad=m.energy_gradient(q,k,beta)[1];fd=np.empty_like(grad)
    for j in range(3*n):
        inc=np.zeros((n,3));inc[j//3,j%3]=1e-6
        plus=m.energy_gradient(m.rotate_increment(q,inc),k,beta)[0]
        minus=m.energy_gradient(m.rotate_increment(q,-inc),k,beta)[0]
        fd[j//3,j%3]=(plus-minus)/2e-6
    energy_error=float(np.max(abs(fd-grad))/np.max(abs(fd)))
    assert energy_error<1e-8
    no=m.geometry(m.rotate_increment(q,rng.normal(scale=0.003,size=(n,3))))[1]
    npv=m.geometry(m.rotate_increment(q,rng.normal(scale=0.003,size=(n,3))))[1]
    kp=k+rng.normal(scale=0.01,size=k.shape);errors=[]
    for order in (1,2):
        jac=m.mechanical_jacobian(q,k,kp,no,npv,beta,zeta,xi,chi,eta,0.002,order)
        fd=np.empty_like(jac)
        for j in range(3*n):
            inc=np.zeros((n,3));inc[j//3,j%3]=1e-6
            plus=m.residual_q(m.rotate_increment(q,inc),k,no,kp,npv,beta,zeta,xi,chi,eta,0.002,order)[0]
            minus=m.residual_q(m.rotate_increment(q,-inc),k,no,kp,npv,beta,zeta,xi,chi,eta,0.002,order)[0]
            fd[:,j]=((plus-minus)/2e-6).ravel()
        error=float(np.max(abs(fd-jac))/np.max(abs(fd)))
        assert error<1e-8;errors.append(error)
    n=16;alpha=33;theta=np.zeros(n)
    theta[1]=brentq(lambda x:x-alpha/n*np.sin(x/2),1e-7,np.pi,xtol=1e-14)
    for i in range(2,n):
        prev=theta[i-1]
        theta[i]=brentq(lambda x:x-prev-alpha/n*np.sin((prev+x)/2),prev-alpha/n,prev+alpha/n,xtol=1e-14)
    q=np.array([m.expq(np.array([-v,0.0,0.0])) for v in theta])
    k=np.zeros((n,2));k[1:,0]=-np.diff(theta)*n
    beta,zeta,xi,chi,eta=m.parameters(eta=1,chi=-alpha)
    nodes=m.geometry(q)[1]
    residual=float(np.max(abs(m.residual_q(q,k,nodes,k,nodes,beta,zeta,xi,chi,eta,0.004,2)[0])))
    assert residual<1e-10
    return dict(energy_gradient_relative_error=energy_error,jacobian_relative_errors=errors,discrete_equilibrium_scaled_residual=residual)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--numerical',action='store_true');args=ap.parse_args()
    manifest=json.loads((PACKET/'manifest.json').read_text(encoding='utf-8'))
    for item in manifest['files']:
        p=PACKET/item['path']
        assert p.stat().st_size==item['bytes'],p
        assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],p
    for item in manifest['dependencies']:
        p=ROOT/item['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],p
    runs={p.stem:json.loads(p.read_text(encoding='utf-8')) for p in (PACKET/'runs').glob('*.json')}
    assert len(runs)==12 and sum(r['success'] for r in runs.values())==10
    assert set(name for name,r in runs.items() if not r['success'])=={'spatial-16-dt002','scout-150-n16'}
    arrays=json.loads((PACKET/'arrays.json').read_text(encoding='utf-8'))['files']
    assert set(a['run'] for a in arrays)==set(runs)
    for a in arrays:
        assert a['arrays']['q']['shape'][1:]==[runs[a['run']]['segments'],4]
        assert a['arrays']['k']['shape'][1:]==[runs[a['run']]['segments'],2]
    clearance={r['run']:r for r in json.loads((PACKET/'clearance.json').read_text(encoding='utf-8'))}
    assert clearance['scout-150-n64']['samples_below_threshold']==67
    assert clearance['spatial-32-dt002']['samples_below_threshold']==0
    eq=json.loads((PACKET/'spurious-equilibrium.json').read_text(encoding='utf-8'))
    assert eq['scaled_mechanical_residual']<1e-10 and eq['relaxation_residual']<1e-10
    result=dict(retained_hashes=len(manifest['files'])+len(manifest['dependencies']),runs=len(runs),completed=10,failed=2,
                numerical=numerical() if args.numerical else 'not requested',
                scope='Retained packet and optional algebra only; no external-array replay or continuum/chaos certificate.')
    print(json.dumps(result))


if __name__=='__main__':main()
