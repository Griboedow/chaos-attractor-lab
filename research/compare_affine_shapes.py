"""Covariance-whitened follow-up: remove nonsingular affine coordinate changes.

Finite empirical distributions, not an attractor-equivalence or novelty test.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from compare_attractor_shapes import fingerprint, distance, align

ROOT=Path(__file__).resolve().parents[1]


def whiten(z):
    z=z-z.mean(0)
    ev,basis=np.linalg.eigh(np.cov(z.T))
    if ev[0] <= 1e-12*ev[-1]:
        raise ValueError('Empirical covariance is effectively singular')
    return z@basis@np.diag(1/np.sqrt(ev))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=Path,default=ROOT/'results/attractor-novelty')
    base=parser.parse_args().base
    refs={};excluded={}
    for path in sorted((base/'reference/trajectories').glob('*.npz')):
        try:
            refs[path.stem]=whiten(np.load(path)['points'])
        except ValueError as e:
            excluded[path.stem]=str(e)
    old=np.load(ROOT/'results/twist-geometry/long-runs.npz')
    new=np.load(ROOT/'results/gaussian-vortex/long-runs.npz')
    candidates={f'Twist-{i}':z for i,z in zip(old['ids'],old['tail'])}
    candidates.update({f'Vortex-{i}':z for i,z in zip(new['ids'][1::2],new['tail'][1::2])})
    rng=np.random.default_rng(20261111)
    descriptors={n:fingerprint(z) for n,z in refs.items()}
    rows={};clouds={}
    for name,z in candidates.items():
        white=whiten(z);fp=fingerprint(white)
        ranks=sorted((float(distance(fp,d)),n) for n,d in descriptors.items())
        a=whiten(z[:len(z)//2]);b=whiten(z[len(z)//2:])
        within=float(distance(fingerprint(a),fingerprint(b,1)))
        self_icp,_=align(a[rng.choice(len(a),1500,replace=False)],b[rng.choice(len(b),1500,replace=False)])
        nearest=refs[ranks[0][1]]
        rms,aligned=align(white[rng.choice(len(white),1500,replace=False)],nearest[rng.choice(len(nearest),1500,replace=False)])
        clouds[name]=aligned
        rows[name]={'within_time_halves':within,'within_time_halves_icp':self_icp,
                    'nearest':[{'reference':n,'fingerprint_distance':d} for d,n in ranks[:5]],
                    'nearest_best_icp_symmetric_rms':rms}
        print(name,json.dumps(rows[name]),flush=True)
    raw=np.load(base/'reference/trajectories/Lorenz.npz')['points']
    affine=np.array([[2.,.7,-.2],[0.,.4,.1],[.3,0.,1.1]])
    error=float(distance(fingerprint(whiten(raw)),fingerprint(whiten(raw@affine+[3,-7,11]))))
    assert error<1e-12
    info={'reference_count':len(refs),'excluded':excluded,'candidates':rows,
          'exact_affine_control':error,
          'meaning':'Whitening removes empirical covariance anisotropy. Affinely related full-rank clouds whiten to orthogonally related clouds. Sampling, time reparameterization, different invariant measures and nonlinear coordinate changes remain outside this test.',
          'icp_limit':'Best local alignment among 48 PCA starts, not a certified global minimum.'}
    (base/'comparison/affine-summary.json').write_text(json.dumps(info,indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(base/'comparison/affine-aligned-candidates.npz',**clouds)


if __name__=='__main__':
    main()
