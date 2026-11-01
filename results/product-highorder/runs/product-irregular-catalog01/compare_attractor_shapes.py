"""Finite sample geometric comparison; no topological or priority inference."""
from pathlib import Path
import argparse,json,itertools
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]

def normalize(z):
 z=np.asarray(z);z=z-z.mean(0);return z/np.sqrt(np.mean((z*z).sum(1)))

def fingerprint(z,seed=0):
 z=normalize(z);rng=np.random.default_rng(seed)
 ii=rng.integers(len(z),size=(60000,2));d=np.linalg.norm(z[ii[:,0]]-z[ii[:,1]],axis=1)
 h=np.histogram(d,bins=np.r_[np.linspace(0,4,65),np.inf])[0].astype(float);h=np.sqrt(h/h.sum())
 r=np.linalg.norm(z,axis=1);g=np.histogram(r,bins=np.r_[np.linspace(0,3,33),np.inf])[0].astype(float);g=np.sqrt(g/g.sum())
 ev=np.linalg.eigvalsh(np.cov(z.T))[::-1];ev/=ev.sum()
 return np.r_[h,ev,g]

def distance(a,b):
 return np.linalg.norm(a[:65]-b[:65])+.5*np.linalg.norm(a[65:68]-b[65:68])+.25*np.linalg.norm(a[68:]-b[68:])

def align(a,b):
 """Best of 48 PCA initializations, allowing mirror images; local ICP."""
 a=normalize(a);b=normalize(b)
 va,aa=np.linalg.eigh(np.cov(a.T));vb,bb=np.linalg.eigh(np.cov(b.T))
 tree=cKDTree(b);best=(np.inf,None)
 for perm in itertools.permutations(range(3)):
  for signs in itertools.product([-1,1],repeat=3):
   rotation=aa@np.eye(3)[:,perm]@np.diag(signs)@bb.T
   z=a@rotation
   for iteration in range(20):
    _,idx=tree.query(z,workers=1);target=b[idx]
    center_a=z.mean(0);center_b=target.mean(0)
    u,_,vt=np.linalg.svd((z-center_a).T@(target-center_b))
    rot=u@vt
    z=(z-center_a)@rot+center_b
   forward=tree.query(z,workers=1)[0];back=cKDTree(z).query(b,workers=1)[0]
   residual=np.sqrt((np.mean(forward**2)+np.mean(back**2))/2)
   if residual<best[0]:best=(float(residual),z)
 return best

def main():
 parser=argparse.ArgumentParser()
 parser.add_argument('--base',type=Path,default=ROOT/'results/attractor-novelty')
 BASE=parser.parse_args().base
 out=BASE/'comparison';out.mkdir(exist_ok=True)
 refs={path.stem:np.load(path)['points'] for path in sorted((BASE/'reference/trajectories').glob('*.npz'))}
 twists=np.load(ROOT/'results/twist-geometry/long-runs.npz')
 new=np.load(ROOT/'results/gaussian-vortex/long-runs.npz')
 candidates={f'Twist-{i}':z for i,z in zip(twists['ids'],twists['tail'])}
 candidates.update({f'Vortex-{i}':z for i,z in zip(new['ids'][1::2],new['tail'][1::2])})
 descriptors={name:fingerprint(z) for name,z in refs.items()}
 results={};alignments={};rng=np.random.default_rng(20261109)
 for name,z in candidates.items():
  fp=fingerprint(z);ranks=sorted([(distance(fp,d),n) for n,d in descriptors.items()])
  split=distance(fingerprint(z[::2]),fingerprint(z[1::2],1))
  # Disjoint temporal halves are stronger than interleaved samples for a flow.
  half=distance(fingerprint(z[:len(z)//2]),fingerprint(z[len(z)//2:],1))
  rows=[]
  for dist,refname in ranks[:5]:
   ref=refs[refname]
   a=z[rng.choice(len(z),1500,replace=False)];b=ref[rng.choice(len(ref),1500,replace=False)]
   residual,aligned=align(a,b)
   rows.append({'reference':refname,'fingerprint_distance':float(dist),'best_icp_symmetric_rms':residual})
   alignments[f'{name}--{refname}--candidate']=aligned
   alignments[f'{name}--{refname}--reference']=normalize(b)
  # Same-system point-cloud residual: a finite-sample scale for the ICP output.
  a=z[:len(z)//2];b=z[len(z)//2:]
  self_residual,_=align(a[rng.choice(len(a),1500,replace=False)],b[rng.choice(len(b),1500,replace=False)])
  results[name]={'within_interleaved':float(split),'within_time_halves':float(half),
                 'within_time_halves_icp':self_residual,'nearest':rows}
  print(name,json.dumps(results[name]),flush=True)
 # Exact rotation/translation/scale positive control; descriptor ignores these.
 a=refs['Lorenz'];theta=.71
 rot=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
 control=distance(fingerprint(a),fingerprint(3.7*(a@rot)+[11,-4,7]))
 assert control<1e-12
 sample=a[:1500]
 icp_control,_=align(sample,3.7*(sample@rot)+[11,-4,7])
 assert icp_control<1e-10
 # Repeat the descriptor ranking with independently drawn 6,000-point samples.
 stability={n:[] for n in candidates}
 for seed in [101,102,103,104]:
  random=np.random.default_rng(seed)
  ref_fps={n:fingerprint(z[random.choice(len(z),6000,replace=False)],seed) for n,z in refs.items()}
  for name,z in candidates.items():
   fp=fingerprint(z[random.choice(len(z),6000,replace=False)],seed)
   dist,nearest=min((distance(fp,d),n) for n,d in ref_fps.items())
   stability[name].append({'seed':seed,'nearest':nearest,'distance':float(dist)})
 within={n:float(distance(fingerprint(z[:len(z)//2]),fingerprint(z[len(z)//2:],1))) for n,z in refs.items()}
 info={'reference_count':len(refs),'reference_names':sorted(refs),'candidates':results,
       'exact_similarity_control':float(control),'exact_similarity_icp_control':icp_control,
       'subsample_rankings':stability,'reference_time_halves':within,
       'distance_definition':'Hellinger-embedding L2(D2,65 bins) + 0.5*L2(covariance eigenvalue fractions) + 0.25*L2(radial Hellinger embedding,33 bins)',
       'scope':'89 eligible catalog defaults, four external atlas examples, four published 3D Henon regimes, nine successful Halvorsen parameter variants, and two archived pumped-kicked-top regimes. Two additional Halvorsen runs failed; one prior kicked-top outcome is a fixed point. These were excluded. Not an exhaustive search over parameters or coordinate changes.',
       'icp_limit':'Local optimization with 48 PCA starts, uniform RMS normalization and reflections allowed. Residual is best alignment found, not a certified global minimum.'}
 (out/'summary.json').write_text(json.dumps(info,indent=2)+'\n')
 np.savez_compressed(out/'aligned-clouds.npz',**alignments)
 np.savez_compressed(out/'fingerprints.npz',**descriptors)

if __name__=='__main__':main()
