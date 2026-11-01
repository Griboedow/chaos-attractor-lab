from paths import OUTPUT, ROOT, RESEARCH
"""A finite geometric audit, including earlier candidates from this repository."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
repo=ROOT;sys.path.insert(0,str(repo/'research'))
from compare_attractor_shapes import fingerprint,distance,align,normalize
from compare_affine_shapes import whiten

root=OUTPUT;out=root/'comparison';out.mkdir(exist_ok=False)
refs={p.stem:np.load(p)['points'] for p in sorted((repo/'results/attractor-novelty/reference/trajectories').glob('*.npz'))}
literature_count=len(refs)
for path,prefix in [('twist-geometry/long-runs.npz','OurTwist'),('gaussian-vortex/long-runs.npz','OurVortex')]:
    d=np.load(repo/'results'/path)
    for i,z in zip(d['ids'],d['tail']):refs[f'{prefix}-{i}']=z
td=np.load(root/'tetra-spectra/data.npz');ad=np.load(root/'anderson-spectra/data.npz')
candidates={f'Tetra-{i}':td['tail'][np.where((td['ids']==i)&(td['h']==.0125))[0][0]] for i in [421,565]}
candidates.update({f'RNA-{i}':ad['tail'][np.where((ad['ids']==i)&(ad['modes']==0))[0][0]] for i in [54,449,2268]})
rng=np.random.default_rng(2026110138);rows={};aligned={};excluded={}
for mode,transform in [('similarity',normalize),('affine',whiten)]:
    rr={}
    for n,z in refs.items():
        try:rr[n]=transform(z)
        except ValueError as e:excluded[n]=str(e)
    fps={n:fingerprint(z) for n,z in rr.items()}
    for name,raw in candidates.items():
        z=transform(raw);fp=fingerprint(z);ranks=sorted((float(distance(fp,v)),n) for n,v in fps.items())
        first,second=transform(raw[:len(raw)//2]),transform(raw[len(raw)//2:]);selfdist=float(distance(fingerprint(first),fingerprint(second,1)))
        a=first[rng.choice(len(first),1200,False)];b=second[rng.choice(len(second),1200,False)];selfrms,_=align(a,b)
        nearest=[]
        for dist,n in ranks[:3]:
            a=z[rng.choice(len(z),1200,False)];b=rr[n][rng.choice(len(rr[n]),1200,False)];rms,aa=align(a,b)
            nearest.append({'reference':n,'fingerprint_distance':dist,'best_local_icp_rms':rms})
            aligned[f'{mode}-{name}-{n}-candidate']=aa;aligned[f'{mode}-{name}-{n}-reference']=normalize(b)
        row={'time_halves_fingerprint':selfdist,'time_halves_icp':selfrms,'nearest':nearest,'top10':[{'reference':n,'distance':v} for v,n in ranks[:10]]};rows[f'{mode}-{name}']=row
        print(json.dumps({'mode':mode,'candidate':name,'self':selfrms,'nearest':nearest[0]}),flush=True)
    np.savez_compressed(out/'aligned.npz',**aligned)
    (out/'report.json').write_text(json.dumps({'known_regime_count':literature_count,'including_earlier_repo_candidates':len(refs),'excluded':excluded,'comparisons':rows,'limitations':'No calibrated novelty probability. A finite cloud comparison does not prove orbit equivalence or inequivalence. ICP uses 48 local starts; it is not a certified global optimization. Different time parameterizations and invariant measures can change these statistics. Reiter tetrahedral maps and Swift-Barany C3 normal forms are literature gaps, not covered by these 108 regimes.'},indent=2))
# Exact affine control with every geometric statistic computed by the same functions.
z=refs['Lorenz'];a=np.array([[2,.3,.7],[0,.4,-.2],[.5,0,1.3]])
error=float(distance(fingerprint(whiten(z)),fingerprint(whiten(z@a+[3,-7,11]))))
assert error<1e-12
(out/'control.json').write_text(json.dumps({'exact_affine_fingerprint_error':error},indent=2))
