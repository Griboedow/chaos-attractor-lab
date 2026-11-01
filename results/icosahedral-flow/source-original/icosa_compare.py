"""Compare polyhedral candidates to a reproducible older map and a finite catalog."""
from pathlib import Path
import sys,json
import numpy as np
repo=Path('D:/Codex/chaos-attractor-lab');sys.path.insert(0,str(repo/'research'))
from compare_attractor_shapes import fingerprint,distance,align,normalize
from compare_affine_shapes import whiten
root=Path(__file__).parent;out=root/'icosa-comparison';out.mkdir(exist_ok=False)
refs={p.stem:np.load(p)['points'] for p in sorted((repo/'results/attractor-novelty/reference/trajectories').glob('*.npz'))}
refs['Reiter1999-Figure1']=np.load(root/'reiter-replay/trajectory.npz')['points']
for path,prefix in [('twist-geometry/long-runs.npz','OurTwist'),('gaussian-vortex/long-runs.npz','OurVortex')]:
    d=np.load(repo/'results'/path)
    for i,z in zip(d['ids'],d['tail']):refs[f'{prefix}-{i}']=z
candidates={f'Icosa-{i}':np.load(root/f'icosa-long/{i}-dense.npy')[::16] for i in [183,43,177]}
candidates.update({f'Tetra-{i}':np.load(repo/f'results/unexpected-search/tetra-render/sample-{i}.npz')['points'] for i in [421,565]})
rng=np.random.default_rng(2026110149);rows={};aligned={}
for mode,transform in [('similarity',normalize),('affine',whiten)]:
    rr={n:transform(z) for n,z in refs.items()};fps={n:fingerprint(z) for n,z in rr.items()}
    for name,raw in candidates.items():
        z=transform(raw);fp=fingerprint(z);ranks=sorted((float(distance(fp,v)),n) for n,v in fps.items())
        first,second=transform(raw[:len(raw)//2]),transform(raw[len(raw)//2:]);selfdist=float(distance(fingerprint(first),fingerprint(second,1)))
        a=first[rng.choice(len(first),1500,False)];b=second[rng.choice(len(second),1500,False)];selfrms,_=align(a,b)
        names=list(dict.fromkeys([n for _,n in ranks[:3]]+['Reiter1999-Figure1']));nearest=[]
        for n in names:
            a=z[rng.choice(len(z),1500,False)];b=rr[n][rng.choice(len(rr[n]),1500,False)];rms,aa=align(a,b)
            nearest.append({'reference':n,'fingerprint_distance':float(distance(fp,fps[n])),'best_local_icp_rms':rms})
            aligned[f'{mode}-{name}-{n}-candidate']=aa;aligned[f'{mode}-{name}-{n}-reference']=normalize(b)
        rows[f'{mode}-{name}']={'time_halves_fingerprint':selfdist,'time_halves_icp':selfrms,'comparisons':nearest,'top10':[{'reference':n,'distance':v} for v,n in ranks[:10]]}
        print(json.dumps({'mode':mode,'candidate':name,'self':selfrms,'nearest':nearest[0],'Reiter':next(a for a in nearest if a['reference']=='Reiter1999-Figure1')}),flush=True)
    np.savez_compressed(out/'aligned.npz',**aligned)
    (out/'report.json').write_text(json.dumps({'reference_count':len(refs),'literature_regimes':109,'comparisons':rows,'limits':'Finite empirical geometry, not historical novelty or dynamical equivalence. Local ICP with 48 starts permits reflections and is not globally certified. Whitening removes affine covariance anisotropy but not nonlinear transforms or invariant-measure differences. Reiter Figure 2 and the full 1997 tetrahedral map remain uncovered.'},indent=2))
