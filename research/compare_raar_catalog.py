"""Finite geometric comparison against the frozen 108-regime reference set."""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np
from compare_attractor_shapes import fingerprint,distance,align,normalize
from compare_affine_shapes import whiten

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(args):
    repo=Path(__file__).resolve().parents[1];out=args.out;out.mkdir(parents=True,exist_ok=False)
    source=repo/'results/raar-circles/render/case-1051-animation.npz'
    with np.load(source) as f:z=f['states']
    paths=sorted((repo/'results/attractor-novelty/reference/trajectories').glob('*.npz'))
    refs={p.stem:np.load(p)['points'] for p in paths};assert len(refs)==108
    result={'candidate':'RAAR R1051','source_sha256':sha(source),'reference_files':{p.name:sha(p) for p in paths},'representations':{},'scope':'Finite sampled geometry only. Neither metric separation nor nearest matches establish dynamical equivalence, inequivalence, chaos or historical uniqueness. The catalog does not exhaust known parameters or nonlinear coordinate changes.'}
    arrays={};rng=np.random.default_rng(2026110121)
    for label,transform in [('similarity',normalize),('affine',whiten)]:
        transformed=transform(z);rr={n:transform(v) for n,v in refs.items()};fps={n:fingerprint(v) for n,v in rr.items()}
        ranked=sorted((float(distance(fingerprint(transformed),fp)),n) for n,fp in fps.items())
        a=transform(z[:500000]);b=transform(z[500000:]);self_fp=float(distance(fingerprint(a),fingerprint(b,1)))
        self_icp,_=align(a[rng.choice(len(a),1500,False)],b[rng.choice(len(b),1500,False)])
        rows=[]
        for dist,name in ranked[:3]:
            ai=rng.choice(len(z),1500,False);bi=rng.choice(len(refs[name]),1500,False)
            rms,aligned=align(transformed[ai],rr[name][bi])
            rows.append(dict(reference=name,fingerprint_distance=dist,best_local_icp_rms=rms))
            arrays[f'{label}--{name}--candidate']=aligned;arrays[f'{label}--{name}--reference']=normalize(rr[name][bi]);arrays[f'{label}--{name}--candidate_indices']=ai;arrays[f'{label}--{name}--reference_indices']=bi
        stability=[]
        for seed in [101,102,103,104]:
            random=np.random.default_rng(seed)
            cf=fingerprint(transform(z[random.choice(len(z),6000,False)]),seed)
            ranks=sorted((float(distance(cf,fingerprint(transform(v[random.choice(len(v),6000,False)]),seed))),name) for name,v in refs.items())
            stability.append(dict(seed=seed,nearest=ranks[0][1],distance=ranks[0][0]))
        result['representations'][label]=dict(nearest=rows,within_time_halves_fingerprint=self_fp,within_time_halves_icp=self_icp,subsample_rankings=stability,all_rankings=[dict(reference=n,distance=d) for d,n in ranked])
        print(label,json.dumps({k:v for k,v in result['representations'][label].items() if k!='all_rankings'}),flush=True)
    affine=np.array([[2.,.7,-.2],[0.,.4,.1],[.3,0.,1.1]])
    control=float(distance(fingerprint(whiten(z)),fingerprint(whiten(z@affine+[3,-7,11]))));assert control<1e-12
    result['exact_affine_control']=control
    result['icp_scope']='Best local alignment of 1500-point samples among 48 PCA starts; reflections allowed. Not a certified global minimum.'
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');np.savez_compressed(out/'aligned-samples.npz',**arrays)
    for name in ['compare_raar_catalog.py','compare_attractor_shapes.py','compare_affine_shapes.py']:shutil.copy2(Path(__file__).with_name(name),out/name)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args())
