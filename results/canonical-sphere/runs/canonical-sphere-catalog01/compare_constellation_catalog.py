"""Finite geometry comparison; neither a novelty test nor a conjugacy test."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from compare_attractor_shapes import fingerprint,distance,align,normalize
from compare_affine_shapes import whiten
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from data_cache import data_path,records


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    root=Path(__file__).resolve().parents[1];args.out.mkdir(parents=True,exist_ok=False)
    names=sorted(n for n in records() if n.startswith('results/attractor-novelty/reference/trajectories/'))
    local={p.relative_to(root).as_posix() for p in (root/'results/attractor-novelty/reference/trajectories').glob('*.npz')}
    names=sorted(set(names)|local);assert len(names)==108,len(names)
    refs={Path(n).stem:np.load(data_path(n))['points'] for n in names}
    result=dict(ids=args.ids,reference_files={n:sha(data_path(n)) for n in names},cases={},scope='108 frozen reference regimes, not all known models or parameters. Sampled geometric distance cannot establish dynamical inequivalence, topology or historical uniqueness. Whitening removes affine distortion; it does not test nonlinear conjugacy.',alignment='1500-point samples; 48 PCA starts, reflected fits allowed, 20 local ICP steps. Best found residual is not a global minimum.')
    for name in [Path(__file__).name,'compare_attractor_shapes.py','compare_affine_shapes.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    arrays={};rng=np.random.default_rng(2026110118);begin=time.time()
    for label,transform in [('similarity',normalize),('affine',whiten)]:
        rr={};excluded={}
        for n,z in refs.items():
            try:rr[n]=transform(z)
            except ValueError as e:excluded[n]=str(e)
        fps={n:fingerprint(z) for n,z in rr.items()}
        for ident in args.ids:
            path=args.source/f'case-{ident}.npz'
            with np.load(path) as f:z=f['states']
            transformed=transform(z);fp=fingerprint(transformed)
            ranking=sorted((float(distance(fp,v)),n) for n,v in fps.items())
            a=transform(z[:len(z)//2]);b=transform(z[len(z)//2:]);self_fp=float(distance(fingerprint(a),fingerprint(b,1)))
            self_rms,_=align(a[rng.choice(len(a),1500,False)],b[rng.choice(len(b),1500,False)])
            nearest=[]
            for d,n in ranking[:3]:
                ai=rng.choice(len(z),1500,False);bi=rng.choice(len(rr[n]),1500,False)
                rms,aligned=align(transformed[ai],rr[n][bi]);prefix=f'{ident}--{label}--{n}'
                arrays[prefix+'--candidate']=aligned;arrays[prefix+'--reference']=normalize(rr[n][bi]);arrays[prefix+'--candidate_indices']=ai;arrays[prefix+'--reference_indices']=bi
                nearest.append(dict(reference=n,fingerprint_distance=d,best_local_icp_rms=rms))
            stability=[]
            for seed in [101,102,103,104]:
                random=np.random.default_rng(seed);cf=fingerprint(transform(z[random.choice(len(z),6000,False)]),seed)
                d,n=min((float(distance(cf,fingerprint(transform(v[random.choice(len(v),6000,False)]),seed))),n) for n,v in refs.items() if n not in excluded)
                stability.append(dict(seed=seed,reference=n,distance=d))
            record=dict(source_sha256=sha(path),within_time_halves_fingerprint=self_fp,within_time_halves_icp=self_rms,nearest=nearest,subsample_rankings=stability,all_rankings=[dict(reference=n,distance=d) for d,n in ranking],excluded=excluded)
            result['cases'].setdefault(str(ident),{})[label]=record
            (args.out/'partial.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
            print(json.dumps(dict(id=ident,representation=label,self_rms=self_rms,nearest=nearest,seconds=time.time()-begin)),flush=True)
    z=refs['Lorenz'];affine=np.array([[2.,.7,-.2],[0.,.4,.1],[.3,0.,1.1]])
    control=float(distance(fingerprint(whiten(z)),fingerprint(whiten(z@affine+[3,-7,11]))));assert control<1e-12
    result['exact_affine_control']=control
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');np.savez_compressed(args.out/'aligned-samples.npz',**arrays)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())
