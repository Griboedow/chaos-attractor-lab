"""Numerical search of regular preimage chains near a repelling cycle."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from snapback_seed import inverse
from cycle_seeds import orbit_eval
from reproduce import raw_parameters


def search(args):
    rec=json.loads((args.source/f'case-{args.case}.json').read_text());raw=raw_parameters(rec)
    certs=json.loads((args.source/'cycles'/f'cycles-{args.case}.json').read_text())['cycles']
    rep=next(v for v in certs if v['proof']['unstable_dimension']==3)
    p=rep['period'];x=np.array(rep['center'],float)
    _,m,orb=orbit_eval(x,raw,p)
    vals,u=np.linalg.eig(m);complex_idx=np.flatnonzero(vals.imag>0)
    if not len(complex_idx):B=u.real
    else:
        j=complex_idx[0];i=np.flatnonzero(abs(vals.imag)<1e-10)[0]
        B=np.column_stack([u[:,i].real,u[:,j].real,u[:,j].imag])
    bi=np.linalg.inv(B);mi=np.linalg.inv(m)
    with np.load(rec['cloud']) as f:z=f['points'].astype(float)
    tree=cKDTree(z);rng=np.random.default_rng(args.case+421)
    population=[(x,[],False)];best=1e10;args.out.mkdir(parents=True,exist_ok=True)
    for depth in range(1,args.depth+1):
        phase=orb[(-depth)%p];out=[]
        try: linear=np.linalg.inv(orbit_eval(phase,raw,1)[1])
        except: linear=np.eye(3)
        for target,chain,excursion in population:
            _,near=tree.query(target,k=90)
            seeds=[phase,phase+linear@(target-orb[(1-depth)%p])]
            for j in near:
                if j>0 and all(np.linalg.norm(z[j-1]-s)>.045 for s in seeds):seeds.append(z[j-1])
                if len(seeds)>=16:break
            seeds.extend(z[rng.choice(len(z),args.random,replace=False)])
            for seed in seeds:
                a=inverse(target,seed,raw,1)
                if a is None or np.linalg.norm(a)>10 or any(np.linalg.norm(a-y[0])<1e-6 for y in out):continue
                if any(np.linalg.norm(a-np.asarray(q['seed']))<1e-8 for q in chain[-p:]) and excursion: continue
                away=excursion or min(np.linalg.norm(a-q) for q in orb)>.012
                nex=chain+[dict(period=1,seed=a.tolist())]
                dist=np.linalg.norm(a-phase)
                out.append((a,nex,away))
                if not away or depth%p:continue
                if dist<best:
                    best=dist;print('near',depth,best,flush=True)
                if dist>.12:continue
                b=a;chain2=list(nex);old=np.max(abs(bi@(b-x)))
                for k in range(60):
                    options=[x+mi@(b-x),x,.5*(x+a)]
                    preds=[s for opt in options if (s:=inverse(b,opt,raw,p)) is not None]
                    if not preds:break
                    b=min(preds,key=lambda s:np.max(abs(bi@(s-x))))
                    diff=np.max(abs(bi@(b-x)))
                    if diff>max(.18,2*old):break
                    old=diff;chain2.append(dict(period=p,seed=b.tolist()))
                    if diff<3e-7:
                        dest=args.out/f'tree-{args.case}.json'
                        dest.write_text(json.dumps(dict(case=args.case,phase=0,period=p,base=x.tolist(),basis=B.tolist(),distance=diff,blocks=chain2),indent=2))
                        print('FOUND',args.case,depth,len(chain2),diff,flush=True);return
        if not out: break
        out.sort(key=lambda t:np.linalg.norm(t[0]-phase))
        # Preserve some distant branches across skipped dangerous regions.
        population=out[:args.beam//2]
        extras=out[args.beam//2:]
        if extras:population += [extras[j] for j in rng.choice(len(extras),min(len(extras),args.beam-len(population)),replace=False)]
        print('depth',depth,len(out),'best',best,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--case',type=int,required=True);p.add_argument('--depth',type=int,default=16)
    p.add_argument('--beam',type=int,default=48);p.add_argument('--random',type=int,default=8)
    search(p.parse_args())
