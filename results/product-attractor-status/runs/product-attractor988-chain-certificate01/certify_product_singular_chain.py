"""Certify a finite preimage chain ending at the normalization singularity.

A root in a specified neighborhood prevents that whole neighborhood from
belonging to any forward-invariant smooth-domain set. It does not exclude a
smaller trapping region, or put the root on the observed invariant set.
"""
import argparse,hashlib,json,sys,time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from product_binary_interval import parameters,model,mm,mv
from certify_product_cycles import det


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    sys.path.append(str(args.mpmath_path));import mpmath as mp
    import verify_product_coverings as full
    mp.mp.dps=mp.iv.dps=args.digits;iv=mp.iv;full.mp=mp
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    exact=lambda x,ctx=iv:ctx.mpf(Q(x).numerator)/Q(x).denominator
    enc=lambda v:dict(decimal=str(v),binary=[list(t) for t in v._mpi_])
    dec=lambda v:mp.nstr(v,args.digits+10)
    save=lambda name,v:(args.out/name).write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n')
    def decode(v):
        e=[]
        for sign,m,exponent,bits in v['binary']:
            assert int(m).bit_length()==bits;e.append(exact(Q((-1 if sign else 1)*m)*Q(2)**exponent))
        return iv.mpf([e[0].a,e[1].b])
    sources={}
    for name in [Path(__file__).name,'product_binary_interval.py','certify_product_cycles.py','verify_product_coverings.py']:
        source=Path(__file__).with_name(name);(args.out/name).write_bytes(source.read_bytes());sources[name]=sha(source)
    raw=json.loads(args.input.read_text());singular=json.loads(args.singular.read_text());nodes=json.loads(args.nodes.read_text())
    assert singular['status']=='passed' and singular['input_sha256']==sha(args.input)
    selected=nodes[args.node] if args.node is not None else min(nodes[1:],key=lambda v:v['cloud_distance'])
    sequence=[];row=selected
    while row['parent'] is not None:sequence.append(row);row=nodes[row['parent']]
    sequence.reverse();assert [r['depth'] for r in sequence]==list(range(1,len(sequence)+1))
    save('input-chain.json',sequence)
    par=parameters(raw['params'],np.array(raw['points']),mp.mp);pari=parameters(raw['params'],np.array(raw['points']),iv)
    independent=full.Model(raw);field=independent.value_jacobian if args.engine=='full' else lambda x:model(x,pari,iv)
    numeric=lambda x:model(x,par,mp.mp)
    S=list(map(decode,singular['krawczyk']));target=[exact(v,mp.mp) for v in singular['center']]
    radius=exact(args.radius);chain=[]
    for proposal in sequence:
        root=mp.findroot(lambda *x:[a-b for a,b in zip(numeric(x)[0],target)],proposal['state'],tol=mp.mpf(10)**(-args.digits+25),maxsteps=40)
        center=list(map(dec,root));xi=list(map(exact,center));box=[v+iv.mpf([-radius.b,radius.b]) for v in xi]
        y0,_,_=field(xi);_,J,r=field(box);assert r.a>0
        Ctext=[[dec(v) for v in row] for row in (mp.matrix(numeric(root)[1])**-1).tolist()];C=[[exact(v) for v in row] for row in Ctext]
        residual=[a-b for a,b in zip(y0,S)];R=[[int(i==j)-v for j,v in enumerate(row)] for i,row in enumerate(mm(C,J))]
        correction=mv(C,residual);rest=mv(R,[iv.mpf([-radius.b,radius.b])]*3);K=[xi[i]-correction[i]+rest[i] for i in range(3)]
        norm=max(sum(abs(v).b for v in row) for row in R)
        assert norm<1 and all(k.a>b.a and k.b<b.b for k,b in zip(K,box)) and abs(det(C)).a>0
        chain.append(dict(depth=proposal['depth'],proposal_id=proposal['id'],center=center,radius=args.radius,box=list(map(enc,box)),target=list(map(enc,S)),
            residual=list(map(enc,residual)),derivative=[[enc(v) for v in row] for row in J],preconditioner=Ctext,remainder=[[enc(v) for v in row] for row in R],
            krawczyk=list(map(enc,K)),contraction=enc(norm),normalization=enc(r),determinant=enc(det(J))))
        S=K;target=list(root)
    with np.load(args.source) as f:reference=f['states']
    _,idx=cKDTree(reference).query(np.array(list(map(float,target))));w=list(map(lambda v:exact(float(v)),reference[idx]))
    distance=iv.sqrt(sum((S[i]-w[i])**2 for i in range(3)));assert distance.b<exact(args.distance_bound).a
    save('chain.json',chain)
    result=dict(status='passed',engine=args.engine,digits=args.digits,sources=sources,input_sha256=sha(args.input),singular_sha256=sha(args.singular),nodes_sha256=sha(args.nodes),reference_sha256=sha(args.source),
        selected_id=selected['id'],depth=len(chain),root_radius=args.radius,reference_index=int(idx),reference_state=reference[idx].tolist(),distance=enc(distance),distance_bound=args.distance_bound,
        chain_sha256=sha(args.out/'chain.json'),minimum_normalization=min(float(decode(v['normalization']).a) for v in chain),seconds=time.time()-begin,
        scope='A certified point z lies within the stated distance of one retained reference state and satisfies F^depth(z)=s, where F is undefined. Any neighborhood containing z cannot be contained in a forward-invariant subset of the smooth domain. This does not disprove all smaller traps or any invariant-set membership hypothesis.')
    save('report.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['distance','sources']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','singular','nodes','source','mpmath-path','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--node',type=int);p.add_argument('--engine',choices=['factorized','full'],default='factorized')
    p.add_argument('--digits',type=int,default=140);p.add_argument('--radius',default='1e-40');p.add_argument('--distance-bound',default='0.00006')
    main(p.parse_args())
