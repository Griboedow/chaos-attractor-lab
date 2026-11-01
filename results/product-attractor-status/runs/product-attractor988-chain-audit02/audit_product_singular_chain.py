"""Exact-rational audit and independent-expression replays of a singular chain."""
import argparse,hashlib,json,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from audit_product_orientation import endpoint,interval,point,add,sub,mul,isum,contains,mm,mv,det


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);read=lambda p:json.loads(p.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name in [Path(__file__).name,'audit_product_orientation.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    singular=read(args.singular);assert singular['status']=='passed';S=list(map(interval,singular['krawczyk']));checks=[]
    for folder in args.runs:
        r=read(folder/'report.json');chain=read(folder/'chain.json');assert r['status']=='passed' and r['singular_sha256']==sha(args.singular)
        assert sha(folder/'chain.json')==r['chain_sha256'] and len(chain)==r['depth']==27
        target=S
        for n,row in enumerate(chain,1):
            assert row['depth']==n and list(map(interval,row['target']))==target
            center=list(map(point,row['center']));radius=Q(row['radius']);B=list(map(interval,row['box']));K=list(map(interval,row['krawczyk']))
            C=[[point(v) for v in rr] for rr in row['preconditioner']];J=[[interval(v) for v in rr] for rr in row['derivative']];R=[[interval(v) for v in rr] for rr in row['remainder']]
            product=mm(C,J);correction=mv(C,list(map(interval,row['residual'])))
            for i in range(3):
                assert contains(B[i],add(center[i],(-radius,radius)))
                for j in range(3):assert contains(R[i][j],sub(point(i==j),product[i][j]))
                bound=add(sub(center[i],correction[i]),isum(mul(R[i][j],(-radius,radius)) for j in range(3)))
                assert contains(K[i],bound) and B[i][0]<K[i][0]<=K[i][1]<B[i][1]
            assert max(sum(max(abs(a),abs(b)) for a,b in rr) for rr in R)<1
            assert det(C)[0]==det(C)[1]!=0 and interval(row['normalization'])[0]>0
            target=K
        w=list(map(lambda v:point(Q(v)),r['reference_state']));distance2=sum(max(abs(a),abs(b))**2 for a,b in [sub(x,y) for x,y in zip(target,w)])
        assert distance2<Q(r['distance_bound'])**2 and interval(r['distance'])[1]**2>=distance2
        dest=args.out/'replays'/folder.name
        cmd=[sys.executable,str(folder/'certify_product_singular_chain.py'),'--input',str(args.input),'--singular',str(args.singular),'--nodes',str(args.nodes),'--source',str(args.source),'--mpmath-path',str(args.mpmath_path),'--node',str(r['selected_id']),'--engine',r['engine'],'--digits',str(r['digits']),'--radius',r['root_radius'],'--distance-bound',r['distance_bound'],'--out',str(dest)]
        subprocess.run(cmd,check=True,capture_output=True,text=True);replay=read(dest/'report.json')
        assert sha(dest/'chain.json')==sha(folder/'chain.json') and replay['distance']==r['distance']
        checks.append(dict(run=folder.name,report_sha256=sha(folder/'report.json'),chain_length=27,exact_krawczyk_arithmetic=True,reference_distance_upper=float(interval(r['distance'])[1]),exact_replay=True))
    # The proposed wide grid contains the previously certified central pole.
    with np.load(args.grid) as f:edges=f['edges'];shape=f['shape'];keys=f['occupied']
    middle=np.array([float((a+b)/2) for a,b in S]);ijk=np.array([np.searchsorted(edges[i,:shape[i]+1],middle[i],side='right')-1 for i in range(3)])
    key=int((ijk[0]*shape[1]+ijk[1])*shape[2]+ijk[2]);assert key in keys
    cell=[(Q(float(edges[i,ijk[i]])),Q(float(edges[i,ijk[i]+1]))) for i in range(3)]
    assert all(a<x[0]<=x[1]<b for (a,b),x in zip(cell,S))
    result=dict(status='passed',runs=checks,wide_grid=dict(sha256=sha(args.grid),central_pole_cell_key=key,indices=ijk.tolist(),contains_certified_pole=True),
        scope='Exact rational Krawczyk and finite-distance checks, frozen-source replays, and exact inclusion of the central singularity in the rejected proposed grid. Does not exclude every smaller trapping region or classify the asymptotic numerical regime.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','singular','nodes','source','grid','mpmath-path','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--runs',type=Path,nargs=2,required=True);main(p.parse_args())
