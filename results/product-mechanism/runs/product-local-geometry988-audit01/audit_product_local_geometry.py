"""Exact endpoint audit and frozen-source replay of local geometry certificates."""
import argparse, hashlib, json, subprocess, sys
from fractions import Fraction as Q
from pathlib import Path
from audit_product_orientation import interval,point,add,sub,mul,isum,contains,mm,mv,det,nonzero,opposite


def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(v):return max(abs(t) for t in v)


def audit(folder,args):
    r=read(folder/'report.json');assert r['status']=='passed'
    for name in ['input','singular','proposals','critical']:assert r[name+'_sha256']==sha(getattr(args,name))
    for name,h in r['sources'].items():assert sha(folder/name)==h
    target=list(map(interval,read(args.singular)['krawczyk']))
    assert len(r['preimages'])==2 and len(r['folds'])==6
    for p in r['preimages']:
        x=list(map(point,p['center']));rad=Q(p['radius']);B=list(map(interval,p['box']));K=list(map(interval,p['krawczyk']))
        C=[[point(v) for v in row] for row in p['preconditioner']];J=[[interval(v) for v in row] for row in p['derivative']]
        R=[[interval(v) for v in row] for row in p['remainder']];CH=mm(C,J);res=list(map(interval,p['residual']));correction=mv(C,res)
        for i in range(3):
            assert contains(interval(p['target'][i]),target[i]);assert contains(B[i],add(x[i],(-rad,rad)))
            for j in range(3):assert contains(R[i][j],sub(point(i==j),CH[i][j]))
            bound=add(sub(x[i],correction[i]),isum(mul(R[i][j],(-rad,rad)) for j in range(3)))
            assert contains(K[i],bound) and B[i][0]<K[i][0]<=K[i][1]<B[i][1]
        assert max(sum(norm(v) for v in row) for row in R)<1 and nonzero(det(C))
        assert interval(p['normalization'])[0]>0 and contains(interval(p['determinant']),det(J)) and nonzero(det(J))
    assert {(p['axis'],p['side']) for p in r['folds']}=={(a,s) for a in range(3) for s in [-1,1]}
    for p in r['folds']:
        B=list(map(interval,p['box']));rad=Q(p['radius']);axis=p['axis']
        for j in range(3):assert contains(B[j],add(point(p['center'][j]),(-rad,rad) if j==axis else point(0)))
        assert opposite(*map(interval,p['endpoint_determinants'])) and interval(p['normalization'])[0]>0
        J=[[interval(v) for v in row] for row in p['derivative']]
        H=[[[interval(v) for v in row] for row in mat] for mat in p['hessian']]
        grad=[]
        for k in range(3):
            grad.append(isum(det([[H[i][j][k] if row==i else J[row][j] for j in range(3)] for row in range(3)]) for i in range(3)))
        G=list(map(interval,p['determinant_gradient']));assert all(contains(g,v) for g,v in zip(G,grad)) and nonzero(G[axis])
        kernel=[sub(mul(J[0][a],J[1][b]),mul(J[0][b],J[1][a])) for a,b in [(1,2),(2,0),(0,1)]]
        n=list(map(interval,p['kernel']));assert all(contains(a,b) for a,b in zip(n,kernel)) and any(nonzero(v) for v in n)
        folded=isum(mul(a,b) for a,b in zip(G,n));assert contains(interval(p['fold_transversality']),folded) and nonzero(folded)
    return dict(run=folder.name,report_sha256=sha(folder/'report.json'),preimages=2,folds=6)


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);args.replays.mkdir(parents=True,exist_ok=False)
    for name in [Path(__file__).name,'audit_product_orientation.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[audit(p,args) for p in args.runs];replays=[]
    for folder in args.runs:
        r=read(folder/'report.json');dest=args.replays/folder.name
        command=[sys.executable,str(folder/'certify_product_local_geometry.py'),'--out',str(dest),'--engine',r['engine'],'--digits',str(r['digits']),'--radius',r['preimages'][0]['radius']]
        for name in ['input','singular','proposals','critical','mpmath-path']:command+=['--'+name,str(getattr(args,name.replace('-','_')).resolve())]
        subprocess.run(command,check=True,capture_output=True,text=True)
        new=read(dest/'report.json')
        for key in ['preimages','folds','independent_expression_errors']:assert new[key]==r[key],key
        replays.append(dict(run=folder.name,exact_payloads=True))
    report=dict(status='passed',certificates=rows,replays=replays,
                scope='Exact rational Krawczyk, rank, determinant-gradient and transversality arithmetic; frozen-source replays of both expressions. Both transcendental enclosures use mpmath.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','singular','proposals','critical','out','replays','mpmath-path']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--runs',type=Path,nargs='+',required=True);main(p.parse_args())
