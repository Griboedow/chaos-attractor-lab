"""Exact structural audit of soft-map window coverings; no image recomputation."""
import argparse,gzip,hashlib,json,math
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import numpy as np

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def interval(v):
    if isinstance(v,str):v=v.strip('[]').split(',')
    a,b=map(Q,v);assert a<=b;return a,b

def partition(path,kind,orientation):
    tree={};depths=Counter();modes=Counter();minimum=None;axis_min=None;count=0
    with gzip.open(path,'rt',encoding='utf-8') as handle:
        for line in handle:
            row=json.loads(line);address=row['path'];assert len(address)%2==0;node=tree
            for k in range(0,len(address),2):
                token=address[k:k+2];assert token[0] in '012' and token[1] in '01';assert kind=='volume' or token[0]!='0';assert '#' not in node;node=node.setdefault(token,{})
            assert not node;node['#']=True;count+=1;depths[len(address)//2]+=1
            image=[interval(v) for v in row['image']];axis=interval(row['minimum_axis'])[0] if isinstance(row['minimum_axis'],str) else Q(row['minimum_axis']);assert axis>0
            axis_min=axis if axis_min is None else min(axis_min,axis)
            if kind=='volume':
                if image[0][0]>1:margin=image[0][0]-1;mode='right-exterior'
                elif image[0][1]<-1:margin=-1-image[0][1];mode='left-exterior'
                else:margin=min(v for a,b in image[1:] for v in [a+1,1-b]);mode='stable-interior'
            else:
                sign=-orientation if kind=='left' else orientation;margin=(image[0][0] if sign==1 else -image[0][1])-1;mode='exit-face'
            assert margin>0;minimum=margin if minimum is None else min(minimum,margin);modes[mode]+=1
    def visit(node):
        if '#' in node:assert set(node)=={'#'};return 1
        keys=sorted(node);assert len(keys)==2 and keys[0][0]==keys[1][0] and keys[0][1]=='0' and keys[1][1]=='1';return sum(visit(v) for v in node.values())
    assert visit(tree)==count and sum(Q(n,2**d) for d,n in depths.items())==1
    return dict(cells=count,maximum_depth=max(depths),minimum_axis_exact=str(axis_min),minimum_strict_margin_exact=str(minimum),minimum_strict_margin=float(minimum),modes=dict(modes))

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);report=json.loads((args.source/'report.json').read_text());raw=json.loads((args.source/'input.json').read_text());assert report['status']=='passed' and sha(args.source/'input.json')==report['input_sha256']
    for name,digest in report['sources'].items():assert sha(args.source/name)==digest
    assert raw['edges']==[[0,0],[0,1],[1,2],[2,3],[3,4],[4,5],[5,0]] and len(raw['centers'])==6 and len(raw['orientation'])==7
    matrices=[];inverses=[]
    for B,r,stored in zip(raw['bases'],raw['radii'],raw['exact_windows']):
        S=[[Q(v)*Q(r[j]) for j,v in enumerate(row)] for row in B];T=[[Q(v) for v in row] for row in stored['inverse']];assert S==[[Q(v) for v in row] for row in stored['matrix']]
        assert all(sum(S[i][k]*T[k][j] for k in range(3))==int(i==j) for i in range(3) for j in range(3));matrices.append(S);inverses.append(T)
    C=inverses[0];d=[Q(b)-Q(a) for a,b in zip(raw['centers'][0],raw['centers'][1])];center=sum(C[0][k]*d[k] for k in range(3));radius=sum(abs(sum(C[0][k]*matrices[1][k][j] for k in range(3))) for j in range(3));assert center+radius<-1 or center-radius>1
    hull=[[Q(raw['centers'][0][i])+sign*sum(abs(v) for v in matrices[0][i]) for sign in [-1,1]] for i in range(3)]
    assert all(Q(raw['window0_hull'][i][0])<=hull[i][0]<=hull[i][1]<=Q(raw['window0_hull'][i][1]) for i in range(3))
    assert sha(args.trap)==raw['trap_sha256']
    with np.load(args.trap) as f:
        assert f['params'].tolist()==raw['params'] and f['points'].tolist()==raw['points'];origin=f['origin'];h=Q(float(f['step']));shape=f['shape'];occupied=set(f['occupied'].tolist())
    keys=[]
    for i in range(int((hull[0][0]-Q(float(origin[0])))//h),int((hull[0][1]-Q(float(origin[0])))//h)+1):
        for j in range(int((hull[1][0]-Q(float(origin[1])))//h),int((hull[1][1]-Q(float(origin[1])))//h)+1):
            for k in range(int((hull[2][0]-Q(float(origin[2])))//h),int((hull[2][1]-Q(float(origin[2])))//h)+1):
                assert 0<=i<shape[0] and 0<=j<shape[1] and 0<=k<shape[2];key=int((i*shape[1]+j)*shape[2]+k);assert key in occupied;keys.append(key)
    assert keys==raw['window0_hull_cube_keys'];checks=[];seen=set()
    for check in report['checks']:
        e=check['edge'];kind=check['kind'];assert e in range(7) and kind in ['volume','left','right'] and (e,kind) not in seen;seen.add((e,kind));assert check['status']=='passed' and check['pending']==0
        path=args.source/f'edge-{e}-{kind}.jsonl.gz';assert sha(path)==check['certificate_sha256'];row=partition(path,kind,raw['orientation'][e]);assert row['cells']==check['accepted'];assert row['maximum_depth']==check['maximum_depth'];assert 2*row['cells']-1==check['examined'];checks.append(dict(edge=e,kind=kind,**row))
    assert len(seen)==21
    lower=2*sum(Q(1,(2*k+1)*3**(2*k+1)) for k in range(8));assert lower/6>Q(11552,100000)
    result=dict(status='passed',source_report_sha256=sha(args.source/'report.json'),input_sha256=report['input_sha256'],complete_graph=True,total_leaves=sum(r['cells'] for r in checks),checks=checks,exact_separation=[str(center-radius),str(center+radius)],window0_in_trap=True,window0_hull_cube_keys=keys,entropy_bound='log(2)/6 > 0.11552 per original iterate, using the covering theorem and the separately certified trap',scope='Exact complete-tiling, inequality, matrix, separation, trapping-membership and hash audit. Enclosures themselves require independent recomputation.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps({k:v for k,v in result.items() if k not in ['checks','exact_separation','window0_hull_cube_keys']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--trap',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
