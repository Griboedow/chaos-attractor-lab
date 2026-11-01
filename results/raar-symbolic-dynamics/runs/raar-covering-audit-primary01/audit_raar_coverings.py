"""Audit complete dyadic partitions and every strict covering inequality."""
import argparse,gzip,hashlib,json,math
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bounds(s):
    a,b=[Decimal(v.strip()) for v in s.strip('[]').split(',')];assert a.is_finite() and b.is_finite() and a<=b;return a,b

def rational_inverse(A):
    a=[row.copy()+[Fraction(i==j) for j in range(3)] for i,row in enumerate(A)]
    for k in range(3):
        pivot=next(i for i in range(k,3) if a[i][k]);a[k],a[pivot]=a[pivot],a[k];v=a[k][k];a[k]=[x/v for x in a[k]]
        for i in range(3):
            if i!=k:v=a[i][k];a[i]=[x-v*y for x,y in zip(a[i],a[k])]
    return [row[3:] for row in a]

def audit_partition(path,kind,orientation):
    tree={};count=0;depths=Counter();modes=Counter();minimum_axis=Decimal('Infinity');minimum_margin=Decimal('Infinity')
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            record=json.loads(line);branch=record['path'];assert len(branch)%2==0;node=tree
            for k in range(0,len(branch),2):
                token=branch[k:k+2];assert token[0] in '012' and token[1] in '01';assert kind=='volume' or token[0]!='0';assert '#' not in node;node=node.setdefault(token,{})
            assert not node;node['#']=True;depths[len(branch)//2]+=1;count+=1
            image=[bounds(v) for v in record['image']];axis=bounds(record['minimum_axis'])[0];assert axis>0;minimum_axis=min(minimum_axis,axis)
            if kind=='volume':
                if image[0][0]>1:margin=image[0][0]-1;mode='right-exterior'
                elif image[0][1]<-1:margin=-1-image[0][1];mode='left-exterior'
                else:
                    margin=min(v for a,b in image[1:] for v in [a+1,1-b]);assert margin>0;mode='stable-interior'
            else:
                sign=(-orientation if kind=='left' else orientation);margin=(image[0][0] if sign==1 else -image[0][1])-1;assert margin>0;mode='exit-face'
            minimum_margin=min(minimum_margin,margin);modes[mode]+=1
    def visit(node):
        if '#' in node:assert set(node)=={'#'};return 1
        keys=sorted(node);assert len(keys)==2 and keys[0][0]==keys[1][0] and keys[0][1]=='0' and keys[1][1]=='1';return sum(visit(v) for v in node.values())
    assert visit(tree)==count and sum(Fraction(n,2**d) for d,n in depths.items())==1
    return dict(cells=count,maximum_depth=max(depths),minimum_axis=str(minimum_axis),minimum_strict_margin=str(minimum_margin),classification_counts=dict(modes))

def main(args):
    repo=Path(__file__).resolve().parents[1];records=[];raw=None;input_digest=None
    for edge in args.edges:
        folder=args.base/(('raar-covering-proof-pilot01' if edge==0 else f'raar-covering-proof-edge{edge}-01') if args.mode=='primary' else f'raar-covering-ad-edge{edge}-02')
        report=json.loads((folder/'report.json').read_text());assert report['status']=='passed' and report['edges']==[edge];assert sha(folder/'input.json')==report['input_sha256'];assert sha(folder/'certify_raar_coverings.py')==report['source_sha256']
        if args.mode=='ad':assert sha(folder/'check_raar_coverings_ad.py')==report['second_source_sha256'] and report['interval_digits']==100
        else:assert report['interval_digits']==70
        if raw is None:raw=json.loads((folder/'input.json').read_text());input_digest=report['input_sha256']
        assert report['input_sha256']==input_digest;assert {r['kind'] for r in report['checks']}=={'volume','left','right'}
        for check in report['checks']:
            assert check['edge']==edge and check['status']=='passed' and check['pending']==0;path=folder/f"edge-{edge}-{check['kind']}.jsonl.gz";assert sha(path)==check['certificate_sha256'];row=audit_partition(path,check['kind'],int(raw['orientation'][edge]));assert row['cells']==check['accepted'] and 2*row['cells']-1==check['examined'] and row['maximum_depth']==check['maximum_depth'];records.append(dict(edge=edge,kind=check['kind'],**row))
    expected=[[0,0],[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,0]];assert raw['edges']==expected and len(raw['centers'])==7 and all(v in [-1,1] for v in raw['orientation'])
    matrices=[[[Fraction(v)*Fraction(r[j]) for j,v in enumerate(row)] for row in B] for B,r in zip(raw['bases'],raw['radii'])];inverses=[rational_inverse(A) for A in matrices];C=inverses[0]
    displacement=[Fraction(b)-Fraction(a) for a,b in zip(raw['centers'][0],raw['centers'][1])];center=sum(C[0][k]*displacement[k] for k in range(3));radius=sum(abs(sum(C[0][k]*matrices[1][k][j] for k in range(3))) for j in range(3));assert center+radius<-1 or center-radius>1
    cycle=json.loads((repo/'results/raar-cycle-structure/cycle-proof/base-unstable-period-7.json').read_text());assert raw['params']==cycle['parameters']
    cyclebox=[[Fraction(v)-Fraction(raw['centers'][0][i]) for v in bounds(s)] for i,s in enumerate(cycle['box'])]
    for i in range(3):
        intervals=[sorted([C[i][k]*v for v in cyclebox[k]]) for k in range(3)];assert sum(v[0] for v in intervals)>-1 and sum(v[1] for v in intervals)<1
    logarithm_lower=2*sum(Fraction(1,(2*k+1)*3**(2*k+1)) for k in range(6));assert logarithm_lower/49>Fraction(1414,100000)
    result=dict(status='passed',mode=args.mode,checked_edges=args.edges,complete_graph=set(args.edges)==set(range(8)),input_sha256=input_digest,checks=records,total_accepted_cells=sum(r['cells'] for r in records),global_minimum_axis=str(min(Decimal(r['minimum_axis']) for r in records)),window_0_1_separation_exact=[str(center-radius),str(center+radius)],certified_cycle_box_inside_window_0=True,entropy_bound='log(2)/49 > 0.01414 per F iterate, conditional on the complete graph and the covering theorem',scope='The audit validates stored strict inequalities, complete closed dyadic tilings, exact window invertibility/separation and source identity. It does not independently recompute the interval image enclosures or certify attraction/novelty.')
    args.out.mkdir(parents=True,exist_ok=False);(args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps({k:v for k,v in result.items() if k not in ['checks','window_0_1_separation_exact']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--mode',choices=['primary','ad'],required=True);p.add_argument('--edges',type=int,nargs='+',default=list(range(8)));p.add_argument('--out',type=Path,required=True);main(p.parse_args())
