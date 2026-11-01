"""Independent tree, union-find and multiplier-product orientation audit."""
import argparse,hashlib,json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bounds(s):
    a,b=map(Q,s.strip('[]').split(','));assert a<=b;return a,b
def multiply(a,b):
    products=[x*y for x in a for y in b];return min(products),max(products)

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);record=json.loads((args.source/'report.json').read_text(encoding='utf-8'));assert record['status']=='passed';assert sha(args.source/'certify_soft_local_injectivity.py')==record['source_sha256'];assert sha(args.source/'spanning-tree.npz')==record['tree_sha256'];assert sha(args.trap)==record['trap_sha256']
    with np.load(args.trap) as f:keys=list(map(int,f['occupied']));shape=list(map(int,f['shape']));assert str(Q(float(f['step'])))==record['step_exact'];assert [str(Q(float(v))) for v in f['origin']]==record['origin_exact']
    with np.load(args.source/'spanning-tree.npz') as f:order=f['visit_order'];parent=f['parent']
    n=len(keys);assert n==record['cubes'] and keys==sorted(set(keys));assert np.array_equal(np.sort(order),np.arange(n));assert len(parent)==n;seen=set();depth={};coordinates={}
    for pos,key in enumerate(keys):
        x,rem=divmod(key,shape[1]*shape[2]);y,z=divmod(rem,shape[2]);assert 0<=x<shape[0] and 0<=y<shape[1] and 0<=z<shape[2];coordinates[pos]=(x,y,z)
    for counter,pos in enumerate(map(int,order)):
        before=int(parent[pos])
        if counter==0:assert before==-1 and keys[pos]==record['tree_root_key'];depth[pos]=0
        else:
            assert before in seen;assert sum(abs(a-b) for a,b in zip(coordinates[pos],coordinates[before]))==1;depth[pos]=depth[before]+1
        seen.add(pos)
    assert len(seen)==n and max(depth.values())==record['tree_maximum_depth'] and record['tree_edges']==n-1
    # Separate connectivity algorithm: union only positive-coordinate neighbors.
    roots=list(range(n));size=[1]*n;lookup={v:k for k,v in coordinates.items()};edge_count=0
    def find(i):
        while roots[i]!=i:roots[i]=roots[roots[i]];i=roots[i]
        return i
    for i,coord in coordinates.items():
        for axis in range(3):
            adjacent=list(coord);adjacent[axis]+=1;j=lookup.get(tuple(adjacent))
            if j is None:continue
            edge_count+=1;a,b=find(i),find(j)
            if a!=b:
                if size[a]<size[b]:a,b=b,a
                roots[b]=a;size[a]+=size[b]
    assert len({find(i) for i in range(n)})==1 and edge_count==record['undirected_face_edges'];orientations=[]
    for label,path in [('covariance',args.cycles),('automatic_derivative',args.control)]:
        assert sha(path)==record['cycle_sources'][label+'_report_sha256'];report=json.loads(path.read_text(encoding='utf-8'));assert report['trap_sha256']==record['trap_sha256']
        for ident,sign in [(0,1),(5,-1)]:
            cycle=next(c for c in report['cases'] if c['id']==ident);assert cycle['period']==(1 if ident==0 else 5);product=(Q(1),Q(1));previous=None
            for bracket in cycle['multiplier_brackets']:
                root=bounds(bracket['interval']);left=bounds(bracket['at_left']);right=bounds(bracket['at_right']);assert left[1]<0<right[0] or right[1]<0<left[0];assert previous is None or previous<root[0];previous=root[1];product=multiply(product,root)
            assert product[0]>0 if sign==1 else product[1]<0
            original=next(r for r in record['orientation_witnesses'] if r['engine']==label and r['cycle_id']==ident);d=list(map(Q,original['determinant_interval_exact']));assert max(d[0],product[0])<=min(d[1],product[1]);assert original['sign']==sign;orientations.append(dict(engine=label,cycle_id=ident,sign=sign,multiplier_product_interval_exact=list(map(str,product))))
    result=dict(status='passed',primary_report_sha256=sha(args.source/'report.json'),tree_complete=True,tree_edges_checked=n-1,union_find_components=1,positive_neighbor_edges=edge_count,orientation_controls=orientations,source_sha256=sha(Path(__file__)),scope='Every exact tree edge and full cube coverage checked; independent union-find confirms connectedness. Independent multiplier-product bounds confirm both determinant signs and overlap the direct determinant enclosures. The attracting-set and root existence results are frozen dependencies; the topological argument is not machine formalized.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps({k:v for k,v in result.items() if k!='orientation_controls'}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--trap',type=Path,required=True);p.add_argument('--cycles',type=Path,required=True);p.add_argument('--control',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
