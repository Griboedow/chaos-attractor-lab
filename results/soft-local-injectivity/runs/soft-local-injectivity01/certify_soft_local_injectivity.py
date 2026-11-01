"""Exact connectivity and orientation witnesses for the posterior-map attractor.

The topological implication is stated in the manuscript, not formalized here.
The previous trapping and periodic-point certificates remain dependencies.
"""
import argparse,hashlib,itertools,json,time
from collections import deque
from fractions import Fraction as Q
from pathlib import Path
import numpy as np

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def box(s):
    a,b=map(Q,s.strip('[]').split(','));assert a<=b;return a,b
def plus(a,b):return a[0]+b[0],a[1]+b[1]
def minus(a):return -a[1],-a[0]
def times(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def determinant(A):
    result=(Q(0),Q(0))
    for p in itertools.permutations(range(3)):
        term=(Q(1),Q(1))
        for i in range(3):term=times(term,A[i][p[i]])
        if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2:term=minus(term)
        result=plus(result,term)
    return result
def coordinates(key,shape):return key//(shape[1]*shape[2]),key//shape[2]%shape[1],key%shape[2]

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    with np.load(args.trap) as f:keys=list(map(int,f['occupied']));shape=list(map(int,f['shape']));origin=list(map(float,f['origin']));step=float(f['step']);params=f['params'].tolist();points=f['points'].tolist()
    assert keys==sorted(set(keys)) and len(keys)==35784 and step==1/256
    index={k:i for i,k in enumerate(keys)};parent=np.full(len(keys),-2,np.int32);parent[0]=-1;order=[];queue=deque([0]);depth=np.zeros(len(keys),np.int32);edges=0
    while queue:
        pos=queue.popleft();order.append(pos);key=keys[pos];coord=coordinates(key,shape)
        for axis,stride in enumerate([shape[1]*shape[2],shape[2],1]):
            for sign in [-1,1]:
                if not 0<=coord[axis]+sign<shape[axis]:continue
                target=index.get(key+sign*stride)
                if target is None:continue
                edges+=1
                if parent[target]==-2:parent[target]=pos;depth[target]=depth[pos]+1;queue.append(target)
    assert len(order)==len(keys) and np.all(parent>=-1)
    np.savez_compressed(args.out/'spanning-tree.npz',visit_order=np.array(order,np.uint32),parent=parent)
    sources={};rows=[]
    for label,path in [('covariance',args.cycles),('automatic_derivative',args.control)]:
        report=json.loads(path.read_text(encoding='utf-8'));assert report['status']=='passed' and report['trap_sha256']==sha(args.trap);sources[label+'_report_sha256']=sha(path)
        for ident,sign in [(0,1),(5,-1)]:
            case=next(c for c in report['cases'] if c['id']==ident);assert case['period']==(1 if ident==0 else 5);d=determinant([[box(v) for v in row] for row in case['monodromy']]);assert (d[0]>0 if sign==1 else d[1]<0)
            # The enclosing periodic-point proof already establishes existence
            # and membership. Check its saved strict root-box inclusion again.
            outer=list(map(box,case['box']));inner=list(map(box,case['krawczyk']));assert all(a[0]<b[0]<=b[1]<a[1] for a,b in zip(outer,inner));assert box(case['contraction_bound'])[1]<1
            rows.append(dict(engine=label,cycle_id=ident,period=case['period'],determinant_interval_exact=list(map(str,d)),determinant_approximate=[float(v) for v in d],sign=sign))
    primary=json.loads(args.cycles.read_text(encoding='utf-8'));assert primary['params']==params and primary['points']==points
    result=dict(status='passed',trap_sha256=sha(args.trap),cycle_sources=sources,cubes=len(keys),grid_shape=shape,origin_exact=list(map(lambda v:str(Q(v)),origin)),step_exact=str(Q(step)),connected_components=1,adjacency='Closed cubes sharing a full square face; exact integer addresses.',undirected_face_edges=edges//2,tree_root_key=keys[0],tree_edges=len(keys)-1,tree_maximum_depth=int(depth.max()),tree_sha256=sha(args.out/'spanning-tree.npz'),orientation_witnesses=rows,seconds=time.time()-begin,source_sha256=sha(Path(__file__)),scope='Exact face-connectedness of the fixed trapping union and opposite determinant signs from both located-cycle interval matrices. Conditional on the frozen trapping/root certificates and invariance of domain, the connected attracting set contains a point where F is not locally injective. No intrinsic collision, fold normal form, transitivity or historical uniqueness claim.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps({k:v for k,v in result.items() if k not in ['orientation_witnesses','origin_exact']}));print(json.dumps([{'engine':r['engine'],'cycle_id':r['cycle_id'],'determinant':r['determinant_approximate']} for r in rows]))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trap',type=Path,required=True);p.add_argument('--cycles',type=Path,required=True);p.add_argument('--control',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
