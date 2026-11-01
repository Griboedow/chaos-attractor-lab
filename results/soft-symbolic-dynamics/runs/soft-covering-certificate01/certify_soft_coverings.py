"""Exhaustive outward one-step window coverings inside the soft-map trap."""
import argparse,gzip,hashlib,itertools,json,time
from fractions import Fraction
from pathlib import Path
import numpy as np
from numba import njit
from soft_interval import evaluate,point,add,sub,mul


def exact_inverse(A):
    a=[list(row)+[Fraction(int(i==j)) for j in range(3)] for i,row in enumerate(A)]
    for j in range(3):
        pivot=next(i for i in range(j,3) if a[i][j]);a[j],a[pivot]=a[pivot],a[j];value=a[j][j];a[j]=[v/value for v in a[j]]
        for i in range(3):
            if i!=j:
                value=a[i][j];a[i]=[a[i][k]-value*a[j][k] for k in range(6)]
    return [row[3:] for row in a]


def outward(value):
    f=float(value);a=np.nextafter(f,-np.inf);b=np.nextafter(f,np.inf)
    assert Fraction(float(a))<=value<=Fraction(float(b));return [a,b]


@njit(cache=False)
def mm(A,B):
    result=np.zeros((3,3,2))
    for i in range(3):
        for j in range(3):
            for k in range(3):result[i,j]=add(result[i,j],mul(A[i,k],B[k,j]))
    return result


@njit(cache=False)
def mv(A,x):
    result=np.zeros((3,2))
    for i in range(3):
        for j in range(3):result[i]=add(result[i],mul(A[i,j],x[j]))
    return result


@njit(cache=False)
def enclose(ci,S,cj,T,mid,rad,beta,tau,points):
    q=np.empty((3,2));delta=np.empty((3,2))
    for j in range(3):q[j]=point(mid[j]);delta[j]=(-rad[j],rad[j])
    offset=mv(S,q);spread=mv(S,delta);C=np.empty((3,2));X=np.empty((3,2))
    for j in range(3):C[j]=add(point(ci[j]),offset[j]);X[j]=add(C[j],spread[j])
    _,J,axis,flag=evaluate(X,beta,tau,points)
    if flag:return np.zeros((3,2)),axis,flag
    center,_,_,flag=evaluate(C,beta,tau,points)
    if flag:return np.zeros((3,2)),axis,flag
    for j in range(3):center[j]=sub(center[j],point(cj[j]))
    Y=mv(T,center);offset=mv(mm(mm(T,J),S),delta)
    for j in range(3):Y[j]=add(Y[j],offset[j])
    return Y,axis,0


def accepted(Y,kind,orientation):
    if kind=='volume':return Y[0,1]<-1 or Y[0,0]>1 or bool(np.all(Y[1:,0]>-1) and np.all(Y[1:,1]<1))
    sign=(-orientation if kind=='left' else orientation)
    return (Y[0,0] if sign==1 else -Y[0,1])>1


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    with np.load(args.input) as f:raw={k:f[k] for k in ['params','points','centers','bases','radii','edges','orientation']}
    beta,tau=raw['params'];centers=raw['centers'];bases=raw['bases'];radii=raw['radii'];S=[];T=[];exact=[]
    for B,r in zip(bases,radii):
        matrix=[[Fraction(float(B[i,j]))*Fraction(float(r[j])) for j in range(3)] for i in range(3)];inverse=exact_inverse(matrix)
        S.append([[outward(v) for v in row] for row in matrix]);T.append([[outward(v) for v in row] for row in inverse]);exact.append(dict(matrix=[[str(v) for v in row] for row in matrix],inverse=[[str(v) for v in row] for row in inverse]))
    S=np.array(S);T=np.array(T);delta=np.tile([-1.,1.],(3,1));difference=np.array([sub(point(centers[1,j]),point(centers[0,j])) for j in range(3)]);a=mv(T[0],difference);b=mv(mm(T[0],S[1]),delta);separation=add(a[0],b[0]);assert separation[1]<-1 or separation[0]>1
    # A rectangular outer bound of N0 suffices here and is checked exactly.
    spread=mv(S[0],delta);hull=np.array([add(point(centers[0,j]),spread[j]) for j in range(3)])
    with np.load(args.trap) as f:assert np.array_equal(f['params'],raw['params']) and np.array_equal(f['points'],raw['points']);origin=f['origin'];shape=f['shape'];h=float(f['step']);occupied=set(f['occupied'].tolist())
    lower=[int((Fraction(float(hull[j,0]))-Fraction(float(origin[j])))//Fraction(h)) for j in range(3)];upper=[int((Fraction(float(hull[j,1]))-Fraction(float(origin[j])))//Fraction(h)) for j in range(3)];membership=[]
    for i in range(lower[0],upper[0]+1):
        for j in range(lower[1],upper[1]+1):
            for k in range(lower[2],upper[2]+1):
                assert 0<=i<shape[0] and 0<=j<shape[1] and 0<=k<shape[2];key=int((i*shape[1]+j)*shape[2]+k);assert key in occupied;membership.append(key)
    input_record=dict(**{k:v.tolist() for k,v in raw.items()},exact_windows=exact,proposal_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),trap_sha256=hashlib.sha256(args.trap.read_bytes()).hexdigest(),window0_hull=hull.tolist(),window0_hull_cube_keys=membership,window01_separation=list(separation))
    (args.out/'input.json').write_text(json.dumps(input_record,indent=2)+'\n',encoding='utf-8',newline='\n')
    sources={}
    for name in [Path(__file__).name,'soft_interval.py','raar_interval_scalar.py']:
        data=Path(__file__).with_name(name).read_bytes();(args.out/name).write_bytes(data);sources[name]=hashlib.sha256(data).hexdigest()
    records=[]
    for edgeid,(i,j) in enumerate(raw['edges']):
        orientation=int(raw['orientation'][edgeid]);assert orientation in [-1,1];scale=np.sum(np.maximum(abs(S[i,:,:,0]),abs(S[i,:,:,1])),axis=0)
        for kind in ['volume','left','right']:
            mid=np.zeros(3);rad=np.ones(3)
            if kind!='volume':mid[0]=-1 if kind=='left' else 1;rad[0]=0
            stack=[(mid,rad,'')];examined=accepted_count=0;maxdepth=0;minimum=1e100;last=time.time();path=args.out/f'edge-{edgeid}-{kind}.jsonl.gz'
            with gzip.open(path,'wt',encoding='utf-8',newline='\n') as handle:
                while stack and examined<args.max_cells:
                    mid,rad,address=stack.pop();examined+=1;maxdepth=max(maxdepth,len(address)//2)
                    Y,axis,flag=enclose(centers[i],S[i],centers[j],T[j],mid,rad,beta,tau,raw['points']);good=not flag and np.isfinite(Y).all() and accepted(Y,kind,orientation)
                    if good:
                        accepted_count+=1;minimum=min(minimum,axis);handle.write(json.dumps(dict(path=address,image=Y.tolist(),minimum_axis=axis),separators=(',',':'))+'\n')
                    else:
                        split=int(np.argmax(scale*rad));assert rad[split]>2.**-40,(edgeid,kind,address);child=rad.copy();child[split]*=.5
                        for sign in [-1,1]:
                            cm=mid.copy();cm[split]+=sign*child[split];stack.append((cm,child.copy(),address+str(split)+('0' if sign<0 else '1')))
                    if time.time()-last>15:print(json.dumps(dict(edge=edgeid,kind=kind,examined=examined,accepted=accepted_count,pending=len(stack),seconds=time.time()-begin)),flush=True);last=time.time()
            row=dict(edge=edgeid,source=int(i),target=int(j),kind=kind,status='passed' if not stack else 'incomplete',examined=examined,accepted=accepted_count,pending=len(stack),maximum_depth=maxdepth,minimum_axis=minimum,certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest());records.append(row);print(json.dumps(row),flush=True)
            if stack:(args.out/f'edge-{edgeid}-{kind}-pending.json').write_text(json.dumps([[m.tolist(),r.tolist(),p] for m,r,p in stack]),encoding='utf-8',newline='\n')
            (args.out/'progress.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8',newline='\n')
    result=dict(status='passed' if all(r['status']=='passed' for r in records) else 'incomplete',input_sha256=hashlib.sha256((args.out/'input.json').read_bytes()).hexdigest(),sources=sources,checks=records,seconds=time.time()-begin,scope='Exhaustive strict covering inequalities for the listed one-step graph. N0 lies in the previously certified trapping union; N0 and N1 are disjoint. Symbolic entropy and membership consequences require the complete graph and a coverage/source audit. No historical novelty, transitivity or equality with the sampled ribbon claim.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--trap',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--max-cells',type=int,default=200000);main(p.parse_args())
