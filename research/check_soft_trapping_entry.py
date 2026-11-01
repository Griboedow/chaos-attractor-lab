"""Outward finite-orbit entry check for the original start and nearby boxes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from soft_interval import evaluate,point,add,sub,mul,down,up


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    with np.load(args.source) as f:beta,tau=f['params'];points=f['points'];origin=f['origin'];shape=f['shape'];h=float(f['step']);occupied=set(f['occupied'].tolist())
    start=np.array([.317,-.211,.419]);rows=[]
    for radius in [0.,2.**-50,2.**-45,2.**-40,2.**-35]:
        box=np.array([add(point(x),(-radius,radius)) for x in start]);history=[];success=False;status='step_limit'
        for k in range(21):
            a=np.array([int(np.floor(down(down(box[j,0]-origin[j])/h))) for j in range(3)])
            b=np.array([int(np.floor(up(up(box[j,1]-origin[j])/h))) for j in range(3)])
            contained=(a>=0).all() and (b<shape).all() and int(np.prod(b-a+1))<1000000
            if contained:
                for ix in range(a[0],b[0]+1):
                    for iy in range(a[1],b[1]+1):
                        for iz in range(a[2],b[2]+1):
                            if int((ix*shape[1]+iy)*shape[2]+iz) not in occupied:contained=False
            history.append(dict(iteration=k,box=box.tolist(),contained=bool(contained)))
            if contained:success=True;status='entered';break
            image,J,axis,flag=evaluate(box,beta,tau,points)
            if flag:status='unresolved_axis';break
            mid=(box[:,0]+box[:,1])*.5;midbox=np.c_[mid,mid];center,_,_,flag=evaluate(midbox,beta,tau,points);assert not flag
            for i in range(3):
                value=tuple(center[i])
                for j in range(3):value=add(value,mul(J[i,j],sub(box[j],point(mid[j]))))
                image[i,0]=max(image[i,0],value[0]);image[i,1]=min(image[i,1],value[1]);assert image[i,0]<=image[i,1]
            box=image
        rows.append(dict(radius=radius,status=status,entered=success,history=history))
    files={}
    for name in [Path(__file__).name,'soft_interval.py','raar_interval_scalar.py']:
        data=Path(__file__).with_name(name).read_bytes();(args.out/name).write_bytes(data);files[name]=hashlib.sha256(data).hexdigest()
    result=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),initial=start.tolist(),cases=rows,sources=files,scope='Finite-orbit interval entry into the proposed union. Attraction additionally requires a completed trapping certificate; this check establishes no chaos or novelty.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps([dict(radius=r['radius'],status=r['status'],iterations=r['history'][-1]['iteration'],last_box=r['history'][-1]['box']) for r in rows]),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
