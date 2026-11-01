"""Numerical screen of published cyclic relaxed DR on three spatial circles.

T=T_12 o T_23 o T_31, evaluated right to left, equations (3b)-(3c) of
Dinh, Jansen and Luke, doi:10.1007/s10957-026-02996-2. A cycle is one step.
This is a known algorithm and a toy feasibility benchmark, not a physical model.
"""
import argparse, hashlib, json, platform, time
from pathlib import Path
import numpy as np
from numba import njit, prange, set_num_threads
from scipy.stats import qmc
from search_raar_circles import project


@njit(cache=True)
def geometry(p):
    centers=np.zeros((3,3));centers[1]=p[3:6];centers[2]=p[9:12]
    normals=np.zeros((3,3));normals[0,2]=1.
    normals[1]=np.array([0.,-np.sin(p[2]),np.cos(p[2])])
    normals[2]=np.array([np.sin(p[7])*np.cos(p[8]),np.sin(p[7])*np.sin(p[8]),np.cos(p[7])])
    return centers,normals,np.array([1.,p[1],p[6]])


@njit(cache=True)
def step(s,beta,centers,normals,radii):
    jac=np.eye(3);minimum=np.inf
    for i,j in ((2,0),(1,2),(0,1)):
        b,jb,lb,flag=project(s,centers[j],normals[j],radii[j])
        if flag:return s,jac,lb,1
        a,ja,la,flag=project(2*b-s,centers[i],normals[i],radii[i])
        minimum=min(minimum,la,lb)
        if flag:return s,jac,minimum,1
        derivative=beta*np.eye(3)+beta*ja@(2*jb-np.eye(3))+(1-2*beta)*jb
        jac=derivative@jac
        s=beta*s+beta*a+(1-2*beta)*b
    return s,jac,minimum,0


@njit(cache=True,parallel=True)
def screen(parameters,warm,count):
    n=len(parameters);rates=np.full((n,10),np.nan);tails=np.full((n,2048,3),np.nan)
    ends=np.zeros((n,3));status=np.zeros((n,2),np.int64);minaxis=np.full(n,np.inf);boundratios=np.zeros(n)
    for k in prange(n):
        p=parameters[k];c,norm,r=geometry(p);beta=p[0]
        s=np.array([.317,-.211,.419]);v=np.array([1.,.371,-.217]);v/=np.linalg.norm(v)
        constants=np.zeros(3)
        for ii,(i,j) in enumerate(((2,0),(1,2),(0,1))):
            constants[ii]=beta*(np.linalg.norm(c[i])+r[i])+abs(1-2*beta)*(np.linalg.norm(c[j])+r[j])
        bound=max(np.linalg.norm(s),(beta*beta*constants[0]+beta*constants[1]+constants[2])/(1-beta**3))
        sums=np.zeros(10)
        for t in range(warm+count):
            s,j,axis,flag=step(s,beta,c,norm,r);minaxis[k]=min(minaxis[k],axis)
            if flag or not np.isfinite(s).all():status[k,0]=1;break
            v=j@v;length=np.linalg.norm(v)
            if not np.isfinite(length) or length<1e-280:status[k,0]=2;break
            v/=length;boundratios[k]=max(boundratios[k],np.linalg.norm(s)/bound)
            if t>=warm:sums[(t-warm)//(count//10)]+=np.log(length);status[k,1]+=1
            if t>=warm+count-2048:tails[k,t-(warm+count-2048)]=s
        if status[k,0]==0:rates[k]=sums/(count//10)
        ends[k]=s
    return rates,tails,ends,status,minaxis,boundratios


def independent(s,p):
    c,n,r=geometry(p)
    def projection(x,i):
        d=x-c[i];d-=n[i]*np.dot(n[i],d)
        return c[i]+r[i]*d/np.linalg.norm(d)
    for i,j in ((2,0),(1,2),(0,1)):
        b=projection(s,j);reflected=2*b-s
        s=p[0]/2*(2*projection(reflected,i)-reflected+s)+(1-p[0])*b
    return s


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    u=qmc.Sobol(12,scramble=True,seed=2026110113).random_base2(args.power)
    p=np.column_stack([.6+.399*u[:,0],.4+2*u[:,1],.1+(np.pi/2-.1)*u[:,2],2.4*u[:,3:6]-1.2,
                       .4+2*u[:,6],np.arccos(2*u[:,7]-1),2*np.pi*u[:,8],2.4*u[:,9:12]-1.2])
    # Known nonchaotic control: all three constraints are the same unit circle.
    control=np.array([[.8,1,0,0,0,0,1,0,0,0,0,0.]])
    p=np.vstack([p,control]);rng=np.random.default_rng(2026110114);errors=[];forms=[]
    for row in p[:64]:
        x=rng.uniform(-2,2,3);c,n,r=geometry(row);y,j,_,flag=step(x,row[0],c,n,r);assert not flag
        forms.append(float(np.max(abs(y-independent(x.copy(),row)))))
        fd=np.column_stack([(independent(x+np.eye(3)[i]*1e-6,row)-independent(x-np.eye(3)[i]*1e-6,row))/2e-6 for i in range(3)])
        errors.append(float(np.max(abs(j-fd))/max(1.,np.max(abs(j)))))
    assert max(forms)<1e-12 and max(errors)<1e-5
    source={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in [Path(__file__).name,'search_raar_circles.py']}
    protocol=dict(seed=2026110113,count=len(p),warm=args.warm,observed=args.count,windows=10,initial=[.317,-.211,.419],
                  parameter_names=['beta','radius_B','tilt_B','cx_B','cy_B','cz_B','radius_C','polar_C','azimuth_C','cx_C','cy_C','cz_C'],
                  source_sha256=source,python=platform.python_version(),numpy=np.__version__,checks=dict(n=64,nested_formula=max(forms),relative_fd_jacobian=max(errors)),
                  source='https://doi.org/10.1007/s10957-026-02996-2',equations='3b and 3c; three circles; right-to-left pair order (3,1),(2,3),(1,2)',
                  scope='Numerical scout, no chaos/novelty certificate. Only full-cycle actual R3 states saved; no phase mixing. All trials retained outside Git.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    for name in source:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    parts=[[] for _ in range(6)];begin=time.time()
    for offset in range(0,len(p),64):
        d=screen(p[offset:offset+64],args.warm,args.count)
        for a,b in zip(parts,d):a.append(b)
        np.savez_compressed(args.out/f'block-{offset:04d}.npz',params=p[offset:offset+64],rates=d[0],tails=d[1],ends=d[2],status=d[3],minaxis=d[4],boundratios=d[5])
        print(json.dumps(dict(completed=min(offset+64,len(p)),total=len(p),seconds=time.time()-begin)),flush=True)
    rates,tails,ends,status,axes,ratios=[np.concatenate(v) for v in parts]
    flags=(status[:,0]==0)&(rates.min(1)>.005);assert ratios.max()<1+1e-12
    np.savez_compressed(args.out/'screen.npz',params=p,rates=rates,tails=tails,ends=ends,status=status,minaxis=axes,boundratios=ratios,flags=flags)
    report=dict(protocol,complete=int((status[:,0]==0).sum()),flagged=int(flags.sum()),flagged_ids=np.flatnonzero(flags).tolist(),max_bound_ratio=float(ratios.max()),seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--power',type=int,default=9);ap.add_argument('--warm',type=int,default=5000);ap.add_argument('--count',type=int,default=15000);main(ap.parse_args())
