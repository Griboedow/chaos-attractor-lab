from paths import OUTPUT, ROOT, RESEARCH
"""Full spectra and exact linear reduction, with independently coded map."""
from pathlib import Path
import sys,json,time
import numpy as np
from numba import njit,prange,set_num_threads
sys.path.insert(0,str(RESEARCH))
from twist_map import orthonormalize
from anderson_extended import step
set_num_threads(4)

@njit
def reduced(z,q,k):
    mean=z.mean();dev=z-mean;den=1+np.dot(dev,dev)
    grad=k/(3*den)-2*k*mean*dev/(den*den)
    out=np.empty_like(q);out[0]=q[1];out[1]=q[2];out[2]=grad@q
    return np.array([z[1],z[2],k*mean/den]),out

@njit(parallel=True)
def integrate(params,starts,modes,warm,observed,window,keep):
    w=np.zeros((len(params),observed//window,3));tails=np.zeros((len(params),keep,3));bounds=np.zeros(len(params))
    for i in prange(len(params)):
        z=starts[i].copy();q=np.eye(3);acc=np.zeros(3)
        for t in range(warm+observed):
            if modes[i]:z,q=reduced(z,q,params[i,0])
            else:
                qq=np.empty((3,3));zz=z.copy()
                for col in range(3):zz,qq[:,col]=step(z,q[:,col],params[i])
                z=zz;q=qq
            q,logs=orthonormalize(q)
            if t>=warm:
                acc+=logs;bounds[i]=max(bounds[i],np.max(np.abs(z)))
                if (t-warm+1)%window==0:w[i,(t-warm+1)//window-1]=acc/window;acc[:]=0
                offset=t-(warm+observed-keep)
                if offset>=0:tails[i,offset]=z
    return w,tails,bounds

def main():
    root=OUTPUT;d=np.load(root/'anderson-v2/data.npz');out=root/'anderson-spectra';out.mkdir(exist_ok=False)
    ids=np.array([41,54,165,449,555,629,682,1178,1616,2268,2332,2524]);p=d['parameters'][ids];starts=d['starts'][ids]
    # Three exact reduced linear cases are added as independently evaluated controls.
    lin=np.array([2268,2332,2524]);lp=d['parameters'][lin];lz=d['starts'][lin]*(lp[:,0,None]-1)/np.sqrt(lp[:,2,None])
    pp=np.vstack([p,lp]);zz=np.vstack([starts,lz]);modes=np.r_[np.zeros(len(ids),np.int64),np.ones(len(lin),np.int64)]
    begin=time.time();w,z,b=integrate(pp,zz,modes,100000,1000000,100000,200000)
    np.savez_compressed(out/'data.npz',ids=np.r_[ids,lin],parameters=pp,starts=zz,modes=modes,windows=w,tail=z,max_abs=b)
    rows=[{'id':int(i),'reduced':bool(mode),'spectrum':a.mean(0).tolist(),'min_leading':float(a[:,0].min()),'max_abs':float(bound)} for i,mode,a,bound in zip(np.r_[ids,lin],modes,w,b)]
    (out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'rows':rows},indent=2));print(json.dumps(rows),flush=True)

if __name__=='__main__':main()
