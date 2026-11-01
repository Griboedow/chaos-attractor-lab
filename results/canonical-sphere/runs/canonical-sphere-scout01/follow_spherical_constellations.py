"""Fresh-start full spectra for every flag in the paired sphere survey."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from numba import njit, prange, set_num_threads
from search_spherical_constellations import step


@njit(cache=True)
def qr_step(J,Q):
    R=J@Q; W=np.zeros((3,3)); logs=np.zeros(3)
    for a in range(3):
        v=R[:,a].copy()
        for repeat in range(2):
            for b in range(a): v-=np.dot(W[:,b],v)*W[:,b]
        length=np.linalg.norm(v)
        if length<1e-280 or not np.isfinite(length): return W,logs,1
        W[:,a]=v/length; logs[a]=np.log(length)
    return W,logs,0


@njit(cache=True,parallel=True)
def follow(params,points,counts,warm,count):
    n=len(params); rates=np.full((n,20,3),np.nan); states=np.full((n,100000,3),np.nan)
    status=np.zeros((n,2),np.int64); minima=np.full((n,2),np.inf); volume=np.zeros((n,20))
    for k in prange(n):
        beta=params[k,0]; p=points[k,:counts[k]]; x=np.array([.317,-.211,.419]); Q=np.eye(3); sums=np.zeros((20,3))
        for it in range(warm+count):
            x,J,gap,r,flag=step(x,beta,p); minima[k,0]=min(minima[k,0],gap); minima[k,1]=min(minima[k,1],r)
            if flag or not np.isfinite(x).all(): status[k,0]=flag if flag else 3; break
            Q,logs,flag=qr_step(J,Q)
            if flag: status[k,0]=4; break
            if it>=warm:
                w=(it-warm)//(count//20); sums[w]+=logs; volume[k,w]+=3*np.log(beta)+2*np.log(abs(1-1/r)); status[k,1]+=1
            if it>=warm+count-100000: states[k,it-(warm+count-100000)]=x
        if status[k,0]==0: rates[k]=sums/(count//20); volume[k]/=count//20
    return rates,states,status,minima,volume


def main(args):
    args.out.mkdir(parents=True,exist_ok=False); set_num_threads(4)
    # Independent LAPACK QR check with positive diagonal convention.
    rng=np.random.default_rng(2026110119); errors=[]
    for k in range(64):
        J=rng.normal(size=(3,3)); Q=np.linalg.qr(rng.normal(size=(3,3)))[0]
        W,logs,flag=qr_step(J,Q); U,R=np.linalg.qr(J@Q); U*=np.sign(np.diag(R))
        assert not flag; errors.append(max(float(np.max(abs(W-U))),float(np.max(abs(logs-np.log(abs(np.diag(R))))))))
    assert max(errors)<1e-11
    with np.load(args.source) as f:
        ids=np.flatnonzero(f['flags']); params=f['params'][ids]; points=f['points'][ids]; counts=f['counts'][ids]; kinds=f['kinds'][ids]
    begin=time.time(); rows=[]
    for offset in range(0,len(ids),4):
        s=slice(offset,offset+4); rates,states,status,minima,volume=follow(params[s],points[s],counts[s],30000,500000)
        for j,i in enumerate(ids[s]):
            k=offset+j; np.savez_compressed(args.out/f'case-{i}.npz',states=states[j],params=params[k],points=points[k,:counts[k]],rates=rates[j],status=status[j],minima=minima[j],kind=kinds[k],volume_rates=volume[j])
            row=dict(id=int(i),kind=int(kinds[k]),status=status[j].tolist(),spectrum=rates[j].mean(0).tolist(),minimum_leading_window=float(rates[j,:,0].min()),minimum_gap=float(minima[j,0]),minimum_radius=float(minima[j,1]))
            if status[j,0]==0:
                error=float(np.max(abs(rates[j].sum(1)-volume[j]))); assert error<1e-8
                ev=np.linalg.eigvalsh(np.cov(states[j].T)); row['variance_fractions']=(ev/ev.sum()).tolist(); row['logdet_qr_error']=error
                scale=max(float(np.linalg.norm(np.ptp(states[j],axis=0))),1.)
                residuals=[float(np.max(np.linalg.norm(states[j,-2048+lag:]-states[j,-2048:-lag],axis=1))) for lag in range(1,513)]
                periods=[lag+1 for lag,v in enumerate(residuals) if v<1e-8*scale]; row['detected_period']=min(periods) if periods else None; row['return_residuals']=residuals
            rows.append(row)
        print(json.dumps(dict(completed=min(offset+4,len(ids)),total=len(ids),seconds=time.time()-begin)),flush=True)
    report=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),warm=30000,observed=500000,windows=20,initial=[.317,-.211,.419],qr_check=dict(samples=64,max_error=max(errors)),cases=rows,scope='Local derivative spectra of a discontinuous nearest-point map; no asymptotic or historical novelty proof. All original flags followed.',seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'search_spherical_constellations.py']: (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--source',type=Path,required=True); p.add_argument('--out',type=Path,required=True); main(p.parse_args())
