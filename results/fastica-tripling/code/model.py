"""Exact homogeneous fourth-cumulant model of a structured 720-row data set."""
from pathlib import Path
import hashlib,json,math,time
import numpy as np
from numba import njit,prange

EXP=np.array([[4,0,0,0],[0,4,0,0],[2,2,0,0],
 [3,0,1,0],[0,3,0,1],
 [2,0,2,0],[0,2,0,2],[2,0,0,2],[0,2,2,0],
 [1,0,3,0],[0,1,0,3],
 [0,0,4,0],[0,0,0,4],[0,0,2,2]],dtype=np.int64)

def coefficients(e):
    return np.array([1,1,-6,24*e,24*e,72*e**2,72*e**2,-36*e**2,-36*e**2,72*e**3,72*e**3,-18*e**4,-18*e**4,-216*e**4],float)

def dataset(e):
    """Literal uniform observations; the analytic model uses exact radicals."""
    radial=[1.,1.,1.,1.,np.sqrt(6.)];noise=[np.sqrt(3.),-np.sqrt(3.),0.,0.,0.,0.];sigma=np.sqrt(1-6*e*e);rows=[]
    for radius in radial:
        for axis,sign in [(0,1),(0,-1),(1,1),(1,-1)]:
            xy=np.zeros(2);xy[axis]=sign*radius;shift=e*(radius*radius-4)*xy
            for a in noise:
                for b in noise:rows.append([*xy,shift[0]+sigma*a,shift[1]+sigma*b])
    return np.array(rows)

@njit(cache=True)
def update(w,coef):
    powers=np.ones((4,5))
    for j in range(4):
        for k in range(1,5):powers[j,k]=powers[j,k-1]*w[j]
    h=np.zeros(4);H=np.zeros((4,4))
    for t in range(len(coef)):
        for j in range(4):
            if EXP[t,j]==0:continue
            term=coef[t]*EXP[t,j]/4
            for k in range(4):term*=powers[k,EXP[t,k]-(k==j)]
            h[j]+=term
            for k in range(4):
                n=EXP[t,k]-(k==j)
                if n==0:continue
                term=coef[t]*EXP[t,j]*n/4
                for l in range(4):term*=powers[l,EXP[t,l]-(l==j)-(l==k)]
                H[j,k]+=term
    hn=np.sqrt(np.dot(h,h))
    if hn<1e-200:return w,H,hn
    f=h/hn;J=(H-np.outer(f,f@H))/hn
    return f,J,hn

@njit(cache=True)
def integrate(w,coef,warmup,steps,window,stride):
    w=w/np.sqrt(np.dot(w,w));a=w.copy();a[3]+=1 if w[3]>=0 else -1;a/=np.sqrt(np.dot(a,a));Q=np.eye(4)[:,:3].copy()
    for j in range(3):Q[:,j]-=2*a*a[j]
    rates=np.zeros((steps//window,3));cloud=np.empty((steps//stride,4));mn=1e100;diagmin=1e100;status=0;completed=0
    for it in range(warmup+steps):
        w,J,norm=update(w,coef);mn=min(mn,norm)
        if norm<1e-15 or not np.isfinite(norm):status=1;break
        V=J@Q
        for j in range(3):
            V[:,j]-=w*np.dot(w,V[:,j])
            for k in range(j):V[:,j]-=Q[:,k]*np.dot(Q[:,k],V[:,j])
            gain=np.sqrt(np.dot(V[:,j],V[:,j]));diagmin=min(diagmin,gain)
            if gain<1e-200 or not np.isfinite(gain):status=2;break
            Q[:,j]=V[:,j]/gain
            if it>=warmup:rates[(it-warmup)//window,j]+=math.log(gain)/window
        if status:break
        if it>=warmup and (it-warmup)%stride==0:cloud[(it-warmup)//stride]=w
        completed=it+1
    return rates,cloud,mn,diagmin,w,status,completed

@njit(parallel=True,cache=True)
def survey(coeffs,w):
    rates=np.empty((len(coeffs),5,3));minima=np.empty((len(coeffs),2));clouds=np.empty((len(coeffs),10000,4));status=np.zeros((len(coeffs),2),np.int64)
    for k in prange(len(coeffs)):
        r,c,m,d,last,s,n=integrate(w,coeffs[k],5000,100000,20000,10);rates[k]=r;clouds[k]=c;minima[k]=np.array([m,d]);status[k]=np.array([s,n])
    return rates,clouds,minima,status

if __name__=='__main__':
    root=Path(__file__).resolve().parent;out=root/'pilot-01';out.mkdir(exist_ok=False);eps=np.arange(1,81)/200;coeffs=np.array([coefficients(e) for e in eps]);w=np.array([.83,.29,-.13,.11]);protocol={'parameter':'epsilon = k/200, k=1..80; exact real model coefficients','start':w.tolist(),'warmup':5000,'observed':100000,'window':20000,'cloud_stride':10,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');t=time.time();rates,clouds,minima,status=survey(coeffs,w)
    np.savez_compressed(out/'survey.npz',epsilon=eps,rates=rates,clouds=clouds,minima=minima,status=status)
    rows=[]
    for k,e in enumerate(eps):
        row={'case':k+1,'epsilon_numerator':k+1,'epsilon_denominator':200,'mean':rates[k].mean(0).tolist(),'rates':rates[k].tolist(),'minimum_update_norm':float(minima[k,0]),'minimum_qr_diagonal':float(minima[k,1]),'status':int(status[k,0]),'completed':int(status[k,1])};rows.append(row)
    report={'wall_seconds':time.time()-t,'rows':rows,'positive_cases':[r['case'] for r in rows if r['status']==0 and min(v[0] for v in r['rates'])>.005]};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'seconds':report['wall_seconds'],'positive_cases':report['positive_cases'],'selected':rows[::10]}),flush=True)
