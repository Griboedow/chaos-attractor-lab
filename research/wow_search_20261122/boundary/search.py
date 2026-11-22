"""Find critically finite Newton maps and classify the basins seen at every kth step."""
from pathlib import Path
import argparse
import json
import math
import os
import numpy as np
from numba import njit
from scipy.stats import qmc

ROOT=Path(__file__).resolve().parents[3]
OUTPUT=ROOT/'results/wow_search_20261122/boundary'
RAW=Path(os.environ.get('WOW_BOUNDARY_CACHE','D:/Codex/chaos-attractor-cache/wow-search-20261122/boundary'))

@njit(cache=True)
def critical_iter(a,d,k):
    z=0j;v=0j
    for _ in range(k):
        if abs(z)>1e12:return complex(np.nan),complex(np.nan)
        zd=z**(d-2);zd1=zd*z
        A=(d-1)*zd1*z-a;B=d*zd1-a
        if abs(B)<1e-17:return complex(np.nan),complex(np.nan)
        derivative=d*(d-1)*zd*(z*zd1-a*z+a)/(B*B)
        v=derivative*v+(A-B)/(B*B);z=A/B
    return z,v

@njit(cache=True,nogil=True)
def candidates(initial,d,k):
    out=np.full(len(initial),complex(np.nan))
    for i in range(len(initial)):
        a=initial[i]
        for it in range(55):
            f,p=critical_iter(a,d,k)
            if not np.isfinite(f.real) or not np.isfinite(p.real) or abs(p)<1e-11 or abs(a)>100:break
            change=f/p
            if abs(f)<1e-10 and abs(change)<2e-11:
                if abs(a)>1e-4:out[i]=a-change
                break
            if abs(change)>2:change*=2/abs(change)
            a-=change
    return out

def orbit(a,d,k):
    z=0j;pts=[]
    for j in range(k):
        pts.append(z);z=((d-1)*z**d-a)/(d*z**(d-1)-a)
    return np.array(pts),z

def roots(d,a):
    coeff=np.zeros(d+1,dtype=complex);coeff[0]=1;coeff[-2]=-a;coeff[-1]=a
    r=np.roots(coeff);return r[np.argsort(np.angle(r)%(2*np.pi))]

@njit(cache=True,nogil=True)
def classify(z,d,a,k,roots,nmax):
    for i in range(nmax+1):
        if abs(z)<2e-7:return d+((-i)%k),i
        for j in range(d):
            if abs(z-roots[j])<2e-7:return j,i
        if abs(z)>1e25:return -1,i
        b=d*z**(d-1)-a
        if abs(b)<1e-24:return -1,i
        z=((d-1)*z**d-a)/b
    return -1,nmax

@njit(cache=True,nogil=True)
def plane(d,a,k,roots,xmin,xmax,ymin,ymax,nx,ny,nmax):
    labels=np.empty((ny,nx),np.int16);counts=np.empty((ny,nx),np.int16)
    for y in range(ny):
        im=ymax-(y+.5)*(ymax-ymin)/ny
        for x in range(nx):
            re=xmin+(x+.5)*(xmax-xmin)/nx
            labels[y,x],counts[y,x]=classify(complex(re,im),d,a,k,roots,nmax)
    return labels,counts

def measure_map(a,d,k):
    r=roots(d,a);cycle,err=orbit(a,d,k)
    if np.max(abs(cycle))>16 or min(abs(cycle[1:]))<.02 or abs(err)>1e-7:return None
    extent=max(1.6,min(5.,max(abs(r))*1.12))
    labels,counts=plane(d,a,k,r,-extent,extent,-extent,extent,91,91,140)
    mass=np.bincount(labels.ravel()+1,minlength=d+k+1)[1:]/labels.size
    cyclemass=mass[d:].sum()
    if min(mass[:d])<.015 or cyclemass<.027 or np.count_nonzero(mass[d:])<k:return None
    edgex=(labels[:,1:]!=labels[:,:-1]).mean()
    edgey=(labels[1:]!=labels[:-1]).mean()
    mix=.5*(edgex+edgey)
    tail=np.quantile(counts,.85)
    # Favor spatial complexity, a substantial non-root basin and several visible phase outcomes.
    score=2.6*mix+min(.2,cyclemass)*2+min(.035,np.quantile(mass[d:],.30))*3+min(36,tail)*.002
    return dict(d=d,k=k,a=[a.real,a.imag],cycle_fractions=mass[d:].tolist(),root_fractions=mass[:d].tolist(),
                cycle_fraction=cyclemass,mix=mix,score=score,extent=extent)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--levels',type=int,default=12);ap.add_argument('--max-seen',type=int,default=40)
    args=ap.parse_args();RAW.mkdir(parents=True,exist_ok=True)
    rows=[];counts=[]
    for d in (4,5,6):
        for k in (3,4,5,6,7,8,9):
            pts=qmc.Sobol(2,scramble=True,seed=83271+100*d+k).random_base2(args.levels)
            seeds=(d-1)+(pts[:,0]-.5)*2.4*d+1j*(pts[:,1]-.5)*2.4*d
            aa=candidates(seeds,d,k);groups={}
            for a in aa[np.isfinite(aa)]:
                if a.imag<0:a=a.conjugate()
                key=tuple(np.round([a.real,a.imag],6))
                groups.setdefault(key,[]).append(a)
            ordered=sorted(groups.values(),key=len,reverse=True)
            valid=0
            for cluster in ordered[:args.max_seen]:
                a=cluster[0];m=measure_map(a,d,k)
                if m:
                    m['starting_seed_hits']=len(cluster);rows.append(m);valid+=1
            counts.append(dict(d=d,k=k,seeds=len(seeds),converged=len(aa[np.isfinite(aa)]),
                               parameter_clusters_including_lower_period=len(groups),checked=min(len(groups),args.max_seen),passed=valid))
            (RAW/'candidate-search.json').write_text(json.dumps(dict(rows=rows,counts=counts),indent=2)+'\n')
            print(d,k,len(groups),valid,'top',sorted([(r['score'],r['d'],r['k'],r['a']) for r in rows],reverse=True)[:1],flush=True)

if __name__=='__main__':main()
