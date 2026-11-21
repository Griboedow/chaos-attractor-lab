"""Native section event sampler for Case 180844, with independent Cartesian event checks."""
import math
from pathlib import Path
import json
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[4]
RESULT = ROOT / 'results/shape_search_20261121/frontiers/flows'
CACHE = Path('D:/Codex/chaos-attractor-cache/shape-search-20261121')
RAW = CACHE / 'frontiers/flows'
OLD = ROOT / 'results/shape_search_20261121/extension/flows'
CASES = {c['identifier']:c for c in json.loads((OLD/'manifest.json').read_text(encoding='utf-8'))['cases']}
STABLE = CASES['stable']
P = np.array(STABLE['parameters'])
assert P[0] == 7

@njit(inline='always')
def velocity(x,y,z,p):
    u=complex(x,y)
    u2=u*u
    u3=u*u2
    u6=u3*u3
    r=x*x+y*y
    f=2/(1+r)
    f52=f*f*math.sqrt(f)
    force=p[6]*f52*u6.conjugate()
    drive=p[11]*(u*u6).real-p[12]*(u*u6).imag
    c=p[1]-r+p[2]*z
    w=p[3]+p[4]*z+p[5]*r
    return x*c-y*w+force.real, x*w+y*c+force.imag, -p[7]*z+p[8]*(r-1)/(1+r)+p[9]*f52*f*drive

@njit(inline='always')
def step(x,y,z,h,p):
    a,b,c=velocity(x,y,z,p)
    d,e,f=velocity(x+h*.5*a,y+h*.5*b,z+h*.5*c,p)
    g,j,k=velocity(x+h*.5*d,y+h*.5*e,z+h*.5*f,p)
    l,m,n=velocity(x+h*g,y+h*j,z+h*k,p)
    return (x+h/6*(a+2*(d+g)+l),y+h/6*(b+2*(e+j)+m),z+h/6*(c+2*(f+k)+n))

@njit(inline='always')
def cubic(u,v,du,dv,s):
    return (1-s)*(1-s)*((1+2*s)*u+s*du)+s*s*((3-2*s)*v+(s-1)*dv)

@njit(cache=True,nogil=True)
def trace(initial,h,steps,p):
    x,y,z=initial
    events=np.empty((max(128,int(steps*h*2)),3))
    count=0
    for i in range(steps):
        xx,yy,zz=step(x,y,z,h,p)
        if y<0 and yy>=0:
            a,b,c=velocity(x,y,z,p)
            d,e,f=velocity(xx,yy,zz,p)
            lo,hi=0.,1.
            for _ in range(36):
                s=(lo+hi)/2
                if cubic(y,yy,h*b,h*e,s)<0:lo=s
                else:hi=s
            s=(lo+hi)/2
            events[count]=((i+s)*h,cubic(x,xx,h*a,h*d,s),cubic(z,zz,h*c,h*f,s))
            count+=1
            if count==len(events):raise ValueError('Too many section events')
        x,y,z=xx,yy,zz
    return events[:count].copy()

def packed(p=P):
    return np.r_[p,np.cos(p[10]),np.sin(p[10])]

def orbit_section():
    source=np.array(STABLE['nodes'][0]);T=STABLE['period'];points=[]
    import sys
    sys.path.insert(0,str(ROOT/'research/shape_search_20261121/flows'))
    from adaptive import polar_ring
    def rhs(t,u):return polar_ring(u,P)
    def event(t,u):return u[1]
    event.direction=1
    for k in range(7):
        phi=k*2*np.pi/7
        c,s=np.cos(phi),np.sin(phi)
        x,y,z=source
        state=(c*x-s*y,s*x+c*y,z)
        sol=solve_ivp(rhs,(0,T),state,events=event,rtol=2e-13,atol=2e-15,method='DOP853')
        assert sol.success
        events=sol.y_events[0][:,[0,2]]
        for a in events:points.append([k,*a])
    return np.array(points)

def raw_trace(path, initial, h, duration):
    if path.exists():return
    events=trace(np.array(initial),h,round(duration/h),packed())
    np.savez_compressed(path,events=events,initial=initial,h=h,duration=duration)
    print(path.name,len(events),flush=True)
