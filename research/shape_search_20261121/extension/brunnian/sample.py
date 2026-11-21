"""Sample three coexisting expanding solenoids in explicit disjoint braided tubes."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numba import njit


@njit(cache=True)
def iterate(degree, contraction, offset, ripple, shift, start, count):
    theta,x,y=start
    out=np.empty((count,3))
    exponent=0.
    for k in range(count+1024):
        x,y=contraction*x+offset*math.cos(theta),contraction*y+offset*math.sin(theta)
        if k>=1024:
            exponent+=math.log(degree*(1+ripple*math.cos(degree*theta)))
        theta=(degree*theta+ripple*math.sin(degree*theta)+shift)%(2*math.pi)
        if k>=1024:
            out[k-1024]=theta,x,y
    return out,exponent/count


def centers(theta,word,strand):
    letters=np.repeat(np.asarray(word),2)
    current=np.arange(3)
    start_slots=[]
    for letter in letters:
        start_slots.append(np.where(current==strand)[0][0])
        k=abs(letter)-1
        current[k:k+2]=current[k:k+2][::-1]
    assert np.array_equal(current,np.arange(3))
    phase=theta%(2*np.pi)*(len(letters)/(2*np.pi))
    sectors=np.minimum(phase.astype(int),len(letters)-1)
    t=phase-sectors
    # The quintic interpolation joins adjacent paths with matching first and second derivatives.
    turn=np.pi*t*t*t*(10+t*(-15+6*t))
    signs=np.sign(letters[sectors])
    left=np.abs(letters[sectors])-1
    slot=np.asarray(start_slots)[sectors]
    active=(slot==left)|(slot==left+1)
    half=slot-left-0.5
    u=np.where(active,left-0.5+half*np.cos(turn),slot-1.)
    v=np.where(active,signs*half*np.sin(turn),0.)
    return u,v


def physical(theta,x,y,row,strand):
    u,v=centers(theta,row["word"],strand)
    r=row["radius"]+u+row["width"]*x
    return np.stack((r*np.cos(theta),r*np.sin(theta),v+row["width"]*y),axis=-1)


def sample(row,n,out):
    clouds=[]
    rates=[]
    starts=[]
    for i,degree in enumerate(row["degrees"]):
        start=[.63891+.3131*i,.142-.13*i,-.17+.12*i]
        z,rate=iterate(degree,row["fiber_contraction"],row["fiber_offset"],row["base_ripple"],row["base_shift"],start,n)
        clouds.append(physical(*z.T,row,i))
        rates.append(rate)
        starts.append(start)
    points=np.concatenate(clouds)
    assert np.isfinite(points).all()
    path=out/f"case-{row['case']}.npz"
    np.savez_compressed(path,points=points,component=np.repeat(np.arange(3,dtype=np.int8),n))
    return dict(case=row["case"],sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                cloud_relative=f"extension/brunnian/{path.name}",shape=list(points.shape),
                steps_per_component=n,burn=1024,initial_states=starts,measured_base_exponents=rates)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("designs",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--cases",nargs="+",type=int)
    ap.add_argument("--steps",type=int,default=40000)
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    rows=json.loads(args.designs.read_text(encoding="utf-8"))["cases"]
    records=[sample(r,args.steps,args.out) for r in rows if args.cases is None or r["case"] in args.cases]
    args.manifest.parent.mkdir(parents=True,exist_ok=True)
    args.manifest.write_text(json.dumps(records,indent=2)+"\n",encoding="utf-8")
    print(f"{len(records)} link configurations with {3*args.steps} states each.")
