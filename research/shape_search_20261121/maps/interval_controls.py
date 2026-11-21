"""Point, derivative, trigonometric and deliberately false certificate controls."""
import argparse
import json
from pathlib import Path
from flint import arb,arb_mat,ctx
from directed_return import rotate as directed_rotate,sincos,trig_constants
from verify_boxes import rotate as ball_rotate,frames,values,column
import numpy as np


def contains(lo,hi,value):
    # Arb midpoint/radius arithmetic encloses the true high-precision value.
    return arb(float(lo))<=value.lower() and arb(float(hi))>=value.upper()


def controls(box,symbolic=None):
    ctx.prec=160
    rng=np.random.default_rng(20261121)
    pi2,cs,eps=trig_constants()
    intervals=[]
    for c in np.r_[rng.uniform(-100,100,128),np.arange(-20,21)*np.pi/2]:
        for width in (0.,2.**-44,.01,1.1):
            intervals.append((c-width,c+width))
    values_checked=0
    for lo,hi in intervals:
        s,co=sincos((lo,hi),pi2,cs,eps)
        for point in np.linspace(lo,hi,9):
            ss,cc=arb(float(point)).sin_cos()
            if not contains(*s,ss) or not contains(*co,cc):
                raise ValueError(f"Failed trigonometric inclusion at {point}")
            values_checked+=2
    params=[arb(v) for v in box["params"]]
    dp=np.array([[v,v] for v in box["params"]])
    for x in rng.uniform(-1.,1.,(192,3)):
        raw=np.empty((3,4,2))
        raw[:,0,:]=x[:,None]
        raw[:,1:,:]=np.eye(3)[:,:,None]
        interval=directed_rotate(raw,dp,pi2,cs,eps,True)
        point,jac=ball_rotate([arb(float(v)) for v in x],params,arb_mat([[int(i==j) for j in range(3)] for i in range(3)]))
        for i in range(3):
            exact=[point[i,0]]+[jac[i,j] for j in range(3)]
            for j,v in enumerate(exact):
                if not contains(*interval[i,j],v):
                    raise ValueError("Independent interval state or tangent mismatch")
    o,_,inv=frames(box)
    y,_=ball_rotate(values(o[0]),params)
    uy=values(inv[0]*(y-o[0]))
    if not any(v>1 or v< -1 for v in uy):
        raise ValueError("The false order-one return passed the center test")
    out=dict(seed=20261121,trigonometric_intervals=len(intervals),point_trigonometric_values=values_checked,
        state_tangent_tests=192,state_or_derivative_coordinates=192*12,
        negative_control="Declaring the four-phase box a one-step return fails at its center.",
        negative_control_normalized_center=[float(v.mid()) for v in uy])
    if symbolic:
        packet,cover=symbolic
        o,m,inv=frames(packet)
        constants=[arb(v) for v in packet["params"][:4]]+[arb(0)]*6
        violations=[]
        for i,f in enumerate(cover["faces"]):
            c=[0.]*3;c[cover["axis"]]=f
            x=o[0]+m[0]*column(c)
            for _ in range(cover["map_steps"]):
                x,_=ball_rotate(values(x),constants)
            v=(inv[0]*(x-o[0]))[cover["axis"],0]
            if (cover["side"][i] and v<cover["faces"][-1]) or (not cover["side"][i] and v>cover["faces"][0]):
                violations.append(dict(face=i,normalized_image_midpoint=float(v.mid())))
        if not violations:
            raise ValueError("A contraction should violate the two-branch conditions")
        out["zero_twist_face_violations"]=violations
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("box",type=Path)
    p.add_argument("--symbolic-box",type=Path);p.add_argument("--symbolic-geometry",type=Path)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    symbolic=None
    if a.symbolic_box:
        symbolic=(json.loads(a.symbolic_box.read_text()),json.loads(a.symbolic_geometry.read_text())[0])
    r=controls(json.loads(a.box.read_text()),symbolic)
    a.output.write_text(json.dumps(r,indent=2)+"\n")
    print(r,flush=True)
