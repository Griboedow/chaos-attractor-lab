"""Isolate stable periodic cycles observed in additional basins of selected maps."""
import argparse
import json
from pathlib import Path
from flint import arb,arb_mat,ctx
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.optimize import root
from periodic_points import iterate
from search import CACHE,RESULTS
from validate import write_json
from verify_boxes import column,matrix,rotate,values
from verify_return import tighten


def cyclic_certificate(p,hit):
    center=np.array(hit["orbit"])
    n=len(center)
    ctx.prec=192
    params=[arb(float(v)) for v in p]
    ident=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    size=3*n
    large=arb_mat([[int(i==j) for j in range(size)] for i in range(size)])
    for exponent in (30,36,41,45):
        rad=2.**-exponent
        full=arb_mat(size,size);res=arb_mat(size,1);W=ident
        phases=[[arb(float(v),rad) for v in u] for u in center]
        if any(not any(phases[i][k]<phases[j][k] or phases[j][k]<phases[i][k] for k in range(3))
               for i in range(n) for j in range(i)):
            continue
        for i,x in enumerate(center):
            j=(i+1)%n
            y,_=rotate([arb(float(v)) for v in x],params)
            _,J=rotate(phases[i],params,ident)
            W=J*W
            for h in range(3):
                res[3*i+h,0]=y[h,0]-arb(float(center[j,h]))
                for l in range(3):
                    full[3*i+h,3*i+l]+=J[h,l]
                full[3*i+h,3*j+h]-=1
        approx=arb_mat([[full[i,j].mid() for j in range(size)] for i in range(size)]).inv()
        pre=arb_mat([[approx[i,j].mid() for j in range(size)] for i in range(size)])
        if pre.det().contains(0):
            continue
        error=large-pre*full;delta=pre*res
        contraction=[sum(error[i,j].abs_upper() for j in range(size)) for i in range(size)]
        if not all(v<1 and v*rad+delta[i,0].abs_upper()<rad for i,v in enumerate(contraction)):
            continue
        trace=sum(W[i,i] for i in range(3))
        c2=sum(W[i,i]*W[j,j]-W[i,j]*W[j,i] for i in range(3) for j in range(i))
        sum_coeff=trace.abs_upper()+c2.abs_upper()+params[0]**(3*n)
        if sum_coeff<1:
            return dict(period=n,native_phase_centers=center.tolist(),native_cube_radius=rad,
                all_boxes_disjoint=True,root_contraction_upper_approximation=max(float(v) for v in contraction),
                characteristic_remainder_norm_upper_approximation=float(sum_coeff.upper()),
                conclusion="unique cyclic product-space root and all multipliers strictly inside the unit circle")
    return None


def certify(p,hit):
    period=hit["period"]
    def residual(x):
        y,J=iterate(x,p,period)
        return y-x,J-np.eye(3)
    center=root(residual,np.array(hit["orbit"][0]),jac=True,tol=1e-12).x
    _,J=iterate(center,p,period)
    if max(abs(np.linalg.eigvals(J)))>=1:
        return None
    P=solve_discrete_lyapunov(J.T,np.eye(3))
    shape=np.linalg.inv(np.linalg.cholesky(P).T)
    ctx.prec=192
    M=matrix(shape.tolist());iv=M.inv()
    o=column(center)
    param=[arb(float(v)) for v in p]
    spectral=np.linalg.norm(np.linalg.solve(shape,J@shape),ord=2)
    factor=float((1+spectral)/2)
    if factor>=1:
        return None
    for exponent in (28,33,38,42,46):
        rad=2.**-exponent
        c=o
        wide=o+M*column([arb(0,rad) for _ in range(3)])
        jac=M
        separated=True
        for j in range(period):
            c,_=rotate(values(c),param)
            wide,jac=rotate(values(wide),param,jac)
            wide=tighten(wide,c,jac,[rad]*3)
            if j+1<period:
                u=values(iv*(wide-o))
                if not any(v>rad or v< -rad for v in u):
                    separated=False;break
        if not separated:
            continue
        K=iv*jac
        S=arb_mat([[arb(factor)**2*int(i==j)-sum(K[k,i]*K[k,j] for k in range(3))
                    for j in range(3)] for i in range(3)])
        principal=[S[0,0],S[0,0]*S[1,1]-S[0,1]*S[1,0],S.det()]
        delta=sum(v.abs_upper()**2 for v in values(iv*(c-o))).sqrt()
        if all(d>0 for d in principal) and delta+factor*rad<rad:
            return dict(period=period,native_center=center.tolist(),linear_chart=shape.tolist(),
                normalized_euclidean_ball_radius=rad,uniform_derivative_norm_upper=factor,
                normalized_center_return_error_upper_approximation=float(delta.upper()),
                sylvester_lower_approximations=[float(v.lower()) for v in principal],
                all_earlier_phases_disjoint=True,conclusion="strict contraction on the closed ellipsoid; one attracting orbit of the stated least period")
    return cyclic_certificate(p,hit)


def build():
    records=[]
    for f in sorted((CACHE/"basins").glob("*.json")):
        data=json.loads(f.read_text())
        row=json.loads((RESULTS/"checks"/f.name).read_text())
        rings=[certify(np.array(row["params"]),h) for h in data["hits"]]
        records.append(dict(case=data["case"],initial_states=data["initial_states"],
            finite_negative_periodic_starts=data["negative_periodic_states"],
            observed_alternate_periods=[h["period"] for h in data["hits"]],
            attracting_cycle_certificates=[r for r in rings if r],
            example_starts=[h["start"] for h in data["hits"]]))
    write_json(RESULTS/"additional-basins.json",dict(precision_bits=192,cases=records))
    print([(r["case"],[h["period"] for h in r["attracting_cycle_certificates"]]) for r in records],flush=True)


if __name__=="__main__":
    build()
