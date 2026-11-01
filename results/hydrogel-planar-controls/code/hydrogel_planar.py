"""Independent planar discretizations of Boiardi et al., arXiv:2604.09486v1.

The nonlinear solver approximates an invariant planar subspace, not the full
three-dimensional rod. Cubic Hermite beam elements supply an independent
linearized continuum control. No rotational drag or inertia is added.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals
from scipy.optimize import brentq
from numba import njit


def parameters(slenderness=50.0, eta=0.9, xi=5e-6, chi=100.0, plane="yz"):
    log = np.log(2 * slenderness / (1 + eta))
    z1 = 4*np.pi/(log + eta/(1+eta))
    z2 = 4*np.pi/(log + 1/(1+eta))
    z3 = 2*np.pi/(log - 0.5)
    b1 = np.pi*eta/(4*slenderness**4)
    b2 = np.pi*eta**3/(4*slenderness**4)
    if plane == "yz":
        b, zn, a, c = b1, z2, 1.0, chi
    else:
        b, zn, a, c = b2, z1, eta**-2, chi*eta**-3
    return dict(slenderness=slenderness, eta=eta, xi=xi, chi=chi,
                plane=plane, beta=b, normal_drag=zn, tangent_drag=z3,
                relaxation=a, drive=c)


def hermite_matrices(n):
    """Exact Gaussian quadrature for mass, curvature and slope bilinear forms."""
    h = 1/n
    size = 2*(n+1)
    mass, stiff, mixed = [np.zeros((size, size)) for _ in range(3)]
    q, w = np.polynomial.legendre.leggauss(4)
    for x, weight in zip((q+1)/2, w/2):
        shape = np.array([1-3*x*x+2*x**3, h*(x-2*x*x+x**3),
                          3*x*x-2*x**3, h*(-x*x+x**3)])
        d1 = np.array([-6*x+6*x*x, h*(1-4*x+3*x*x),
                       6*x-6*x*x, h*(-2*x+3*x*x)])/h
        d2 = np.array([-6+12*x, h*(-4+6*x),
                       6-12*x, h*(-2+6*x)])/h**2
        for el in range(n):
            ids = np.arange(2*el, 2*el+4)
            mass[np.ix_(ids, ids)] += weight*h*np.outer(shape, shape)
            stiff[np.ix_(ids, ids)] += weight*h*np.outer(d2, d2)
            mixed[np.ix_(ids, ids)] += weight*h*np.outer(d2, d1)
    return mass[2:,2:], stiff[2:,2:], mixed[2:,2:]


def hermite_eigenvalues(n, p, matrices=None):
    m, k, h = hermite_matrices(n) if matrices is None else matrices
    diffusivity = p["beta"]/(p["xi"]*p["normal_drag"])
    k, h = diffusivity*k, diffusivity*h
    a, c = p["relaxation"], p["drive"]
    z, eye = np.zeros_like(m), np.eye(len(m))
    left = np.block([[z, eye], [-a*k-c*h, -a*m-k]])
    right = np.block([[eye, z], [z, m]])
    values = eigvals(left, right)
    return values[np.argsort(values.real)[::-1]]


@njit(cache=True)
def angle_mass(theta, normal_drag, tangent_drag, xi):
    """Exact segment integral of translational RFT dissipation.

    Piecewise constant tangent; every segment's position derivative is linear
    in its local coordinate. Integrals of 1, f and f^2 are 1, 1/2 and 1/3.
    """
    n = len(theta)
    ds = 1.0/n
    sn, cs = np.sin(theta), np.cos(theta)
    tails = np.zeros((n+1, 3))
    for i in range(n-1, -1, -1):
        tails[i,0] = tails[i+1,0] + normal_drag*cs[i]**2 + tangent_drag*sn[i]**2
        tails[i,1] = tails[i+1,1] + (tangent_drag-normal_drag)*sn[i]*cs[i]
        tails[i,2] = tails[i+1,2] + normal_drag*sn[i]**2 + tangent_drag*cs[i]**2
    mass = np.empty((n,n))
    for a in range(n):
        for b in range(a,n):
            tail = tails[b+1]
            value = (cs[a]*cs[b]*tail[0] -
                     (cs[a]*sn[b]+sn[a]*cs[b])*tail[1] +
                     sn[a]*sn[b]*tail[2])
            value += normal_drag/3 if a == b else normal_drag*np.cos(theta[a]-theta[b])/2
            mass[a,b] = xi*ds**3*value
            mass[b,a] = mass[a,b]
    return mass


@njit(cache=True)
def angle_rhs_raw(state, beta, zn, zt, xi, a, c):
    n = len(state)//2
    theta, kappa = state[:n], state[n:]
    strain, orient = np.empty(n), np.empty(n)
    strain[0] = 2*n*theta[0] - kappa[0]
    orient[0] = theta[0]/2
    for i in range(1,n):
        strain[i] = n*(theta[i]-theta[i-1])-kappa[i]
        orient[i] = (theta[i]+theta[i-1])/2
    force = -beta*strain
    for i in range(n-1):
        force[i] += beta*strain[i+1]
    result = np.empty(2*n)
    result[:n] = np.linalg.solve(angle_mass(theta,zn,zt,xi), force)
    result[n:] = -a*kappa-c*np.sin(orient)
    return result


def rhs_parameters(p):
    return tuple(p[k] for k in ["beta","normal_drag","tangent_drag","xi","relaxation","drive"])


def angle_linear_matrix(n,p):
    derivative = np.eye(n)*n
    derivative[0,0] = 2*n
    for i in range(1,n):
        derivative[i,i-1] = -n
    weight = np.ones(n)/n
    weight[0] /= 2
    orient = np.eye(n)*0.5
    for i in range(1,n):
        orient[i,i-1] = 0.5
    b, zn, zt, xi, a, c = rhs_parameters(p)
    mass = angle_mass(np.zeros(n),zn,zt,xi)
    stiffness = b*derivative.T@(weight[:,None]*derivative)
    coupling = b*derivative.T*weight
    return np.block([[-np.linalg.solve(mass,stiffness),np.linalg.solve(mass,coupling)],
                     [-c*orient,-a*np.eye(n)]])


def leading_pair(values):
    lead = values[np.argmax(values.real)]
    return [float(lead.real),float(abs(lead.imag))]


def linear_controls():
    rows = []
    for n in (8,16,32,64):
        matrices = hermite_matrices(n)
        for plane in ("xz","yz"):
            base = parameters(plane=plane)
            def growth(chi,method):
                p = parameters(chi=chi,plane=plane)
                val = (hermite_eigenvalues(n,p,matrices) if method == "hermite"
                       else eigvals(angle_linear_matrix(n,p)))
                return max(val.real)
            thresholds = {method:float(brentq(lambda x:growth(x,method),25,70,xtol=1e-7))
                          for method in ("hermite","angle")}
            rows.append(dict(segments=n,plane=plane,critical_positive_chi=thresholds,
                             chi100_hermite=leading_pair(hermite_eigenvalues(n,base,matrices)),
                             chi100_angle=leading_pair(eigvals(angle_linear_matrix(n,base)))))
    # Passive cantilever exact first bending root and curvature relaxation.
    p = parameters(chi=0)
    exact = -p["beta"]/(p["xi"]*p["normal_drag"])*1.875104068711961**4
    passive = dict(exact_first_bending_eigenvalue=exact,curvature_relaxation=-1.0,
                   hermite=leading_pair(hermite_eigenvalues(32,p)),
                   angle=leading_pair(eigvals(angle_linear_matrix(64,p))))
    # An arbitrary small perturbation, not an eigenvector, checks all RHS signs.
    rng = np.random.default_rng(260409486)
    n = 12
    direction = rng.normal(size=2*n)
    p = parameters(chi=100)
    exact_jvp = angle_linear_matrix(n,p)@direction
    rhs_jvp = angle_rhs_raw(1e-7*direction,*rhs_parameters(p))/1e-7
    jvp_error = float(np.linalg.norm(rhs_jvp-exact_jvp)/np.linalg.norm(exact_jvp))
    # Independently assemble drag from local Jacobians and Gauss quadrature.
    theta = rng.normal(scale=0.4,size=n)
    mass = np.zeros((n,n))
    q,w = np.polynomial.legendre.leggauss(3)
    for i in range(n):
        normal = np.array([np.cos(theta[i]),-np.sin(theta[i])])
        tangent = np.array([np.sin(theta[i]),np.cos(theta[i])])
        drag = p["xi"]*(p["normal_drag"]*np.outer(normal,normal)+
                         p["tangent_drag"]*np.outer(tangent,tangent))
        for f,weight in zip((q+1)/2,w/2):
            jac = np.zeros((2,n))
            jac[:,:i] = np.array([np.cos(theta[:i]),-np.sin(theta[:i])])/n
            jac[:,i] = f*normal/n
            mass += weight/n*(jac.T@drag@jac)
    direct = angle_mass(theta,p["normal_drag"],p["tangent_drag"],p["xi"])
    drag_error = float(np.max(np.abs(direct-mass))/np.max(np.abs(mass)))
    assert jvp_error < 1e-10 and drag_error < 1e-12
    assert np.linalg.eigvalsh(direct).min() > 0
    return dict(rows=rows,passive=passive,linearization_relative_error=jvp_error,
                independent_drag_relative_error=drag_error,
                drag_min_eigenvalue=float(np.linalg.eigvalsh(direct).min()))


def run_planar(n,p,end=100.0,rtol=1e-7,method="BDF",output=None):
    positions = (np.arange(n)+0.5)/n
    state = np.r_[0.01*positions*(2-positions),np.zeros(n)]
    start = time.perf_counter()
    rhs = lambda t,y: angle_rhs_raw(y,*rhs_parameters(p))
    times = np.linspace(0,end,round(end*50)+1)
    def section(t,y):
        return np.mean(np.sin(y[:n]))
    section.direction = 1
    sol = solve_ivp(rhs,(0,end),state,method=method,rtol=rtol,atol=rtol*0.01,
                    t_eval=times,max_step=0.1,events=section)
    theta = sol.y[:n]
    tip = np.vstack([np.mean(np.sin(theta),axis=0),np.mean(np.cos(theta),axis=0)])
    late = sol.t >= end/2
    # Root-refined solver events avoid false return spread from sampled traces.
    event_mask = sol.t_events[0] >= end/2
    roots = sol.t_events[0][event_mask].tolist()
    states = sol.y_events[0][event_mask]
    periods = np.diff(roots)
    returns = [float(np.max(np.abs(b-a))) for a,b in zip(states[:-1],states[1:])]
    report = dict(parameters=p,segments=n,duration=end,rtol=rtol,atol=rtol*0.01,
                  method=method,success=sol.success,message=sol.message,
                  nfev=sol.nfev,njev=sol.njev,nlu=sol.nlu,elapsed_seconds=time.perf_counter()-start,
                  initial_condition="theta(s)=0.01*s*(2-s) at segment centers; kappa=0",
                  sample_interval=0.02,tip_range_late=[[float(v.min()),float(v.max())] for v in tip[:,late]],
                  angle_abs_max_late=float(np.max(np.abs(theta[:,late]))),
                  curvature_abs_max_late=float(np.max(np.abs(sol.y[n:,late]))),
                  crossing_times=roots,period_mean=float(np.mean(periods)) if len(periods) else None,
                  period_std=float(np.std(periods)) if len(periods) else None,
                  full_state_return_max_last10=max(returns[-10:]) if returns else None,
                  final_state=sol.y[:,-1].tolist())
    if output is not None:
        np.savez_compressed(output,t=sol.t,y=sol.y,tip=tip)
        report["external_array"] = str(output)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=["linear","run"])
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--segments",type=int,default=32)
    parser.add_argument("--plane",choices=["xz","yz"],default="yz")
    parser.add_argument("--chi",type=float,default=100)
    parser.add_argument("--slenderness",type=float,default=50)
    parser.add_argument("--eta",type=float,default=0.9)
    parser.add_argument("--xi",type=float,default=5e-6)
    parser.add_argument("--duration",type=float,default=100)
    parser.add_argument("--rtol",type=float,default=1e-7)
    parser.add_argument("--method",choices=["BDF","Radau"],default="BDF")
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    metadata = dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    if args.mode == "linear":
        report = linear_controls()
    else:
        report = run_planar(args.segments,parameters(args.slenderness,args.eta,args.xi,args.chi,args.plane),
                            args.duration,args.rtol,args.method,args.output/"trajectory.npz")
    report["environment"] = metadata
    (args.output/"report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
