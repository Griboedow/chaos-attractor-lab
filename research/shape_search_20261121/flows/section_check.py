"""Independent finite-difference derivative of a first-return section near one period."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.integrate import solve_ivp,simpson
from numba import njit
from field import velocity_tangent,ring_tangent


@njit(cache=True)
def field_only(t,s,p,mode):
    y=np.empty(6);y[:3]=s;y[3:]=0.
    if mode==0:return velocity_tangent(y,p.reshape(3,9))[:3]
    return ring_tangent(y,p)[:3]


def run(case,source,inputs,scales):
    info=json.loads((source/f'case-{case}.json').read_text())
    z=np.load(source/f'case-{case}.npz')['points']
    scale=np.linalg.norm(np.std(z,axis=0))
    start=np.array(info['refined_start'])
    T=info['period'] if 'period' in info else info['final']['period']
    p=np.array(info.get('parameters',inputs.get(str(case),{}).get('parameters'))).ravel()
    mode=0 if info.get('family','ring')=='cartesian' else 1
    normal=field_only(0.,start,p,mode)
    normal=normal/np.linalg.norm(normal)
    unit=np.eye(3)[np.argmin(abs(normal))]
    first=np.cross(unit,normal);first/=np.linalg.norm(first)
    basis=np.column_stack((first,np.cross(normal,first)))
    def section(t,y):
        return np.dot(y-start,normal)
    section.direction=1.
    def F(t,y):
        return field_only(t,y,p,mode)
    rows=[]
    for method,rtol in [('RK45',2e-11),('DOP853',3e-12)]:
        def ret(loc):
            sol=solve_ivp(F,(0,1.025*T),loc,method=method,rtol=rtol,atol=rtol/100,events=section)
            if not sol.success:raise RuntimeError(sol.message)
            roots=sol.t_events[0];times=np.flatnonzero((roots>.975*T)&(roots<1.025*T))
            if len(times)==0:raise RuntimeError('No return close to the tested period')
            j=times[np.argmin(abs(roots[times]-T))]
            return basis.T@(sol.y_events[0][j]-start),float(roots[j])
        center,tc=ret(start)
        for rel in scales:
            h=rel*scale
            columns=[];times=[]
            for column in basis.T:
                above,ta=ret(start+h*column);below,tb=ret(start-h*column)
                columns.append((above-below)/(2*h));times.extend([ta,tb])
            derivative=np.column_stack(columns)
            vals=np.linalg.eigvals(derivative)
            rows.append(dict(method=method,rtol=rtol,relative_displacement=rel,section_jacobian=derivative.tolist(),
                             multipliers=[[float(a.real),float(a.imag)] for a in vals],
                             spectral_radius=float(max(abs(vals))),
                             fixed_point_condition_number_2=float(np.linalg.cond(derivative-np.eye(2))),
                             center_residual_rms=float(np.linalg.norm(center)/scale),
                             center_period=tc,max_time_deviation=float(max(abs(np.array(times)-T)))))
    spline=solve_ivp(F,(0,T),start,method='RK45',rtol=2e-11,atol=2e-13,dense_output=True)
    ts=np.linspace(0,T,8193);trajectory=spline.sol(ts)
    if mode==0:
        log_volume=-p.reshape(3,9)[:,0].sum()*T
    else:
        u=trajectory[0]+1j*trajectory[1];ss=abs(u)**2;h=2/(1+ss);n=int(p[0])
        div=2*p[1]-4*ss+2*p[2]*trajectory[2]-p[7]-(n-2)*p[6]*h**((n-2)/2)*np.real(u**n)/(1+ss)
        log_volume=float(simpson(div,x=ts))
    return dict(case=case,period=T,normal=normal.tolist(),basis=basis.tolist(),tests=rows,
                analytic_divergence_integral=log_volume,liouville_product=float(np.exp(log_volume)),
                full_variational_transverse_multipliers=info.get('transverse_multipliers'),
                continuation_multiplier=info.get('final',{}).get('dominant'))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sources',nargs='+',type=Path,required=True)
    ap.add_argument('--ids',nargs='+',type=int,required=True)
    ap.add_argument('--inputs',type=Path)
    ap.add_argument('--scales',type=float,nargs='+',default=[1e-4,3e-5,1e-5,3e-6])
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    for case in a.ids:
        source=next(p for p in a.sources if (p/f'case-{case}.json').exists())
        inputs=json.loads(a.inputs.read_text()) if a.inputs else {}
        val=run(case,source,inputs,a.scales)
        (a.out/f'case-{case}.json').write_text(json.dumps(val,indent=2)+'\n')
        print(case,[(r['method'],r['spectral_radius']) for r in val['tests']],flush=True)


if __name__=='__main__':
    main()
