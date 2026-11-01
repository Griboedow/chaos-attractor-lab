"""Mechanism control: replace Yogi's sign by a declared smooth tanh switch.

Positive widths define a modified map, not the published optimizer. Every width
uses fresh declared starts; no smoothed result is a library reproduction.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from yogi_map import step_tangent


@njit
def smooth(z,p,width):
    x,m,v=z;b,c,e=p
    denom=v+x*x+e*e
    arg=(v-x*x)/(width*denom);switch=np.tanh(arg);slope=1-switch*switch
    sx=-2*x*(2*v+e*e)/(width*denom*denom)
    sv=(2*x*x+e*e)/(width*denom*denom)
    M=b*m+(1-b)*x;V=v-(1-c)*switch*x*x
    vx=-(1-c)*(2*x*switch+x*x*slope*sx)
    vv=1-(1-c)*x*x*slope*sv
    root=np.sqrt(V);d=root+e
    J=np.array([[1-(1-b)/d+M*vx/(2*root*d*d),-b/d,M*vv/(2*root*d*d)],
                [1-b,b,0.],[vx,0.,vv]])
    return np.array([x-M/d,M,V]),J,b*vv


@njit(parallel=True)
def survey(parameters,starts,widths,warm,observed,window,keep,Q):
    n=len(widths);windows=np.full((n,observed//window,3),np.nan);det_windows=np.full((n,observed//window),np.nan)
    tails=np.full((n,keep,3),np.nan);ends=starts.copy();status=np.zeros(n,dtype=np.int64);counts=np.zeros((n,2),dtype=np.int64)
    for i in prange(n):
        z=starts[i].copy();q=Q.copy();acc=np.zeros(3);logdet=0.;p=parameters[i];h=widths[i]
        for t in range(warm+observed):
            if h==0:
                z,q=step_tangent(z,q,p);det=p[0]
            else:
                z,J,det=smooth(z,p,h);q=J@q
            if not np.isfinite(z).all() or not np.isfinite(q).all() or det==0:
                status[i]=t+1;break
            logs=np.empty(3)
            for j in range(3):
                for _ in range(2):
                    for k in range(j):q[:,j]-=np.dot(q[:,j],q[:,k])*q[:,k]
                norm=np.linalg.norm(q[:,j])
                if norm==0 or not np.isfinite(norm):status[i]=-(t+1);break
                q[:,j]/=norm;logs[j]=np.log(norm)
            if status[i]!=0:break
            if t>=warm:
                acc+=logs;logdet+=np.log(abs(det));counts[i,0]+=int(det<0)
                if h>0:counts[i,1]+=int(abs(z[2]-z[0]*z[0])<h*(z[2]+z[0]*z[0]+p[2]*p[2]))
                k=t-warm+1
                if k%window==0:windows[i,k//window-1]=acc/window;det_windows[i,k//window-1]=logdet/window;acc[:]=0.;logdet=0.
                if k>observed-keep:tails[i,k-(observed-keep)-1]=z
        ends[i]=z
    return windows,det_windows,tails,ends,status,counts


def derivative_controls():
    rng=np.random.default_rng(2026110160);errors=[];det_errors=[]
    for _ in range(400):
        p=np.r_[rng.uniform(.1,.99,2),10**rng.uniform(-6,-1)];h=10**rng.uniform(-6,0)
        x=rng.choice([-1.,1.])*10**rng.uniform(-2,1);v=x*x*np.exp(h*rng.uniform(-2,2));m=rng.uniform(-1,1)
        z=np.array([x,m,v]);value,J,det=smooth(z,p,h)
        def independent(w):
            x,m,v=w;b,c,e=p;S=np.tanh((v-x*x)/(h*(v+x*x+e*e)))
            M=b*m+(1-b)*x;V=v-(1-c)*S*x*x
            return np.array([x-M/(np.sqrt(V)+e),M,V])
        complex_j=np.column_stack([np.imag(independent(z.astype(complex)+1e-25j*np.eye(3)[k]))/1e-25 for k in range(3)])
        errors.append(float(np.max(abs(J-complex_j))/(1+np.max(abs(complex_j)))))
        det_errors.append(float(abs(np.linalg.det(J)-det)/(1+abs(det))))
        assert np.max(abs(value-independent(z)))<1e-12*(1+np.max(abs(value)))
    assert max(errors)<1e-8 and max(det_errors)<1e-7
    return dict(cases=len(errors),complex_step_max_scaled_error=max(errors),determinant_max_scaled_error=max(det_errors))


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    d=np.load(args.source);width_values=np.r_[0.,np.logspace(-6,0,13)];parameters=[];starts=[];widths=[];labels=[]
    for ident in [670,125]:
        k=list(d['ids']).index(ident)
        for h in width_values:
            for kind,start in [('cold',np.array([.73,0.,0.])),('original_late_state',d['ends'][k])]:
                parameters.append(d['parameters'][k]);starts.append(start);widths.append(h);labels.append(dict(grid_id=ident,width=float(h),start_kind=kind))
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        declared=len(labels),widths=width_values.tolist(),grid_ids=[670,125],start_kinds=['cold','original_late_state'],
        smooth_switch='tanh((v-x^2)/(h*(v+x^2+e^2))) for h>0; original sign(v-x^2) at h=0.',
        warmup=100000,observed=1000000,window=50000,tail=4096,
        design='Every width starts independently from the stated cold state or a fixed original-map endpoint, not continuation from the preceding width.',
        scope='Mechanism comparison in an explicitly modified map. Positive widths are not published Yogi or an upstream-library control. Finite changes do not identify a bifurcation or prove equivalence.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for n in ['check_yogi_smoothing.py','yogi_map.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    controls=derivative_controls();(out/'derivative-controls.json').write_text(json.dumps(controls,indent=2)+'\n')
    Q=np.linalg.qr(np.array([[1.,2.,3.],[4.,5.,7.],[2.,6.,8.]]))[0]
    before=time.perf_counter();w,dets,tails,ends,status,counts=survey(np.array(parameters),np.array(starts),np.array(widths),100000,1000000,50000,4096,Q)
    np.savez_compressed(out/'histories.npz',parameters=parameters,starts=starts,widths=widths,windows=w,determinant_windows=dets,tails=tails,ends=ends,status=status,counts=counts)
    rows=[]
    for i,label in enumerate(labels):
        if status[i]:rows.append(dict(**label,status=int(status[i])));continue
        z=tails[i];period=next((lag for lag in range(1,129) if np.max(np.linalg.norm(z[lag:]-z[:-lag],axis=1))<1e-9*(1+np.max(abs(z)))),None)
        residual=float(np.max(abs(w[i].sum(1)-dets[i])));assert residual<1e-7
        row=dict(**label,status=0,spectrum=w[i].mean(0).tolist(),window_minima=w[i].min(0).tolist(),
            period_up_to128=period,determinant_sum_residual=residual,negative_determinant_fraction=float(counts[i,0]/1000000),
            post_update_switch_layer_fraction=float(counts[i,1]/1000000),tail_extrema=[z.min(0).tolist(),z.max(0).tolist()])
        rows.append(row);print(label['grid_id'],label['width'],label['start_kind'],row['spectrum'],period,flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,completed=int(np.sum(status==0)),declared=len(labels),wall_seconds=time.perf_counter()-before),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
