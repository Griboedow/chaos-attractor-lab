"""Flow-normal quotient cocycle to resolve a nearly neutral motor spectrum.

Orthogonal projection removes the exact autonomous flow direction at every
QR boundary. Non-flow rates are compared with the full-spectrum calculation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import solve_ivp
from motor_rows import rhs
from check_motor_asymmetry import augmented


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);began=time.perf_counter()
    source=np.load(args.source);p=source['parameters'];q=source['states'][-1].copy();d=len(q);n=(d+1)//3
    f=rhs(0,q,p,n);frame=np.column_stack([f/np.linalg.norm(f),np.eye(d)])
    # SVD of the orthogonal projector gives a complete transverse frame.
    u,_,_=np.linalg.svd(np.eye(d)-np.outer(f,f)/np.dot(f,f))
    basis=np.column_stack([f/np.linalg.norm(f),u[:,:d-1]])
    logs=np.zeros(d-1);trace=0.;windows=[];states=[];times=[];previous=logs.copy();flow_error=0.;start_speed=None
    protocol=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),parameters=p.tolist(),
        initial=q.tolist(),warm=args.warm,duration=args.duration,rtol=args.rtol,max_step=args.max_step,QR_interval=.5,
        method='DOP853 full variational matrix; replace its first column by the exact unit flow before QR. The remaining columns represent the quotient tangent flow.',
        identity='sum(transverse_log_growth)+log(final_speed/observed_initial_speed)=integrated_divergence.',
        selection='Case2742 had two positive finite-time full-spectrum estimates but poorly resolved autonomous neutral direction. This control tests that specific ambiguity.',
        scope='Finite floating-point quotient spectrum, not a proof of hyperchaos or historical novelty.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_motor_transverse_spectrum.py','motor_rows.py','check_motor_asymmetry.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    warm_steps=int(round(2*args.warm));steps=int(round(2*args.duration));window_steps=steps//12
    assert steps%12==0
    for k in range(warm_steps+steps):
        old_speed=np.linalg.norm(rhs(0,q,p,n))
        sol=solve_ivp(lambda t,z:augmented(z,p,n),(0,.5),np.r_[q,basis.ravel(),0.],method='DOP853',
            rtol=args.rtol,atol=args.rtol*.01,max_step=args.max_step)
        assert sol.success
        q=sol.y[:d,-1];matrix=sol.y[d:-1,-1].reshape(d,d);f=rhs(0,q,p,n);speed=np.linalg.norm(f)
        flow_error=max(flow_error,float(np.max(abs(matrix[:,0]-f/old_speed))))
        matrix[:,0]=f/speed
        basis,R=np.linalg.qr(matrix)
        # Keep the first column positively oriented along the true flow.
        basis[:,0]=f/speed
        if k>=warm_steps:
            if start_speed is None:start_speed=old_speed
            logs+=np.log(abs(np.diag(R)[1:]));trace+=sol.y[-1,-1]
            states.append(q.copy());times.append((k-warm_steps+1)*.5)
            observed=k-warm_steps+1
            if observed%window_steps==0:
                windows.append(dict(end=observed*.5,spectrum=((logs-previous)/(window_steps*.5)).tolist()))
                previous=logs.copy()
                (out/'progress.json').write_text(json.dumps(dict(observed_tau=observed*.5,windows=windows,wall_seconds=time.perf_counter()-began),indent=2)+'\n')
                print(observed*.5,(logs/(observed*.5))[:3],flush=True)
    residual=abs(logs.sum()+np.log(speed/start_speed)-trace)/args.duration
    assert residual<1e-7 and flow_error<1e-7
    np.savez_compressed(out/'history.npz',time=times,states=states,parameters=p,final_basis=basis,log_sums=logs,integrated_trace=trace)
    result=dict(spectrum=(logs/args.duration).tolist(),windows=windows,flow_identity_max_error=flow_error,
        divergence_identity_residual_per_tau=float(residual),initial_speed=float(start_speed),final_speed=float(speed),
        samples=len(states),sampling='QR boundaries every.5; coarse diagnostics, not a dense curve rendering.',wall_seconds=time.perf_counter()-began)
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--warm',type=float,default=1000);p.add_argument('--duration',type=float,default=60000)
    p.add_argument('--rtol',type=float,default=2e-11);p.add_argument('--max-step',type=float,default=.01);main(p.parse_args())
