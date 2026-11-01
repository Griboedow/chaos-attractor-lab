"""Declared physical-sign survey of the published 2022 four-state wing model."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.stats import qmc
from wing_dsa import field_jac,source_equivalence,SOURCE


def write(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


class GPU:
    def __init__(self):
        self.device=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config)
        self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        self.program=cl.Program(self.context,Path(__file__).with_name('wing_dsa.cl').read_text()).build()
        self.advance=cl.Kernel(self.program,'advance');self.evaluate=cl.Kernel(self.program,'evaluate')

    def scan(self,p,initial,warm,duration,rtol=1e-9,maxstep=.03,wall=1200):
        p=np.ascontiguousarray(p);n=len(p);s=np.zeros((n,8));s[:,:4]=initial;s[:,4:]=np.array([1,2,3,4])/np.sqrt(30)
        clocks=np.zeros((n,2));clocks[:,1]=.001;logs=np.zeros((n,4));stats=np.zeros((n,8));stats[:,::2]=np.inf;stats[:,1::2]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32);arrays=[s,p,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        values=np.empty_like(s);vb=cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,values.nbytes)
        self.evaluate(self.queue,(n,),None,buffers[0],buffers[1],vb);cl.enqueue_copy(self.queue,values,vb).wait()
        error=0.
        for i in range(n):
            f,j=field_jac(initial[i],p[i]);expected=np.r_[f,j@s[i,4:]]
            error=max(error,float(np.max(abs(values[i]-expected)/(1+abs(expected)))))
        assert error<1e-11
        began=time.perf_counter();last=began;calls=0
        while True:
            self.advance(self.queue,(n,),None,*buffers,np.int32(128),np.float64(warm),np.float64(duration),np.float64(rtol),np.float64(rtol*.001),np.float64(maxstep),np.int32(10000000)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]);cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<warm+duration)
                if time.perf_counter()-last>25:
                    print('active',int(active.sum()),'wall',round(time.perf_counter()-began,1),flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-began>wall:
                    status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,initial=initial,states=s,clocks=clocks,logs=logs,stats=stats,status=status,counts=counts,
            windows=windows,rate=windows.mean(axis=1),wall_seconds=time.perf_counter()-began,initial_field_error=error)


def controls(gpu,out):
    records=[];v=np.array([1,2,3,4])/np.sqrt(30)
    for p in [np.array([.5,.5,4.,1.]),np.array([2.4,.5,.15,1.]),np.array([.05,.03,20.,1.])]:
        initial=np.array([.2,-.15,.1,-.12])
        def aug(t,q):return np.r_[source_equivalence(q[:4],p),field_jac(q[:4],p)[1]@q[4:]]
        cpu=solve_ivp(aug,(0,2),np.r_[initial,v],method='Radau',rtol=1e-11,atol=1e-13,max_step=.003)
        assert cpu.success
        gpu_result=gpu.scan(p[None],initial[None],0,2,rtol=1e-11,maxstep=.005)
        norm=np.linalg.norm(cpu.y[4:,-1]);errors=dict(state=float(np.max(abs(gpu_result['states'][0,:4]-cpu.y[:4,-1]))),
            tangent=float(np.max(abs(gpu_result['states'][0,4:]-cpu.y[4:,-1]/norm))),rate=float(abs(gpu_result['rate'][0]-np.log(norm)/2)))
        assert gpu_result['status'][0]==0 and max(errors.values())<1e-7
        records.append(dict(parameters=p.tolist(),errors=errors))
    write(out/'independent-controls.json',records)


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source_report=json.loads((args.control/'analysis.json').read_text());assert len(source_report['controls'])==3
    lower=np.array([.03,.01,.001]);upper=np.array([30.,.99,100.]);u=qmc.Sobol(3,scramble=True,seed=2026110103).random_base2(args.power)
    tuples=np.column_stack([np.exp(np.log(lower)+u*np.log(upper/lower)),np.ones(len(u))])
    eigen=np.array([max(np.linalg.eigvals(field_jac(np.zeros(4),p)[1]).real) for p in tuples]);mask=eigen>.001
    tuples=tuples[mask];parameters=np.repeat(tuples,2,axis=0)
    starts=np.tile([[.1,0,0,0],[1,.5,-.3,.2]],(len(tuples),1))
    protocol=dict(source=SOURCE,equations='4 and9; no clipping, forcing or added state',parameters=['r3/omega_n','kappa','mu/k','clock_factor'],
        lower=lower.tolist(),upper=upper.tolist(),clock_factor=1.,seed=2026110103,power=args.power,selection='Origin largest real eigenvalue>.001; does not search hidden attractors with a stable origin.',
        lanes=len(parameters),starts=starts[:2].tolist(),warm=500.,observation=2000.,rtol=1e-9,atol=1e-12,maxstep=.03,
        retention='Every one of four windows>.005; exploratory weak threshold>.001',
        time='tau=sqrt(k/I)*t; nonlinear drag sets the coordinate amplitude unit',
        source_control_sha256=hashlib.sha256((args.control/'analysis.json').read_bytes()).hexdigest(),
        limitations='Mathematical physical-sign parameter family, not fitted insect or robot predictions. The source quadratic-drag coefficient is not recovered; physical angular amplitudes are not claimed.',
        status_codes={'0':'complete','1':'magnitude cutoff1e7','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved'})
    write(out/'protocol.json',protocol);np.savez_compressed(out/'declared-inputs.npz',parameters=parameters,initial=starts,sobol_ids=np.repeat(np.flatnonzero(mask),2),start_id=np.tile([0,1],len(tuples)))
    for name in ['wing_dsa.py','wing_dsa.cl','search_wing_dsa.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    gpu=GPU();controls(gpu,out);result=gpu.scan(parameters,starts,500,2000)
    np.savez_compressed(out/'complete-or-budgeted.npz',**result)
    keep=(result['status']==0)&(result['windows'].min(axis=1)>.005)
    write(out/'analysis.json',dict(lanes=len(parameters),statuses={str(i):int((result['status']==i).sum()) for i in range(5)},retained_ids=np.flatnonzero(keep).tolist(),
        weak_ids=np.flatnonzero((result['status']==0)&(result['windows'].min(axis=1)>.001)).tolist(),wall_seconds=result['wall_seconds'],initial_field_error=result['initial_field_error'],sha256=hashlib.sha256((out/'complete-or-budgeted.npz').read_bytes()).hexdigest()))
    print('complete',int((result['status']==0).sum()),'retained',np.flatnonzero(keep).tolist(),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--control',type=Path,required=True);p.add_argument('--power',type=int,default=13);main(p.parse_args())
