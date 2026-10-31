"""Declared adaptive FP64 scout of the source ADSI three-state kinetics."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.stats import qmc
from rho_adsi import BASE,REFERENCE,SOURCE,TIME_RATE,field_jac,independent_field


def write(path,x):
    path.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')


class GPU:
    def __init__(self):
        self.device=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config)
        self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        self.program=cl.Program(self.context,Path(__file__).with_name('rho_adsi.cl').read_text()).build()
        self.advance=cl.Kernel(self.program,'advance');self.evaluate=cl.Kernel(self.program,'evaluate')

    def scan(self,p,initial,warm,duration,rtol=1e-9,wall=1200.,maxstep=.05):
        p=np.ascontiguousarray(p);n=len(p)
        states=np.zeros((n,8));states[:,:3]=initial;states[:,3:6]=np.array([1.,2.,3.])/np.sqrt(14.)
        clocks=np.zeros((n,2));clocks[:,1]=.001;logs=np.zeros((n,4));stats=np.zeros((n,8))
        stats[:,:6:2]=np.inf;stats[:,1:6:2]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[states,p,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        values=np.empty_like(states);vb=cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,values.nbytes)
        self.evaluate(self.queue,(n,),None,buffers[0],buffers[1],vb);cl.enqueue_copy(self.queue,values,vb).wait()
        worst=0.
        for i in range(n):
            f,j=field_jac(initial[i],p[i]);v=states[i,3:6]
            expected=np.r_[f,j@v,np.trace(j),0.]
            worst=max(worst,float(np.max(abs(values[i]-expected)/(1+abs(expected)))))
        assert worst<1e-11
        started=time.perf_counter();last=started;calls=0
        while True:
            self.advance(self.queue,(n,),None,*buffers,np.int32(128),np.float64(warm),np.float64(duration),
                np.float64(rtol),np.float64(rtol*.001),np.float64(maxstep),np.int32(10000000)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]);cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<warm+duration)
                if time.perf_counter()-last>25:
                    print('active',int(active.sum()),'wall',round(time.perf_counter()-started,1),flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-started>wall:
                    status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,initial=initial,states=states,clocks=clocks,logs=logs,stats=stats,status=status,counts=counts,
                    rate=windows.mean(axis=1),windows=windows,wall_seconds=time.perf_counter()-started,initial_field_error=worst)


def controls(gpu,out):
    records=[]
    for p in [REFERENCE,np.array([.3,3.,.2,4.]),np.array([3.,.3,5.,.1])]:
        x=np.array([.5,3.,1.]);f,j=field_jac(x,p)
        jf=np.column_stack([independent_field(x+1e-25j*np.eye(3)[i],p).imag/1e-25 for i in range(3)])
        algebra=float(max(np.max(abs(j-jf)),np.max(abs(f-independent_field(x,p)))))
        assert algebra<1e-10
        def augmented(t,z):
            f,j=field_jac(z[:3],p);return np.r_[f,j@z[3:]]
        v=np.array([1.,2.,3.])/np.sqrt(14.)
        cpu=solve_ivp(augmented,(0,2),np.r_[x,v],method='Radau',rtol=1e-11,atol=1e-13,max_step=.002)
        assert cpu.success
        result=gpu.scan(p[None],x[None],0,2,rtol=1e-11,maxstep=.005)
        norm=np.linalg.norm(cpu.y[3:,-1]);reference=cpu.y[3:,-1]/norm
        errors=dict(algebra=algebra,state=float(np.max(abs(result['states'][0,:3]-cpu.y[:3,-1]))),
                    tangent=float(np.max(abs(result['states'][0,3:6]-reference))),
                    rate=float(abs(result['rate'][0]-np.log(norm)/2)))
        assert result['status'][0]==0 and max(errors.values())<1e-7
        records.append(dict(parameters=p.tolist(),errors=errors))
    write(out/'independent-controls.json',records)
    # Reproduce the source kinetics before exploring any parameter changes.
    source=solve_ivp(lambda t,x:independent_field(x,REFERENCE),(0,500),[.5,3.,1.],
                    method='DOP853',rtol=1e-11,atol=1e-13,max_step=.025,dense_output=True)
    assert source.success
    grid=np.arange(250,500,.01);states=source.sol(grid).T
    np.savez_compressed(out/'source-control.npz',time_tau=grid,states=states,parameters=REFERENCE)
    write(out/'source-control.json',dict(parameters=REFERENCE.tolist(),warm=250,observation=250,
        minimum=states.min(axis=0).tolist(),maximum=states.max(axis=0).tolist(),
        scope='Published homogeneous kinetics; not a reproduction of the paper spatial PDE patterns.'))


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    u=qmc.Sobol(4,scramble=True,seed=2026103143).random_base2(args.power)
    lo=np.array([.05,.02,.02,.02]);hi=np.array([10.,20.,30.,30.])
    tuples=np.vstack([REFERENCE,np.exp(np.log(lo)+u*np.log(hi/lo))]);params=np.repeat(tuples,2,axis=0)
    starts=np.tile([[.5,3.,1.],[3.,3.,8.]],(len(tuples),1))
    protocol=dict(source=SOURCE,equation='7, with all diffusion terms zero',table='S2',base_k0_to_k10=BASE.tolist(),
        parameters=['GEF','GAP','inhibitor_speed','exchange_speed'],lower=lo.tolist(),upper=hi.tolist(),
        variation='inhibitor_speed multiplies k7,k8,k10; exchange_speed multiplies k5,k6; all other rates fixed.',
        concentration_unit='micromolar',time='tau=.025*seconds',physical_seconds_per_tau=1/TIME_RATE,
        seed=2026103143,power=args.power,tuples=len(tuples),lanes=len(params),warm=args.warm,observation=args.duration,
        rtol=1e-9,atol=1e-12,max_step=.05,method='FP64 adaptive Dormand-Prince 5(4), analytic tangent, normalization every accepted step',
        warning='Exploratory parameter family, not experimentally calibrated. No novelty or asymptotic chaos claim.',
        status_codes={'0':'completed','1':'magnitude cutoff 1e6','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved','5':'negative concentration below -1e-8'})
    write(out/'protocol.json',protocol)
    for name in ['rho_adsi.py','rho_adsi.cl','search_rho_adsi.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    np.savez_compressed(out/'declared-inputs.npz',parameters=params,initial=starts)
    gpu=GPU();controls(gpu,out)
    result=gpu.scan(params,starts,args.warm,args.duration,wall=args.wall)
    np.savez_compressed(out/'complete-or-budgeted.npz',**result)
    keep=(result['status']==0)&(result['windows'].min(axis=1)>.005)
    write(out/'analysis.json',dict(completed=int((result['status']==0).sum()),total=len(params),
        statuses={str(i):int((result['status']==i).sum()) for i in range(6)},retained_ids=np.flatnonzero(keep).tolist(),
        rates=result['rate'][keep].tolist(),retention='all four leading-rate windows > .005 per tau; complete lanes only',
        wall_seconds=result['wall_seconds'],initial_field_error=result['initial_field_error'],
        data_sha256=hashlib.sha256((out/'complete-or-budgeted.npz').read_bytes()).hexdigest()))
    print('complete',int((result['status']==0).sum()),'retained',int(keep.sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--power',type=int,default=10);p.add_argument('--warm',type=float,default=200.)
    p.add_argument('--duration',type=float,default=800.);p.add_argument('--wall',type=float,default=1200.)
    main(p.parse_args())
