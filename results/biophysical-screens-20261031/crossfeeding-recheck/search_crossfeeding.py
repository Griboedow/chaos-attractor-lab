"""Declared source-family chemostat survey with independent controls and budget."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pyopencl as cl
from scipy.integrate import solve_ivp
from scipy.stats import qmc
from crossfeeding import REFERENCE,REFERENCE_START,SOURCE,AUTHOR_COMMIT,PARAMETERS,field_jac,independent_field,balance_weights


def write(path,x):path.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf8')


class GPU:
    def __init__(self):
        self.device=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config)
        self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        self.program=cl.Program(self.context,Path(__file__).with_name('crossfeeding.cl').read_text()).build()
        self.advance=cl.Kernel(self.program,'advance');self.evaluate=cl.Kernel(self.program,'evaluate')

    def scan(self,p,initial,warm,duration,rtol=1e-9,wall=1200.,maxstep=.05):
        p=np.ascontiguousarray(p);n=len(p)
        states=np.zeros((n,16));states[:,:5]=initial;states[:,5:10]=np.arange(1.,6.)/np.sqrt(55.)
        clocks=np.zeros((n,2));clocks[:,1]=.001;logs=np.zeros((n,4));stats=np.zeros((n,16))
        stats[:,:10:2]=np.inf;stats[:,1:10:2]=-np.inf
        stats[:,12]=[balance_weights(k)@x for k,x in zip(p,initial)]
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[states,p,clocks,logs,stats,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        values=np.empty_like(states);vb=cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,values.nbytes)
        self.evaluate(self.queue,(n,),None,buffers[0],buffers[1],vb);cl.enqueue_copy(self.queue,values,vb).wait()
        worst=0.
        for i in range(n):
            f,j=field_jac(initial[i],p[i]);expected=np.r_[f,j@states[i,5:10],np.trace(j),np.zeros(5)]
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
    for x in [REFERENCE_START,np.array([.7,.4,30.,2.,5.]),np.array([.4,.7,2.,30.,5.])]:
        p=REFERENCE;f,j=field_jac(x,p)
        jf=np.column_stack([independent_field(x+1e-25j*np.eye(5)[i],p).imag/1e-25 for i in range(5)])
        algebra=float(max(np.max(abs(j-jf)),np.max(abs(f-independent_field(x,p)))))
        assert algebra<1e-10
        def augmented(t,z):
            f,j=field_jac(z[:5],p);return np.r_[f,j@z[5:]]
        v=np.arange(1.,6.)/np.sqrt(55.)
        cpu=solve_ivp(augmented,(0,2),np.r_[x,v],method='Radau',rtol=1e-11,atol=1e-13,max_step=.002)
        assert cpu.success
        result=gpu.scan(p[None],x[None],0,2,rtol=1e-11,maxstep=.005)
        norm=np.linalg.norm(cpu.y[5:,-1]);reference=cpu.y[5:,-1]/norm
        errors=dict(algebra=algebra,state=float(np.max(abs(result['states'][0,:5]-cpu.y[:5,-1]))),
                    tangent=float(np.max(abs(result['states'][0,5:10]-reference))),
                    rate=float(abs(result['rate'][0]-np.log(norm)/2)))
        assert result['status'][0]==0 and max(errors.values())<1e-7
        records.append(dict(parameters=p.tolist(),initial=x.tolist(),errors=errors))
    write(out/'independent-controls.json',records)


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    seed=2026103148
    low=np.array([.02,.001,.001,5.,.5,.3,.3,1.05,1.05,.3])
    high=np.array([.3,20.,20.,50.,2.,3.,3.,10.,10.,3.])
    u=qmc.Sobol(10,scramble=True,seed=seed).random_base2(args.power)
    design=np.exp(np.log(low)+u*np.log(high/low))
    tuples=np.tile(REFERENCE,(len(design)+1,1))
    tuples[1:,12:16]=design[:,:4];tuples[1:,6]=design[:,4]
    tuples[1:,1]=10*design[:,5];tuples[1:,7]=10*design[:,6]
    tuples[1:,5]=10*design[:,7];tuples[1:,11]=10*design[:,8]
    tuples[1:,2]=10*design[:,9];tuples[1:,8]=10/design[:,9]
    params=np.repeat(tuples,2,axis=0)
    starts=np.zeros((len(params),5));starts[:,:2]=np.tile([[.1,.11],[.11,.1]],(len(tuples),1));starts[:,2:]=params[:,13:16]
    write(out/'protocol.json',dict(source=SOURCE,author_commit=AUTHOR_COMMIT,parameters=PARAMETERS,
        seed=seed,power=args.power,design_parameters=['D','R1_in','R2_in','R3_in','mu2_max','k13/10','k23/10','q11/10','q22/10','k12/10=10/k21'],
        lower=low.tolist(),upper=high.tolist(),tuples=len(tuples),lanes=len(params),
        fixed='mu1_max=1; q13=q23=1; q12=q21=10. First tuple is the exact author source example.',
        initial='(.1,.11) and (.11,.1) populations; resources at their inlet values.',
        warm=args.warm,observation=args.duration,time='tau=D*t, author t in hours',
        rtol=1e-9,atol=1e-12,maxstep=.05,method='FP64 Dormand-Prince5(4), exact piecewise analytic tangent, normalization every accepted step.',
        error_norm='Componentwise physical state tolerance; tangent error components divided by atol+rtol*max(norm(v_old),norm(v_new)), because the tangent is normalized.',
        switching='Literal continuous min growth law; no smoothing or clipping. The tangent is continuous across transversal branch switches.',
        budget='Exact positive linear identity B_tau=B_in-B; numerical budget residual monitored at every accepted step.',
        retention='All four leading windows>.005/tau; exploratory weak threshold>.001/tau, complete lanes only.',
        status_codes={'0':'completed','1':'magnitude cutoff1e6','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved','5':'negative state below-1e-8','6':'population below1e-18: tiny-abundance unresolved','7':'scaled exact-budget residual>1e-6'},
        limitation='Exploratory source-parameter family, not a calibrated biological prediction or a novelty claim. Finite-time growth alone does not classify chaos.'))
    np.savez_compressed(out/'declared-inputs.npz',design=design,parameters=params,initial=starts)
    for name in ['crossfeeding.py','crossfeeding.cl','search_crossfeeding.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    gpu=GPU();controls(gpu,out);result=gpu.scan(params,starts,args.warm,args.duration,wall=args.wall)
    path=out/'complete-or-budgeted.npz';np.savez_compressed(path,**result)
    primary=(result['status']==0)&(result['windows'].min(axis=1)>.005)
    weak=(result['status']==0)&(result['windows'].min(axis=1)>.001)
    write(out/'analysis.json',dict(total=len(params),statuses={str(i):int((result['status']==i).sum()) for i in range(8)},
        primary_ids=np.flatnonzero(primary).tolist(),primary_rates=result['rate'][primary].tolist(),
        exploratory_weak_ids=np.flatnonzero(weak).tolist(),exploratory_weak_rates=result['rate'][weak].tolist(),
        wall_seconds=result['wall_seconds'],initial_field_error=result['initial_field_error'],
        maximum_scaled_budget_error=float(result['stats'][:,11].max()),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    print('completed',int((result['status']==0).sum()),'primary',int(primary.sum()),'weak',int(weak.sum()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--power',type=int,default=12);p.add_argument('--warm',type=float,default=200.)
    p.add_argument('--duration',type=float,default=1000.);p.add_argument('--wall',type=float,default=1200.)
    main(p.parse_args())
