"""FP64 exploratory integration of the literal logarithmic ice-cloud model."""
from pathlib import Path
import time
import numpy as np
import pyopencl as cl


class IceCloudGPU:
    def __init__(self):
        devices=[d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU and d.double_fp_config]
        if not devices:raise RuntimeError('FP64 GPU required')
        self.device=devices[0];self.context=cl.Context([self.device]);self.queue=cl.CommandQueue(self.context)
        program=cl.Program(self.context,Path(__file__).with_name('ice_cloud.cl').read_text()).build()
        self.advance=cl.Kernel(program,'advance');self.evaluate=cl.Kernel(program,'evaluate')

    def scan(self,parameters,initial,warm=200.,duration=800.,rtol=1e-9,max_step=.03,max_attempts=10000000,max_wall=1200.):
        p=np.asarray(parameters,dtype=float);q=np.asarray(initial,dtype=float);n=len(p)
        assert p.shape==(n,10) and q.shape==(n,3) and n>0 and np.isfinite(p).all() and np.isfinite(q).all()
        assert np.min(p)>0 and warm>=0 and min(duration,rtol,max_step,max_wall)>0
        pars=np.zeros((n,16));pars[:,:10]=p;states=np.zeros((n,16));states[:,:3]=q
        tangent=np.sin(np.arange(3)+1.);tangent/=np.linalg.norm(tangent);states[:,3:6]=tangent
        clocks=np.zeros((n,2));clocks[:,1]=1e-4;logs=np.zeros((n,4));ranges=np.zeros((n,16));ranges[:,:6:2]=np.inf;ranges[:,1:6:2]=-np.inf
        status=np.zeros(n,np.int32);counts=np.zeros((n,2),np.int32)
        arrays=[states,pars,clocks,logs,ranges,status,counts]
        buffers=[cl.Buffer(self.context,cl.mem_flags.READ_WRITE|cl.mem_flags.COPY_HOST_PTR,hostbuf=a) for a in arrays]
        began=time.perf_counter();last=began;calls=0
        while True:
            self.advance(self.queue,(n,),None,*buffers,np.int32(64),np.float64(warm),np.float64(duration),
                np.float64(rtol),np.float64(rtol*.001),np.float64(max_step),np.int32(max_attempts)).wait();calls+=1
            if calls%8==0:
                cl.enqueue_copy(self.queue,clocks,buffers[2]).wait();cl.enqueue_copy(self.queue,status,buffers[5]).wait()
                active=(status==0)&(clocks[:,0]<warm+duration)
                if time.perf_counter()-last>25:
                    print(f'{active.sum()} active; mean progress {clocks[:,0].mean()/(warm+duration):.4f}; wall {time.perf_counter()-began:.1f}s',flush=True);last=time.perf_counter()
                if not active.any():break
                if time.perf_counter()-began>max_wall:
                    status[active]=4;cl.enqueue_copy(self.queue,buffers[5],status).wait();break
        for a,b in zip(arrays,buffers):cl.enqueue_copy(self.queue,a,b).wait()
        windows=logs/(duration/4);windows[status!=0]=np.nan
        return dict(parameters=p,initial_states=q,initial_tangent=tangent,final_states=states[:,:3].copy(),
            final_tangents=states[:,3:6].copy(),lle=windows.mean(1),lle_windows=windows,status=status,ranges=ranges[:,:6],
            raw_log_sums=logs,elapsed=clocks[:,0].copy(),step_counts=counts,
            settings=dict(warm=warm,duration=duration,rtol=rtol,atol=rtol*.001,max_step=max_step,max_attempts=max_attempts,max_wall=max_wall),
            wall_seconds=time.perf_counter()-began,device=self.device.name,
            method='DP5(4) exact piecewise analytic tangent; tangent error uses its vector norm, independent of current component orientation.',
            status_codes={'0':'completed','1':'log concentration cutoff abs650 or saturation outside0..5','2':'numerical failure','3':'step budget unresolved','4':'wall budget unresolved'})

    def values(self,q,p,tangents):
        n=len(q);states=np.zeros((n,16));states[:,:3]=q;states[:,3:6]=tangents;pars=np.zeros((n,16));pars[:,:10]=p;out=np.empty_like(states)
        a=cl.Buffer(self.context,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=states)
        b=cl.Buffer(self.context,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=pars)
        c=cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,out.nbytes)
        self.evaluate(self.queue,(n,),None,a,b,c).wait();cl.enqueue_copy(self.queue,out,c).wait();return out[:,:6]
