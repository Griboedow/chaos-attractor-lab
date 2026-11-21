"""Short state/tangent controls across a real OpenCL device and independent CPU algebra."""
import argparse
import json
from pathlib import Path
import numpy as np
import pyopencl as cl
from model import replay
from screen import parameters


def main(args):
    dev=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type & cl.device_type.GPU)
    ctx=cl.Context([dev]); q=cl.CommandQueue(ctx)
    src=Path(__file__).with_name('screen.cl').read_text()
    reports=[]
    for geom in ['single','linked','smooth']:
        linked=geom!='single'
        rng=np.random.default_rng(6161)
        data,_=parameters(351,args.size,geom)
        data[1,:,:3]=rng.uniform(-1.5,1.5,(args.size,3))
        metrics=np.empty((14 if linked else 12,args.size),np.float32)
        tails=np.empty((4,args.size,4),np.float32)
        read=cl.Buffer(ctx,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=data)
        stats=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,metrics.nbytes)
        write=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,tails.nbytes)
        flags=['-DSAVE_CLOUD']+(['-DLINKED'] if linked else [])+(['-DSMOOTH_LINK'] if geom=='smooth' else [])
        cl.Program(ctx,src).build(options=flags).screen(q,(args.size,),None,read,stats,np.int32(args.size),np.int32(0),np.int32(1),write,np.int32(4))
        cl.enqueue_copy(q,metrics,stats).wait()
        cl.enqueue_copy(q,tails,write).wait()
        states=[];exps=[]
        for j in range(args.size):
            raw=data[:,j,:].astype(float)
            cloud,rates,*_=replay(raw,raw[1,:3].copy(),0,1,4,4,False)
            states.append(np.max(np.abs(cloud-tails[:,j,:3]),axis=1))
            exps.append(np.abs(rates[:,0]-metrics[:4,j]))
        states=np.asarray(states);exps=np.asarray(exps)
        reports.append(dict(geometry=geom,size=args.size,seed=6161,parameters_batch=351,device=dev.name,driver=dev.driver_version,opencl=dev.version,numpy=np.__version__,pyopencl=cl.VERSION_TEXT,state_max_abs_by_step=states.max(axis=0).tolist(),state_median_abs_by_step=np.median(states,axis=0).tolist(),one_step_log_tangent_max_abs=float(exps[:,0].max()),all_four_step_log_tangent_max_abs=float(exps.max())))
    args.out.write_text(json.dumps(reports,indent=2)+'\n',encoding='utf8')
    print(json.dumps(reports))
    assert max(r['state_max_abs_by_step'][0] for r in reports)<2e-5
    assert max(r['one_step_log_tangent_max_abs'] for r in reports)<2e-3


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--size',type=int,default=256)
    main(p.parse_args())
