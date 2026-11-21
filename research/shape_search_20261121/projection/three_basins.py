"""Map the separated finite-time regimes in Case 3014945 over initial states."""
import argparse
import json
from pathlib import Path
import numpy as np
import pyopencl as cl
from PIL import Image,ImageDraw
from third_growth import kernel_for,run_gpu
from reproduce import raw_parameters
from model import replay


def classify(r):
    r=np.asarray(r)
    if not np.isfinite(r).all():return 0
    if np.min(r)>.001:return 3
    if r[:,0].min()>.03 and r[:,1].max()<-.025:return 1
    if r[:,0].max()<-.008:return 2
    return 0


def terminal_class(points):
    x=points[:,0]
    if x.min()>-1.8 and x.max()<-.8:return 3
    if x.min()>-.25 and x.max()<.25:return 1
    return 0


def main(args):
    rec=json.loads(args.record.read_text());raw=raw_parameters(rec)
    dev=next(d for p in cl.get_platforms() for d in p.get_devices() if d.type&cl.device_type.GPU)
    ctx=cl.Context([dev]);q=cl.CommandQueue(ctx)
    kernel=cl.Program(ctx,Path(__file__).with_name('screen.cl').read_text()).build(options=['-DSECOND_TANGENT','-DTHIRD_TANGENT','-DSAVE_CLOUD']).screen
    axis=np.linspace(-2,2,args.side,dtype='f4')
    x,y=np.meshgrid(axis,axis)
    grid=np.column_stack([x.ravel(),y.ravel(),np.zeros(x.size)])
    rng=np.random.default_rng(3014945)
    volume=rng.uniform(-2,2,(args.starts,3)).astype('f4')
    starts=np.vstack((grid,volume)).astype('f4');all_rates=[];locations=[]
    for i in range(0,len(starts),2048):
        chunk=starts[i:i+2048]
        data=np.broadcast_to(raw[:,None,:],(len(raw),len(chunk),4)).astype('f4').copy()
        data[1,:,:3]=chunk
        n=len(chunk);stats=np.empty((20,n),'f4');tail=np.empty((128,n,4),'f4')
        a=cl.Buffer(ctx,cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR,hostbuf=data)
        b=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,stats.nbytes)
        c=cl.Buffer(ctx,cl.mem_flags.WRITE_ONLY,tail.nbytes)
        kernel(q,(n,),None,a,b,np.int32(n),np.int32(8192),np.int32(2048),c,np.int32(len(tail)))
        cl.enqueue_copy(q,stats,b).wait();cl.enqueue_copy(q,tail,c).wait()
        all_rates.append(np.stack([stats[:4].T,stats[12:16].T,stats[16:20].T],axis=2))
        locations.extend(terminal_class(z) for z in tail.transpose(1,0,2))
    rates=np.concatenate(all_rates)
    spectral=np.array([classify(v) for v in rates]);labels=np.array(locations)
    args.cache.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.cache,initial=starts,qr_windows=rates,classification=labels,spectral_class=spectral,grid_side=args.side)
    # An unbiased volume subsample provides a directly interpretable agreement rate.
    checks=[]
    for i in rng.choice(np.arange(len(grid),len(starts)),size=args.checks,replace=False):
        tail,rr,_,_,_=replay(raw,starts[i].astype(float),32768,16384,4,128,True)
        checks.append(dict(index=int(i-len(grid)),state=starts[i].tolist(),gpu=int(labels[i]),cpu=terminal_class(tail),gpu_spectral=int(spectral[i]),cpu_spectral=classify(rr),cpu_min=np.min(rr,axis=0).tolist(),cpu_max=np.max(rr,axis=0).tolist()))
    palette=[(74,73,82),(231,164,74),(170,90,152),(67,167,166)]
    field=np.asarray(palette,dtype='u1')[labels[:len(grid)]].reshape(args.side,args.side,3)[::-1]
    im=Image.new('RGB',(470,485),(249,247,243));dr=ImageDraw.Draw(im)
    im.paste(Image.fromarray(field).resize((384,384),Image.Resampling.NEAREST),(55,40))
    dr.text((55,15),'Initial-state basins in the z = 0 plane',fill='black')
    for j,t in enumerate([-2,-1,0,1,2]):
        dr.text((54+96*j,428),str(t),fill='black',anchor='mt')
        dr.text((48,424-96*j),str(t),fill='black',anchor='rm')
    dr.text((248,441),'Initial x',fill='black',anchor='mt')
    dr.text((5,223),'y',fill='black')
    for j,(i,name) in enumerate([(1,'Companion basin'),(3,'Primary basin'),(0,'Outside bands')]):
        pos=15+j*115;dr.rectangle((pos,465,pos+9,474),fill=palette[i]);dr.text((pos+14,463),name,fill='black')
    im.save(args.out/'basins-3014945.png',optimize=True)
    out=dict(case=rec['case'],grid_side=args.side,grid_plane='z=0',coordinate_range=[-2,2],random_seed=3014945,volume_samples=len(volume),gpu_burn=8192,gpu_windows=4,gpu_window=2048,classes=['outside bands','companion','unused','primary'],terminal_test='All last 128 x values in (-1.8,-0.8) for primary or (-0.25,0.25) for companion, enclosing the disjoint ranges of the long portrait samples.',spectral_classes=['unresolved','one positive','contracting','three positive'],gpu_volume_spectral_counts=np.bincount(spectral[len(grid):],minlength=4).tolist(),cpu_spectral_agree=sum(c['gpu_spectral']==c['cpu_spectral'] for c in checks),grid_counts=np.bincount(labels[:len(grid)],minlength=4).tolist(),volume_counts=np.bincount(labels[len(grid):],minlength=4).tolist(),cpu_random_checks=len(checks),cpu_agree=sum(c['gpu']==c['cpu'] for c in checks),cpu_burn=32768,cpu_window=16384,checks=checks,cache_relative='projection/'+args.cache.name)
    (args.out/'basins-3014945.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(out['volume_counts'],out['cpu_agree'],len(checks),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--record',type=Path,required=True);p.add_argument('--cache',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--side',type=int,default=193);p.add_argument('--starts',type=int,default=32768);p.add_argument('--checks',type=int,default=240)
    main(p.parse_args())
