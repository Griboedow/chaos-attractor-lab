"""Compare rigid camera views of unaltered Yogi state clouds."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from yogi_map import trajectory
from render_twist import camera,density,colorize


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);d=np.load(args.source)
    if args.variant=='tanh':
        from yogi_tanh_map import trajectory as advance_cloud
    else:advance_cloud=trajectory
    font=ImageFont.truetype(findfont('DejaVu Sans'),22);large=ImageFont.truetype(findfont('DejaVu Sans'),30)
    records=[]
    for ident in args.ids:
        k=list(d['ids']).index(ident);p=d['parameters'][k];start=d['ends'][k]
        z=advance_cloud(p,start,10000,2000000);path=out/f'case-{ident}.npy';np.save(path,z)
        center=z.mean(0);ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
        if np.linalg.det(basis)<0:basis[:,-1]*=-1
        q=(z-center)@basis;shift=(q.max(0)+q.min(0))/2;q-=shift;center+=basis@shift
        halfwidth=float(np.max(np.linalg.norm(q,axis=1))*1.04)
        canvas=Image.new('RGB',(2000,1140),'#071320');draw=ImageDraw.Draw(canvas)
        draw.text((30,15),f'YOGI {args.variant.upper()} {ident} / EIGHT RIGID VIEWS',font=large,fill='#e8f0f5')
        draw.text((30,62),'2,000,000 actual states; identical units and spatial radius in every panel',font=font,fill='#9ab0c1')
        views=[]
        for i,(el,az) in enumerate((el,az) for el in [.25,.85] for az in [0.,np.pi/4,np.pi/2,3*np.pi/4]):
            projection=camera(az,el);h=density(q,projection,480,halfwidth);assert h.sum()==len(z);im,vmax=colorize(h)
            left=10+500*(i%4);top=110+500*(i//4);canvas.paste(im,(left,top))
            draw.text((left+10,top+455),f'{i}: az={az:.2f}, el={el:.2f}',font=font,fill='#a3b8ca')
            im.save(out/f'case-{ident}-view-{i}.png');np.savez_compressed(out/f'case-{ident}-view-{i}.npz',histogram=h)
            views.append(dict(index=i,azimuth=az,elevation=el,projection=projection.tolist(),log_density_vmax=vmax))
        canvas.save(out/f'case-{ident}-cameras.png')
        last=z[-2048:];returns=[float(np.max(np.linalg.norm(last[lag:]-last[:-lag],axis=1))) for lag in range(1,65)]
        records.append(dict(id=ident,variant=args.variant,source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),parameters=p.tolist(),
            start=start.tolist(),warmup=10000,observed=len(z),trajectory_file=path.name,trajectory_sha256=digest(path),
            center=center.tolist(),pca_basis_columns=basis.tolist(),halfwidth=halfwidth,views=views,late_return_residuals=returns,
            scope='Rigid orthographic camera comparison of the complete state, no axis rescaling or nonlinear deformation. All declared rendering points retained.'))
        (out/'analysis.json').write_text(json.dumps(dict(runs=records),indent=2)+'\n');print('Rendered',ident,'minimum late return',min(returns),flush=True)
    for name in ['preview_yogi.py','yogi_map.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    if args.variant=='tanh':(out/'yogi_tanh_map.py').write_bytes(Path(__file__).with_name('yogi_tanh_map.py').read_bytes())


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True)
    p.add_argument('--variant',choices=['sign','tanh'],default='sign');main(p.parse_args())
