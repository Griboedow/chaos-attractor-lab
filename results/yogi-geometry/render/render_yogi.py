"""Publication-sized Yogi portraits and a rigid rotation of fixed state clouds."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from yogi_map import trajectory
from render_twist import camera,density,colorize


def font(size):return ImageFont.truetype(findfont('DejaVu Sans'),size)


def layout(images,rows,scale):
    width=int(2800*scale);height=int(1680*scale)
    canvas=Image.new('RGB',(width,height),'#071320');draw=ImageDraw.Draw(canvas)
    def label(x,y,s,size,color='#e6edf5'):
        draw.text((int(x*scale),int(y*scale)),s,font=font(int(size*scale)),fill=color)
    label(65,30,'YOGI / GEOMETRY INSIDE AN OPTIMIZER',59)
    label(65,112,'One scalar quadratic objective. Three state variables. Actual iterates.',35,'#9eb5c9')
    for k,(im,row) in enumerate(zip(images,rows)):
        left=50+1400*k
        canvas.paste(im,(int(left*scale),int(270*scale)))
        label(left+20,188,row['title'],43)
        label(left+20,244,f"beta1 = {row['parameters'][0]:.6f}     beta2 = {row['parameters'][1]:.6f}     e = {row['parameters'][2]:.2g}",26,'#9eb5c9')
    label(65,1590,'Normalized state (x, m, v) / rigid cameras / brightness shows visit density',28,'#9eb5c9')
    label(65,1632,'Numerical candidates in a published map. Historical uniqueness remains unresolved.',26,'#9eb5c9')
    return canvas


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);preview=json.loads(args.preview.read_text())['runs']
    ids=[44,7590];titles=['PAIRED FANS / 44','SUNBURST / 7590'];azimuths=[np.pi/4,np.pi/4]
    rows=[];portraits=[];samples=[];histograms=[]
    for ident,title,az in zip(ids,titles,azimuths):
        i=list(source['ids']).index(ident);p=source['parameters'][i];start=source['ends'][i]
        z=trajectory(p,start,10000,8000000);path=out/f'case-{ident}-full.npy';np.save(path,z)
        basis=np.array(next(r for r in preview if r['id']==(7591 if ident==7590 else ident))['pca_basis_columns'])
        assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12 and np.linalg.det(basis)>0
        q=z@basis;shift=(q.max(0)+q.min(0))/2;q-=shift;center=basis@shift
        halfwidth=float(np.max(np.linalg.norm(q,axis=1))*1.055)
        projection=camera(az,.25);h=density(q,projection,1300,halfwidth);assert int(h.sum())==len(z)
        im,vmax=colorize(h);portraits.append(im);histograms.append(h)
        # Retain an explicit, contiguous, unthinned tail. Only its camera moves.
        sample=z[-1000000:].copy();np.savez_compressed(out/f'case-{ident}-animation.npz',states=sample,parameters=p)
        samples.append((sample-center)@basis)
        rows.append(dict(id=ident,title=title,parameters=p.tolist(),start=start.tolist(),warmup=10000,observed=len(z),
            source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),trajectory_file=path.name,
            trajectory_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),center=center.tolist(),basis=basis.tolist(),
            basis_source=str(args.preview),basis_source_case=(7591 if ident==7590 else ident),basis_source_sha256=hashlib.sha256(args.preview.read_bytes()).hexdigest(),
            halfwidth=halfwidth,azimuth=az,elevation=.25,projection=projection.tolist(),resolution=1300,log_density_vmax=vmax,
            gamma=.72,palette='magma',spatially_clipped=0,smoothing=False,
            animation=dict(points=len(sample),selection='Last 1000000 consecutive states of the declared trajectory; no thinning.',
                frames=72,frame_ms=100,resolution=650,log_density_vmax=float(np.log1p(180)))))
        print('Dense portrait',ident,flush=True)
    still=layout(portraits,rows,1.);still.save(out/'yogi-geometry.png')
    np.savez_compressed(out/'portrait-histograms.npz',ids=ids,histograms=histograms)
    frames=[];frame_counts=[]
    for frame,offset in enumerate(np.linspace(0,2*np.pi,72,endpoint=False)):
        images=[];counts=[]
        for sample,row in zip(samples,rows):
            h=density(sample,camera(row['azimuth']+offset,.25),650,row['halfwidth']);assert int(h.sum())==len(sample)
            im,_=colorize(h,vmax=np.log1p(180));images.append(im);counts.append(int(h.sum()))
        frames.append(layout(images,rows,.5));frame_counts.append(counts)
        if frame%18==0:frames[-1].save(out/f'rotation-frame-{frame:02d}.png');print('Rotation frame',frame,flush=True)
    frames[0].save(out/'yogi-rotation.gif',save_all=True,append_images=frames[1:],duration=100,loop=0,optimize=False)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,animation_frame_point_counts=frame_counts,
        scope='All images depict actual map states with rigid camera transformations only. Descriptive shape labels do not assert new dynamical classes.'),indent=2)+'\n')
    for name in ['render_yogi.py','yogi_map.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--preview',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
