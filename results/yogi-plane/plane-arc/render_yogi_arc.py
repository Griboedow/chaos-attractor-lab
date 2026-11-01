"""Render the G0670 state cloud with a rigid camera and uniform spatial units."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from matplotlib.font_manager import findfont
from yogi_map import trajectory
from render_twist import camera, colorize


def histogram(points, projection, width, height, units_per_pixel):
    xy = points @ projection.T
    ix = np.floor(xy[:,0]/units_per_pixel+width/2).astype(np.int32)
    iy = np.floor(xy[:,1]/units_per_pixel+height/2).astype(np.int32)
    assert np.all((ix>=0)&(ix<width)&(iy>=0)&(iy<height))
    h = np.bincount(iy*width+ix,minlength=width*height).reshape(height,width)
    assert h.sum() == len(points)
    return h


def compose(picture, parameters, scale):
    canvas=Image.new('RGB',(int(2400*scale),int(1260*scale)),'#071320')
    draw=ImageDraw.Draw(canvas)
    def text(x,y,label,size,color='#e8f0f5'):
        font=ImageFont.truetype(findfont('DejaVu Sans'),int(size*scale))
        draw.text((int(x*scale),int(y*scale)),label,font=font,fill=color)
    text(65,35,'YOGI / FEATHERED ARC',68)
    text(65,126,'One-variable quadratic optimization / full three-variable memory state',32,'#9ab0c1')
    canvas.paste(picture,(int(50*scale),int(230*scale)))
    b,c,e=parameters
    text(65,1085,f'G0670     beta1 = {b:.8f}     beta2 = {c:.8f}     e = {e:.8g}',29,'#b8c9d7')
    text(65,1140,'Actual iterates; orthographic view; equal spatial units on both image axes.',28,'#9ab0c1')
    text(65,1186,'Numerical chaos candidate in a published algorithm. Historical novelty remains unresolved.',27,'#9ab0c1')
    return canvas


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);i=list(source['ids']).index(670)
    p=source['parameters'][i];start=source['ends'][i]
    z=trajectory(p,start,10000,8000000)
    full=out/'case-670-full.npy';np.save(full,z)
    preview=json.loads(args.preview.read_text())['runs'];basis=np.array(next(r for r in preview if r['id']==670)['pca_basis_columns'])
    assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12 and np.linalg.det(basis)>0
    q=z@basis;shift=(q.max(0)+q.min(0))/2;q-=shift
    # A sphere containing all declared states bounds every rigid view. The
    # vertical bound uses the actual camera elevation and axis-aligned limits.
    radius=float(np.max(np.linalg.norm(q,axis=1)));ext=np.max(abs(q),axis=0)
    el=.25;vertical=float(np.sin(el)*np.hypot(ext[0],ext[1])+np.cos(el)*ext[2])
    width,height=2300,800
    unit=1.055*max(2*radius/width,2*vertical/height)
    projection=camera(0.,el);h=histogram(q,projection,width,height,unit)
    im,vmax=colorize(h);compose(im,p,1).save(out/'yogi-feathered-arc.png')
    np.savez_compressed(out/'portrait-histogram.npz',histogram=h)
    sample=z[-1000000:].copy();np.savez_compressed(out/'animation.npz',states=sample,parameters=p)
    sample=sample@basis-shift
    frames=[];counts=[]
    for frame,az in enumerate(np.linspace(0,2*np.pi,72,endpoint=False)):
        hh=histogram(sample,camera(az,el),width//2,height//2,unit*2)
        im,_=colorize(hh,vmax=np.log1p(180))
        frames.append(compose(im,p,.5));counts.append(int(hh.sum()))
        if frame%18==0:frames[-1].save(out/f'rotation-frame-{frame:02d}.png')
    frames[0].save(out/'yogi-feathered-arc.gif',save_all=True,append_images=frames[1:],duration=100,loop=0,optimize=False)
    last=z[-2048:];returns=[float(np.max(np.linalg.norm(last[lag:]-last[:-lag],axis=1))) for lag in range(1,513)]
    metadata=dict(grid_id=670,parameters=p.tolist(),start=start.tolist(),warmup=10000,observed=len(z),
        source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        preview=str(args.preview),preview_sha256=hashlib.sha256(args.preview.read_bytes()).hexdigest(),
        full_trajectory_file=full.name,full_trajectory_sha256=hashlib.sha256(full.read_bytes()).hexdigest(),
        basis=basis.tolist(),center=(basis@shift).tolist(),units_per_pixel=unit,
        resolution=[width,height],projection=projection.tolist(),elevation=el,azimuth=0.,
        density_vmax=vmax,gamma=.72,palette='magma',smoothing=False,clipped_points=0,
        animation=dict(points=1000000,selection='Last consecutive million states, unthinned',frames=72,
            frame_ms=100,frame_point_counts=counts,resolution=[width//2,height//2],units_per_pixel=unit*2,
            density_vmax=float(np.log1p(180))),late_return_residuals_1_to512=returns,
        scope='Rigid cameras of actual normalized (x,m,v) states. Equal spatial units on image axes; rectangular field of view discards no points. Shape label is descriptive, not a dynamical classification.')
    (out/'analysis.json').write_text(json.dumps(metadata,indent=2)+'\n')
    for name in ['render_yogi_arc.py','yogi_map.py','render_twist.py']:
        (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print('Rendered G0670; minimum late return',min(returns),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','preview','output']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
