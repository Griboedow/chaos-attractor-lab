"""Rigid portraits and camera rotations for the existing tanh Yogi variant."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from yogi_tanh_map import trajectory
from render_twist import camera,density,colorize


def layout(images,rows,scale):
    canvas=Image.new('RGB',(int(2800*scale),int(1770*scale)),'#071320');draw=ImageDraw.Draw(canvas)
    def label(x,y,s,size,color='#e6edf5'):
        draw.text((int(x*scale),int(y*scale)),s,font=ImageFont.truetype(findfont('DejaVu Sans'),int(size*scale)),fill=color)
    label(65,30,'YOGI / THE SMOOTH MEMORY MAP',62)
    label(65,115,'An existing TensorFlow Addons option, on one scalar quadratic objective.',34,'#9eb5c9')
    for j,(im,r) in enumerate(zip(images,rows)):
        left=50+1400*j;canvas.paste(im,(int(left*scale),int(290*scale)))
        label(left+20,192,r['title'],42)
        label(left+20,250,f"b={r['parameters'][0]:.6f}  c={r['parameters'][1]:.6f}  e={r['parameters'][2]:.3g}  k={r['parameters'][3]:.5g}",26,'#9eb5c9')
        label(left+20,1602,r['initialization_caption'],26,'#f5c181')
    label(65,1660,'Actual normalized (x, m, v) states / rigid cameras / equal spatial units / visit-density color',28,'#9eb5c9')
    label(65,1710,'Finite numerical evidence. Neither the equation nor a new attractor class is claimed as novel.',27,'#9eb5c9')
    return canvas


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    preview=json.loads(args.preview.read_text())['runs'];rows=[];portraits=[];samples=[];histograms=[]
    for ident,title,az,caption in [(2850,'FOLDED RIBBONS / S2850',np.pi/4,'Autonomous cold start; the library cold start settles.'),(39,'FEATHERED CRESCENT / S0039',0.,'Bias-corrected start; reproduced in TensorFlow.')]:
        i=list(d['ids']).index(ident);p=d['parameters'][i];start=d['ends'][i];z=trajectory(p,start,10000,8000000)
        assert np.isfinite(z).all();path=out/f'case-{ident}-full.npy';np.save(path,z)
        basis=np.array(next(r for r in preview if r['id']==ident)['pca_basis_columns']);assert np.linalg.det(basis)>0 and np.max(abs(basis.T@basis-np.eye(3)))<1e-12
        q=z@basis;shift=(q.min(0)+q.max(0))/2;q-=shift;center=basis@shift;half=float(np.linalg.norm(q,axis=1).max()*1.055)
        h=density(q,camera(az,.25),1300,half);assert h.sum()==len(z);im,vmax=colorize(h);portraits.append(im);histograms.append(h)
        sample=z[-1000000:].copy();np.savez_compressed(out/f'case-{ident}-animation.npz',states=sample,parameters=p);samples.append((sample-center)@basis)
        last=z[-2048:];returns=[float(np.max(np.linalg.norm(last[lag:]-last[:-lag],axis=1))) for lag in range(1,513)]
        rows.append(dict(id=ident,title=title,initialization_caption=caption,parameters=p.tolist(),start=start.tolist(),source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
            warmup=10000,observed=len(z),trajectory_file=path.name,trajectory_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),basis=basis.tolist(),center=center.tolist(),
            basis_source=str(args.preview),basis_source_sha256=hashlib.sha256(args.preview.read_bytes()).hexdigest(),halfwidth=half,azimuth=az,elevation=.25,
            projection=camera(az,.25).tolist(),resolution=1300,log_density_vmax=vmax,smoothing=False,clipped=0,late_return_residuals_1_to512=returns,
            animation=dict(points=1000000,selection='Last consecutive million states, unthinned',frames=72,frame_ms=100,resolution=650,log_density_vmax=float(np.log1p(180)))))
        print('Dense portrait',ident,flush=True)
    layout(portraits,rows,1.).save(out/'yogi-tanh-geometry.png');np.savez_compressed(out/'portrait-histograms.npz',histograms=histograms,ids=[2850,39])
    frames=[];counts=[]
    for frame,offset in enumerate(np.linspace(0,2*np.pi,72,endpoint=False)):
        pictures=[];n=[]
        for sample,r in zip(samples,rows):
            h=density(sample,camera(r['azimuth']+offset,.25),650,r['halfwidth']);assert h.sum()==len(sample)
            im,_=colorize(h,vmax=np.log1p(180));pictures.append(im);n.append(int(h.sum()))
        frames.append(layout(pictures,rows,.5));counts.append(n)
        if frame%18==0:frames[-1].save(out/f'rotation-frame-{frame:02}.png');print('Frame',frame,flush=True)
    frames[0].save(out/'yogi-tanh-rotation.gif',save_all=True,append_images=frames[1:],duration=100,loop=0,optimize=False)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,frame_point_counts=counts,scope='Actual states and rigid cameras, no geometric deformation. Labels are visual descriptions, not novelty or dynamical classifications.'),indent=2)+'\n')
    for n in ['render_yogi_tanh.py','yogi_tanh_map.py','render_twist.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['source','preview','output']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
