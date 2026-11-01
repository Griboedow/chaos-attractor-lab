"""Dense actual-state portraits and rigid rotations of two RK45 candidates."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from gallery_rk45_oscillator import trajectory
from render_twist import camera,colorize
from render_yogi_arc import histogram

def layout(images,rows,scale):
    canvas=Image.new('RGB',(int(2800*scale),int(1320*scale)),'#071320');draw=ImageDraw.Draw(canvas)
    def label(x,y,text,size,color='#e6edf5'):
        draw.text((round(x*scale),round(y*scale)),text,font=ImageFont.truetype(findfont('DejaVu Sans'),round(size*scale)),fill=color)
    label(65,30,'RK45 / GEOMETRY INSIDE THE SOLVER',60)
    label(65,115,'A linear spring decays. Its adaptive numerical solver can keep moving near tolerance.',32,'#a7bbcc')
    for j,(im,r) in enumerate(zip(images,rows)):
        left=50+1400*j;canvas.paste(im,(round(left*scale),round(290*scale)))
        label(left+20,192,r['title'],40)
        label(left+20,250,f"damping={r['params'][0]:.9f}   rtol={r['params'][1]:.6g}",28,'#a7bbcc')
        label(left+20,1155,'Standard acceptance, rejection and step-size rules.',27,'#f5c181')
    label(65,1210,'Actual (position / atol, velocity / atol, proposed step) / rigid cameras / visit-density color',28,'#a7bbcc')
    label(65,1260,'Numerical-controller candidates. The physical ODE is stable; historical uniqueness remains unresolved.',26,'#a7bbcc')
    return canvas

def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    preview=json.loads(args.preview.read_text())['render'];rows=[];portraits=[];samples=[];histograms=[]
    for ident,title in [(335,'TWIN ROSETTES / D335'),(448,'SPLIT FANS / D448')]:
        p=d['params'][ident];start=d['ends'][ident]
        z=trajectory(p,start,100000,8000000);assert np.all(np.isfinite(z))
        path=out/f'case-{ident}-full.npy';np.save(path,z)
        basis=np.array(next(r for r in preview if r['id']==ident)['basis'])
        assert np.linalg.det(basis)>0 and np.max(np.abs(basis.T@basis-np.eye(3)))<1e-12
        q=z@basis;shift=(q.min(0)+q.max(0))/2;q-=shift;center=basis@shift
        radius=float(np.linalg.norm(q,axis=1).max());ext=np.max(np.abs(q),axis=0);elevation=.4;azimuth=.6
        vertical=float(np.sin(elevation)*np.hypot(ext[0],ext[2])+np.cos(elevation)*ext[1]);unit=1.045*max(2*radius/1300,2*vertical/850)
        h=histogram(q,camera(azimuth,elevation),1300,850,unit);assert int(h.sum())==len(z)
        im,vmax=colorize(h);portraits.append(im);histograms.append(h)
        sample=z[-1000000:].copy();np.savez_compressed(out/f'case-{ident}-animation.npz',states=sample,params=p);samples.append((sample-center)@basis)
        last=z[-2048:];returns=[float(np.max(np.linalg.norm(last[k:]-last[:-k],axis=1))) for k in range(1,513)]
        row=dict(id=ident,title=title,params=p.tolist(),start=start.tolist(),warm=100000,observed=len(z),source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                 trajectory_file=path.name,trajectory_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),basis=basis.tolist(),center=center.tolist(),basis_source=str(args.preview),basis_source_sha256=hashlib.sha256(args.preview.read_bytes()).hexdigest(),
                 units_per_pixel=unit,azimuth=azimuth,elevation=elevation,resolution=[1300,850],log_density_vmax=vmax,smoothing=False,clipped=0,late_return_residuals_1_to512=returns,
                 animation=dict(points=len(sample),selection='Last consecutive million states, unthinned',frames=72,frame_ms=100,resolution=[650,425],log_density_vmax=float(np.log1p(180))))
        rows.append(row);print('Dense portrait',ident,flush=True)
    layout(portraits,rows,1.).save(out/'rk45-oscillator-geometry.png');np.savez_compressed(out/'portrait-histograms.npz',histograms=histograms,ids=[335,448])
    frames=[];counts=[]
    for frame,offset in enumerate(np.linspace(0,2*np.pi,72,endpoint=False)):
        pictures=[];pointcounts=[]
        for sample,r in zip(samples,rows):
            h=histogram(sample,camera(r['azimuth']+offset,r['elevation']),650,425,2*r['units_per_pixel']);assert int(h.sum())==len(sample)
            im,_=colorize(h,vmax=np.log1p(180));pictures.append(im);pointcounts.append(int(h.sum()))
        frames.append(layout(pictures,rows,.5));counts.append(pointcounts)
        if frame%18==0:frames[-1].save(out/f'rotation-frame-{frame:02}.png');print('Frame',frame,flush=True)
    frames[0].save(out/'rk45-oscillator-rotation.gif',save_all=True,append_images=frames[1:],duration=100,loop=0,optimize=False)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,frame_point_counts=counts,scope='Unchanged actual state coordinates, proper orthogonal cameras, equal units. Names describe appearance only.'),indent=2),encoding='utf-8')
    for n in ['render_rk45_oscillator.py','gallery_rk45_oscillator.py','search_rk45_oscillator.py','render_twist.py','render_yogi_arc.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--preview',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
