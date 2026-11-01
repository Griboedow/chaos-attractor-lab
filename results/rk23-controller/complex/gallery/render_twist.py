"""Render actual map iterates with rigid orthographic cameras and density color."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import colormaps
from PIL import Image, ImageDraw
from twist_map import trajectory


def camera(azimuth, elevation=.35):
    a, e = azimuth, elevation
    return np.array([[np.cos(a),0,np.sin(a)],
                     [np.sin(a)*np.sin(e),np.cos(e),-np.cos(a)*np.sin(e)]])


def density(points, projection, resolution, halfwidth):
    xy = points @ projection.T
    h, _, _ = np.histogram2d(xy[:,0],xy[:,1],bins=resolution,
                           range=[[-halfwidth,halfwidth]]*2)
    return h.T[::-1].astype(np.uint32)


def colorize(h, vmax=None):
    levels=np.log1p(h.astype(float))
    if vmax is None:
        vmax=float(np.quantile(levels[levels>0],.9995))
    intensity=np.minimum(levels/vmax,1)**.72
    rgb=colormaps['magma'](intensity,bytes=True)[:,:,:3]
    rgb[h==0]=0
    return Image.fromarray(rgb),vmax


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--directory',type=Path,default=Path('results/twist-geometry'))
    ap.add_argument('--assets',type=Path,default=Path('docs/assets'))
    ap.add_argument('--points',type=int,default=8000000)
    ap.add_argument('--resolution',type=int,default=2400)
    args=ap.parse_args()
    args.assets.mkdir(parents=True,exist_ok=True)
    candidates=json.loads((args.directory/'candidates.json').read_text())['candidates'][:2]
    metadata=[]
    for candidate,label,az in zip(candidates,['folded-sail','medusa'],[.4,0.]):
        p=np.array(candidate['parameters'])
        z=trajectory(p,np.array([.1,.2,.3]),100000,args.points)
        digest=hashlib.sha256(z.astype('<f8').tobytes()).hexdigest()
        center=z.mean(0)
        eigenvalues,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
        if np.linalg.det(basis)<0:
            basis[:,-1]*=-1
        coordinates=(z-center)@basis
        # Center the fixed spatial envelope, not its uneven probability mass.
        # This is a translation of the camera only, with no shape deformation.
        shift=(coordinates.max(0)+coordinates.min(0))/2
        coordinates-=shift
        center+=basis@shift
        halfwidth=float(np.max(np.linalg.norm(coordinates,axis=1))*1.06)
        projection=camera(az)
        histogram=density(coordinates,projection,args.resolution,halfwidth)
        still,vmax=colorize(histogram)
        filename=f'twist-{label}.png'
        still.save(args.assets/filename)
        # A fixed unthinned subset; only the camera moves in the animation.
        sample=coordinates[-800000:]
        frames=[]
        for angle in np.linspace(az,az+2*np.pi,72,endpoint=False):
            h=density(sample,camera(angle),540,halfwidth)
            frame,_=colorize(h,vmax=np.log1p(120))
            canvas=Image.new('RGB',(540,580),'black');canvas.paste(frame,(0,20))
            draw=ImageDraw.Draw(canvas)
            draw.text((15,6),f'{label.replace("-"," ").title()} / {candidate["id"]}',fill='#cbd7e5')
            draw.text((15,562),'Fixed data; rotating orthographic camera',fill='#8494a8')
            frames.append(canvas)
        frames[0].save(args.assets/f'twist-{label}.gif',save_all=True,
                       append_images=frames[1:],duration=100,loop=0,optimize=False)
        np.savez_compressed(args.directory/f'render-{candidate["id"]}.npz',histogram=histogram)
        row={'id':candidate['id'],'label':label,'image':filename,'warmup':100000,
             'observed':args.points,'start':[.1,.2,.3],'trajectory_sha256_le_float64':digest,
             'center':center.tolist(),'pca_basis_columns':basis.tolist(),
             'pca_eigenvalues_ascending':eigenvalues.tolist(),'halfwidth':halfwidth,
             'azimuth':az,'elevation':.35,'projection_after_pca':projection.tolist(),
             'resolution':args.resolution,'palette':'magma','log_density_vmax':vmax,
             'intensity_gamma':.72,'smoothing':False,'clipped_spatial_points':int(args.points-histogram.sum()),
             'animation':{'points':len(sample),'frames':72,'resolution':540,'log_density_vmax':float(np.log1p(120))}}
        metadata.append(row)
        print(json.dumps({'rendered':label,'sha256':digest,'clipped_spatial_points':row['clipped_spatial_points']}),flush=True)
    (args.directory/'render.json').write_text(json.dumps(metadata,indent=2)+'\n')


if __name__=='__main__':
    main()
