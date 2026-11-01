"""Recompute published kicked-top examples and compare their geometry fairly."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import set_num_threads
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from twist_map import spectra,trajectory
from render_twist import camera,density,colorize


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    source=json.loads((args.reference/'candidates.json').read_text())['candidates']
    labels=['1988 / A = -5, b = 0.5','1988 / A = 35, b = 0.1','1988 / A = 15, b = 0.5']
    parameters=[np.array([b,0.,-np.pi/2,0.,0.,0.,a,0.,0.,0.]) for b,a in [(.5,-5.),(.1,35.),(.5,15.)]]
    for ident in [1709,3254,307]:
        parameters.append(np.array(next(c['parameters'] for c in source if c['id']==ident)))
        labels.append(f'Selected twist / {ident}')
    parameters=np.array(parameters);starts=np.tile([.1,.2,.3],(6,1))
    windows,tails,ends,means,radii=spectra(parameters,starts,100000,1000000,100000,20000)
    np.savez_compressed(out/'spectra.npz',parameters=parameters,starts=starts,windows=windows,tails=tails,ends=ends,means=means,radii=radii)
    bg='#071320';fg='#e8f0f5';muted='#9ab0c1';plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg})
    fig,axes=plt.subplots(2,3,figsize=(15,11),facecolor=bg);records=[]
    for i,(p,label,ax) in enumerate(zip(parameters,labels,axes.ravel())):
        z=trajectory(p,starts[i],100000,1000000);center=z.mean(0)
        ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
        if np.linalg.det(basis)<0:basis[:,-1]*=-1
        q=(z-center)@basis;shift=(q.max(0)+q.min(0))/2;q-=shift;center+=basis@shift
        halfwidth=float(np.max(np.linalg.norm(q,axis=1))*1.04);projection=camera(0. if i==4 else .4)
        h=density(q,projection,650,halfwidth);assert h.sum()==len(z)
        picture,vmax=colorize(h);picture.save(out/f'panel-{i}.png')
        np.savez_compressed(out/f'panel-{i}.npz',histogram=h)
        ax.imshow(picture);ax.set_axis_off();ax.set_title(label,color=fg,fontsize=15,pad=9)
        exponents=windows[i].mean(0)
        ax.text(.02,.03,'λ: '+', '.join(f'{v:.3f}' for v in exponents),transform=ax.transAxes,color=muted,fontsize=11)
        records.append(dict(label=label,parameters=p.tolist(),start=starts[i].tolist(),warmup=100000,observed=1000000,
            spectrum=exponents.tolist(),leading_window_range=[float(windows[i,:,0].min()),float(windows[i,:,0].max())],
            trace_identity_error=float(np.max(abs(windows[i].sum(1)-3*np.log(p[0])))),
            trajectory_sha256_le_float64=hashlib.sha256(z.astype('<f8').tobytes()).hexdigest(),center=center.tolist(),
            pca_basis_columns=basis.tolist(),projection_after_pca=projection.tolist(),halfwidth=halfwidth,histogram_resolution=650,
            log_density_vmax=vmax,clipped_points=0))
        print(label,exponents,flush=True)
    fig.text(.04,.96,'GEOMETRY WITH A KNOWN PREDECESSOR',fontsize=23,weight='bold')
    fig.text(.04,.918,'Top: published 1988 examples recomputed here. Bottom: selected points of the broader three-axis family.',fontsize=12,color=muted)
    fig.text(.04,.055,'One million actual iterates per panel; rigid PCA and orthographic projection, equal scale on both image axes.',fontsize=11,color=muted)
    fig.text(.04,.025,'Each panel fits its own spatial extent and log-density color scale. Appearance alone does not identify a dynamical class.',fontsize=10,color=muted)
    fig.subplots_adjust(left=.025,right=.975,top=.87,bottom=.09,wspace=.025,hspace=.18)
    fig.savefig(out/'twist-kicked-top-comparison.png',dpi=160,facecolor=bg);plt.close(fig)
    (out/'analysis.json').write_text(json.dumps(dict(source='Peplowski and Stefanski1988,Eq.(3),Figs.1-3; scalar b case.',
        coordinate_change='The top row uses the explicitly affine-conjugate twist-map slice. This change is an isotropic scale, orthogonal transformation and translation.',
        scope='Own recomputations of known examples, not reproductions of copyrighted figures; no topological equivalence is inferred from visual similarity.',runs=records),indent=2)+'\n')
    for name in ['compare_twist_kicked_top.py','twist_map.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
