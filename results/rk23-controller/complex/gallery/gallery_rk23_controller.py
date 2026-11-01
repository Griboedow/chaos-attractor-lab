"""Actual complex-state portraits and amplitude/step quotient diagnostics."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from follow_rk23_controller import cloud
from render_twist import camera,density,colorize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    ids=[396,11,719,456,860,1012];rows=[]
    for page in range(2):
        fig,axes=plt.subplots(3,2,figsize=(12,15),facecolor='#071320')
        for row,ident in enumerate(ids[page*3:page*3+3]):
            i=list(d['ids']).index(ident);p=d['params'][i]
            z,logs=cloud(p,d['ends'][i],500000)
            assert np.isfinite(z).all()
            np.savez_compressed(args.out/f'case-{ident}.npz',states=z,log_quotient=logs,params=p,start=d['ends'][i])
            center=z.mean(0);ev,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
            if np.linalg.det(basis)<0:basis[:,-1]*=-1
            q=(z-center)@basis;shift=(q.min(0)+q.max(0))/2;q-=shift;center+=basis@shift
            half=float(np.linalg.norm(q,axis=1).max()*1.03);proj=camera(.35,.35)
            hist=density(q,proj,700,half);assert int(hist.sum())==len(z)
            im,vmax=colorize(hist);ax=axes[row,0];ax.set_axis_off();ax.imshow(im,origin='upper')
            title=f'R{ident}: angle={p[0]:.6f}, rtol={p[1]:.3g}'
            ax.set_title(title+'\nActual (Re y / atol, Im y / atol, proposed h)',fontsize=10,color='#e8f0f5')
            ax=axes[row,1];ax.set_facecolor('#071320');ax.scatter(np.exp(logs[:,0]),np.exp(logs[:,1]),s=.01,c='#efc379',rasterized=True,alpha=.2)
            ax.set_xlabel('Proposed step h',color='#b5c8d8');ax.set_ylabel('Amplitude |y| / atol',color='#b5c8d8');ax.tick_params(colors='#b5c8d8')
            for spine in ax.spines.values():spine.set_color('#3a5267')
            ax.set_title(f'Quotient growth rates {d["rates"][i,:,0].mean():.5f}, {d["rates"][i,:,1].mean():.5f}\nPhase direction is neutral by exact rotational symmetry',fontsize=10,color='#e8f0f5')
            rows.append(dict(id=ident,center=center.tolist(),basis=basis.tolist(),projection=proj.tolist(),halfwidth=half,vmax=vmax,points=len(z),clipped=0))
        fig.text(.04,.965,f'RK23 / STABLE LINEAR ODE, NONLINEAR STEP CONTROLLER / {page+1}',fontsize=17,color='#e8f0f5',weight='bold')
        fig.text(.04,.025,'500,000 consecutive states. Rigid left-hand camera and equal spatial units; no log-coordinate deformation.\nThe complex ODE itself decays to zero. This experiment concerns the numerical algorithm near its tolerance.\nTranscribed SciPy rule; upstream and arbitrary-precision checks still required. No novelty claim.',fontsize=10,color='#b5c8d8')
        fig.subplots_adjust(top=.93,bottom=.09,left=.04,right=.97,hspace=.27,wspace=.22)
        fig.savefig(args.out/f'gallery-{page+1:02}.png',dpi=140,facecolor='#071320');plt.close(fig)
    report=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids,
                selection='R396 largest growth; R11 strong near-absolute tolerance; R719 second angle window; R456/R860 larger relative tolerance; R1012 smaller relative tolerance.',render=rows)
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for n in ['gallery_rk23_controller.py','follow_rk23_controller.py','search_rk23_controller.py','render_twist.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    print('Rendered six cases',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
