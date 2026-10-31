"""Reintegrate cited reference flows for a qualitative silhouette comparison."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from render_motor_four import H
from audit_motor_permutation_geometry import full
from render_motor_four_ribbon import style


@njit(cache=True)
def rhs(t,q,kind):
    x,y,z=q
    if kind==0:return np.array([-1.27*x-4*y-4*z-y*y,-1.27*y-4*z-4*x-z*z,-1.27*z-4*x-4*y-x*x])
    if kind==1:return np.array([10*(y-x),x*(180.78-z)-y,x*y-(8/3)*z])
    return np.array([y,z,x*x*x-x-y-.3783*z])


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=[dict(name='Halvorsen / Sprott',parameter='a = 1.27',url='https://sprott.physics.wisc.edu/chaos/symmetry.htm',
        interpretation='Cyclic quadratic flow; symmetry and period-doubling precedent.',warm=1000.,observation=1000.,step=.02),
        dict(name='Lorenz: a narrow chaotic regime',parameter='sigma = 10, beta = 8/3, rho = 180.78',url='https://doi.org/10.1038/s41598-017-01083-x',
        interpretation='Published partially predictable chaos example, not the familiar rho=28 butterfly.',warm=200.,observation=200.,step=.002),
        dict(name='Cubic jerk: a narrow chaotic regime',parameter='b = 0.3783',url='https://doi.org/10.1038/s41598-017-01083-x',
        interpretation='Equation15 and Figure9 of Wernecke etal2017; published braided-chaos example.',warm=1000.,observation=1000.,step=.02)]
    clouds=[]
    for k,record in enumerate(sources):
        grid=np.linspace(record['warm'],record['warm']+record['observation'],1+round(record['observation']/record['step']))
        sol=solve_ivp(rhs,(0,grid[-1]),[.1,0,0],args=(k,),method='DOP853',rtol=2e-11,atol=2e-13,max_step=record['step'],t_eval=grid)
        assert sol.success and np.max(abs(sol.y))<1000
        path=out/f'reference-{k}.npz';np.savez_compressed(path,time=grid,states=sol.y.T)
        record.update(initial=[.1,0,0],method='DOP853,rtol2e-11,atol2e-13,maxstep=sample spacing.',
            samples=len(grid),file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        clouds.append(sol.y.T)
    d=np.load(args.source);x=d['states'][-30001::2];motor=full(x,4)[:,8:]@H.T
    clouds=[motor]+clouds
    names=['Four motor groups: current candidate']+[r['name'] for r in sources]
    subtitles=['Orthonormal displacement contrasts; 11D full state']+[r['parameter'] for r in sources]
    colors=['#f4b879','#4bd4cd','#a2a8ff','#eb8daf']
    bg='#071320';fg='#e9f0f5';muted='#99afc0'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg})
    fig=plt.figure(figsize=(15,11),facecolor=bg)
    for k,(cloud,title,subtitle,color) in enumerate(zip(clouds,names,subtitles,colors)):
        ax=fig.add_subplot(2,2,k+1,projection='3d',facecolor=bg)
        ax.plot(*cloud.T,color=color,lw=.24,alpha=.5)
        style(ax,cloud,27,-45,zoom=1.00)
        ax.set_title(title+'\n'+subtitle,fontsize=12,pad=1)
    fig.suptitle('VISUAL NOVELTY NEEDS A STRONGER COMPARISON',x=.04,ha='left',y=.98,fontsize=23,weight='bold')
    fig.text(.04,.92,'Known systems already produce cyclic shapes and very thin chaotic braids.\nSeparate coordinate systems and panel scales; the same camera and equal axis units within each panel.',fontsize=12,color=muted)
    fig.text(.04,.037,'Reference equations: Sprott (1997/2004) and Wernecke, Sandor & Gros (2017). All curves are numerical integrations.\nA visual comparison neither proves equivalence nor establishes a new topology. The motor regime remains under investigation.',fontsize=11,color=muted)
    fig.subplots_adjust(left=.02,right=.98,top=.855,bottom=.11,wspace=.01,hspace=.08)
    fig.savefig(out/'known-comparison.png',dpi=150,facecolor=bg);plt.close(fig)
    (out/'render.json').write_text(json.dumps(dict(motor_source=str(args.source),motor_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        motor_display='Final600tau of source12ktau history, stride2/.04tau.',reference_runs=sources,
        coordinate_rule='Motor:fixed orthonormal H X. Reference systems:unmodified x,y,z. Per-panel bounds/equal axis units; no deformation or fitted embedding.',
        limitation='Qualitative appearance only. Reference systems are known; no new scientific result claimed for their reintegration.'),indent=2)+'\n',encoding='utf8')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print('completed three reference integrations and comparison',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
