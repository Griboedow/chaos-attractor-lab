"""Numerical recurrence, cycle multipliers and fixed-view motor-slice figures."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from check_motor_asymmetry import augmented

H=np.array([[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]],float)/2


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    scout=np.load(args.scout);runs=json.loads((args.source/'analysis.json').read_text())['runs']
    result=[];plots=[]
    for r in runs:
        f=args.source/r['file'];d=np.load(f);p=d['parameters'];s=d['section_states'][-800:];scale=1+np.max(abs(s),axis=0)
        errors=[float(np.max(abs(s[lag:]-s[:-lag])/scale)) for lag in range(1,65)]
        close=[i+1 for i,e in enumerate(errors) if e<1e-7];lag=min(close) if close else None
        record=dict(id=r['id'],b=float(p[2]),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),spectrum=r['spectrum'],
            first_numerical_period_lag=lag,maximum_scaled_recurrence_errors=errors,section_count=len(d['section_times']))
        if lag:
            t=d['section_times'];period=float(np.mean(t[-100:][lag:]-t[-100:][:-lag]));seed=d['section_states'][-1]
            sol=solve_ivp(lambda t,z:augmented(z,p,4),(0,period),np.r_[seed,np.eye(11).ravel(),0.],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.01)
            assert sol.success
            monodromy=sol.y[11:-1,-1].reshape(11,11);eig=np.linalg.eigvals(monodromy)
            eig=eig[np.argsort(abs(eig))[::-1]];neutral=int(np.argmin(abs(eig-1)))
            record.update(period=period,cycle_return_error=float(np.max(abs(sol.y[:11,-1]-seed))),
                monodromy_multipliers=[[float(z.real),float(z.imag)] for z in eig],
                largest_nonflow_multiplier_abs=float(np.max(abs(np.delete(eig,neutral)))))
            np.savez_compressed(out/f'cycle-{r["id"]:04}.npz',seed=seed,parameters=p,period=period,monodromy=monodromy,multipliers=eig)
        # Relabel only existing rows, with a deterministic rule, before the fixed linear contrast.
        x=np.column_stack([d['states'][:,8:],-d['states'][:,8:].sum(axis=1)])
        distinguished=int(np.argmin(x.mean(axis=0)));other=[i for i in range(4) if i!=distinguished]
        order=[other[0],other[1],distinguished,other[2]]
        planar=np.column_stack([x[:,other[0]]-x[:,other[1]],(x[:,other[0]]+x[:,other[1]]-2*x[:,other[2]])/np.sqrt(3)])
        circulation=float(np.sum(planar[:-1,0]*planar[1:,1]-planar[:-1,1]*planar[1:,0]))
        if circulation<0:order[0],order[1]=order[1],order[0]
        mask=d['time']>=d['time'][-1]-1200;points=x[mask][:,order]@H.T
        record.update(row_order=order,distinguished_row=distinguished,net_circulation_sign=int(np.sign(circulation)),plotted_points=len(points))
        np.savez_compressed(out/f'plot-{r["id"]:04}.npz',contrasts=points,parameters=p,row_order=order)
        plots.append(points);result.append(record)
    bg='#071320';fg='#dbe6ef';muted='#94a8ba';gold='#efb266';cyan='#49d4c8'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,'ytick.color':muted})
    fig=plt.figure(figsize=(18,10),facecolor=bg);gs=fig.add_gridspec(2,4,left=.04,right=.98,bottom=.09,top=.78,wspace=.12,hspace=.35)
    ax=fig.add_subplot(gs[0,0],facecolor=bg)
    for start in range(2):
        ids=np.arange(start,len(scout['parameters']),2);ax.scatter(scout['parameters'][ids,2],scout['lle'][ids],s=8,alpha=.6,color=[gold,cyan][start],label=f'fresh start {start+1}')
    ax.axhline(0,color=muted,lw=.6);ax.set_xlim(1.45,1.61);ax.set_ylim(-.003,.3)
    ax.set_xlabel('rate asymmetry b');ax.set_ylabel('leading finite-time rate');ax.legend(fontsize=8,facecolor=bg,labelcolor=fg,edgecolor='#273b4d')
    ax.set_title('The local parameter slice',fontsize=11,pad=15)
    for spine in ax.spines.values():spine.set_color('#273b4d')
    allp=np.concatenate(plots);center=(allp.max(0)+allp.min(0))/2;radius=np.ptp(allp,axis=0).max()*.54
    for index,(r,points) in enumerate(zip(result,plots),1):
        ax=fig.add_subplot(gs[index//4,index%4],projection='3d',facecolor=bg)
        color=cyan if r['first_numerical_period_lag'] else gold
        ax.plot(*points.T,color=color,lw=.45 if r['first_numerical_period_lag'] else .35,alpha=.8 if r['first_numerical_period_lag'] else .68)
        for setter,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):setter(c-radius,c+radius)
        ax.view_init(27,-45);ax.set_box_aspect((1,1,1));ax.set_axis_off()
        word='numerical cycle' if r['first_numerical_period_lag'] else 'positive leading exponent'
        ax.set_title(f'b = {r["b"]:g}  |  {word}\nleading rate {r["spectrum"][0]:.5f}',fontsize=10,pad=5)
    fig.text(.045,.95,'HOW THE MOTOR BRAID CHANGES',fontsize=23)
    fig.text(.045,.9,'Fixed mechanical parameters; only the rate asymmetry changes. Seven independent full-spectrum checks, 6000 time units each.',fontsize=12,color=muted)
    fig.text(.045,.855,'The three green traces close numerically in all 11 coordinates. The four gold traces retain a positive leading exponent.',fontsize=11,color=muted)
    fig.text(.045,.04,'Last 1200 time units; fixed orthonormal displacement contrasts, common bounds and camera. Rows relabeled by a declared rule.\nFinite numerical evidence; this figure does not determine the bifurcation mechanism or establish a new attractor type.',fontsize=10,color=muted)
    path=out/'motor-four-b-slice.png';fig.savefig(path,dpi=150,facecolor=bg);plt.close(fig)
    report=dict(scout_sha256=hashlib.sha256(args.scout.read_bytes()).hexdigest(),cases=result,contrast_matrix=H.tolist(),
        recurrence='Last800 solver-root x0=0 upward events, all11 coordinates, scale1+maxabs, lags1..64, maximum<1e-7.',
        row_relabeling='Smallest mean displacement moved to slot2; remaining rows in sorted order, swap first two if their net planar circulation is negative. No coordinate fitting or deformation.',
        figure_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),limitations='Finite numerical cycles and multipliers; no interval validation, bifurcation identification or topological novelty.')
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print([(r['b'],r.get('largest_nonflow_multiplier_abs')) for r in result],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--scout',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
