"""Render a numerical cycle fold beside the motor chaos slice."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from continue_motor_braid_cycles import correct
from motor_rows import rhs
from render_motor_four import H
from audit_motor_permutation_geometry import full


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    data=np.load(args.branch/'cycles.npz');record=json.loads((args.branch/'analysis.json').read_text())
    fold=np.load(args.fold/'fold.npz');rho,r,_=fold['parameters']
    b=np.array([x['b'] for x in record['records']]);unstable=np.array([x['unstable_count'] for x in record['records']])
    select=np.flatnonzero(unstable==1);index=int(select[np.argmin(abs(b[select]-1.555))])
    guess=data['coordinates'][index].copy();guess[11]=155.5
    y,result=correct(guess,rho,r,fixed_b=1.555)
    np.savez_compressed(out/'unstable-cycle.npz',coordinates=y,seed=result['seed'],period=result['period'],parameters=result['parameters'],monodromy=result['monodromy'])
    stable=np.load(args.stable);cycles=[(stable['seed'],float(stable['period']),stable['parameters']),(result['seed'],result['period'],result['parameters'])]
    traces=[]
    for q,period,p in cycles:
        times=np.linspace(0,period,2049)
        sol=solve_ivp(lambda t,x:rhs(t,x,p,4),(0,period),q,method='DOP853',rtol=2e-12,atol=2e-14,max_step=.005,t_eval=times)
        assert sol.success and np.max(abs(sol.y[:,-1]-q))<1e-8
        traces.append(full(sol.y.T,4)[:,8:]@H.T)
    bg='#071320';fg='#e9f0f5';muted='#99afc0';cyan='#49d4c8';gold='#efb266'
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':fg,'axes.labelcolor':muted,'xtick.color':muted,'ytick.color':muted})
    fig=plt.figure(figsize=(15,8),facecolor=bg);gs=fig.add_gridspec(1,3,left=.065,right=.98,bottom=.24,top=.72,wspace=.34)
    for col,(values,ylabel,title) in enumerate([(data['periods'],'period / original time unit','Two cycles meet at a fold'),
        (np.array([x['largest_nonflow_abs'] for x in record['records']]),'largest non-flow multiplier magnitude','Stability changes at +1')]):
        ax=fig.add_subplot(gs[0,col],facecolor=bg)
        for group,color,label in [(0,cyan,'stable cycle'),(1,gold,'saddle cycle')]:
            mask=unstable==group;ax.plot(b[mask],values[mask],'.-',color=color,lw=1.4,ms=3,label=label)
        ax.axvline(fold['parameters'][2],color=muted,ls=':',lw=1)
        ax.plot(float(fold['parameters'][2]),float(fold['period']) if col==0 else 1,marker='D',color=fg,ms=4,label='fold estimate')
        if col==1:ax.axhline(1,color=muted,lw=.8);ax.set_ylim(.5,2.4)
        ax.set_xlim(1.5525,1.5605);ax.set_xlabel('rate asymmetry b');ax.set_ylabel(ylabel);ax.set_title(title,fontsize=12,pad=15)
        ax.ticklabel_format(axis='x',style='plain',useOffset=False);ax.tick_params(axis='x',rotation=25)
        for spine in ax.spines.values():spine.set_color('#31475a')
        if col==0:ax.legend(facecolor=bg,edgecolor='#31475a',labelcolor=fg,fontsize=9)
    ax=fig.add_subplot(gs[0,2],projection='3d',facecolor=bg)
    for trace,color in zip(traces,[cyan,gold]):ax.plot(*trace.T,color=color,lw=1.1,alpha=.8)
    cloud=np.concatenate(traces);center=(cloud.min(0)+cloud.max(0))/2;radius=np.ptp(cloud,axis=0).max()*.55
    for setter,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):setter(c-radius,c+radius)
    ax.set_box_aspect((1,1,1));ax.view_init(27,-45);ax.set_axis_off();ax.set_title('Separate stable and saddle cycles\nat the same b = 1.555',fontsize=12,pad=15)
    fig.text(.055,.93,'A CYCLE PAIR BESIDE THE CHAOTIC REGIME',fontsize=23,weight='bold')
    fig.text(.055,.87,'Numerical shooting and branch continuation in all 11 state variables.',fontsize=13,color=muted)
    fig.text(.055,.81,f"Fold estimate: b = {float(fold['parameters'][2]):.10f}   |   period = {float(fold['period']):.8f}",fontsize=12,color=fg)
    fig.text(.055,.09,'The two colored loops are separate periodic solutions, not one chaotic trajectory. Fixed physical displacement contrasts.\nA floating-point fold calculation does not establish a global route to chaos or historical novelty.',fontsize=11,color=muted)
    fig.savefig(out/'motor-cycle-fold.png',dpi=160,facecolor=bg);plt.close(fig)
    report=dict(branch_sha256=hashlib.sha256((args.branch/'cycles.npz').read_bytes()).hexdigest(),fold_sha256=hashlib.sha256((args.fold/'fold.npz').read_bytes()).hexdigest(),
        stable_sha256=hashlib.sha256(args.stable.read_bytes()).hexdigest(),unstable_b=1.555,unstable_period=result['period'],unstable_shooting_residual=result['residual'],
        unstable_multipliers=[[float(z.real),float(z.imag)] for z in np.linalg.eigvals(result['monodromy'])],
        figure_sha256=hashlib.sha256((out/'motor-cycle-fold.png').read_bytes()).hexdigest(),
        limitation='Numerical cycle branch only. No certified fold, chaos-organizing connection or novelty claim.')
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print('rendered',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--branch',type=Path,required=True);p.add_argument('--fold',type=Path,required=True);p.add_argument('--stable',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
