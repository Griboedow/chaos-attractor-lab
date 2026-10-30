"""Render rotating physical projections of two archived eight-row motor regimes.

This is a read-only visualization of saved trajectories, not a new integration
or a proof calculation. The periodic illustration is the floating-point seed
of the separately certified orbit; the irregular trajectory remains a candidate.
"""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from PIL import Image, ImageChops

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'docs/assets/motor-regimes-rotation')
    parser.add_argument('--preview-only',action='store_true')
    args=parser.parse_args()
    stem=args.output.resolve()
    stem.parent.mkdir(parents=True,exist_ok=True)
    inputs={
        'periodic':ROOT/'results/motor-activity-continuation/probe-destinations/case-16/periodic/cycle.npz',
        'certificate':ROOT/'results/motor-five-cycle-proof/recheck-01.json',
        'irregular':ROOT/'results/motor-activity-continuation/cpu/run-01/block-16.npz',
        'spectrum':ROOT/'results/motor-activity-continuation/cpu/run-01/analysis.json',
    }
    cycle=np.load(inputs['periodic'])
    irregular=np.load(inputs['irregular'])
    certificate=json.loads(inputs['certificate'].read_text())
    spectrum=json.loads(inputs['spectrum'].read_text())
    assert certificate['status']=='attracting_five_cluster_cycle_certified'
    assert spectrum['status']=='complete' and spectrum['observed_tau']==20000
    assert np.all(np.asarray(spectrum['spectrum'])[:6]>0)
    assert np.all(np.isfinite(cycle['full_states'])) and np.all(np.isfinite(irregular['states']))
    center=np.array([math.ldexp(float(v['midpoint'][0]),v['midpoint'][1]) for v in certificate['initial_state_box']])
    seed_difference=float(np.max(abs(cycle['quotient_initial_state']-center)))
    closure=float(np.max(abs(cycle['full_states'][-1]-cycle['full_states'][0])))
    assert seed_difference<1e-10 and closure<1e-9
    indices=[16,17,8]  # Original zero-based X0, X1, s0; no fitted embedding.
    curves=[cycle['full_states'][:,indices],irregular['states'][:,indices]]
    combined=np.vstack(curves)
    lower,upper=combined.min(axis=0),combined.max(axis=0)
    margin=.10*(upper-lower)
    lower-=margin;upper+=margin
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'text.color':'#23354a','axes.labelcolor':'#344a60',
                         'xtick.color':'#526579','ytick.color':'#526579'})
    fig=plt.figure(figsize=(11.6,6.4),dpi=100,facecolor='#f8fafc')
    fig.text(.5,.95,'Motor-row dynamics',ha='center',fontsize=22,fontweight='bold')
    fig.text(.5,.905,'Rotating views of the same physical coordinates',ha='center',fontsize=12,color='#596d80')
    axes=[]
    for index in range(2):
        ax=fig.add_axes([.035+.50*index,.20,.43,.58],projection='3d',facecolor='#f8fafc')
        ax.set(xlim=(lower[0],upper[0]),ylim=(lower[1],upper[1]),zlim=(lower[2],upper[2]),
               xlabel=r'$X_0$',ylabel=r'$X_1$',zlabel=r'$s_0$')
        ax.set_box_aspect(upper-lower)
        ax.set_proj_type('ortho')
        ax.tick_params(labelsize=8,pad=0)
        for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
            axis.set_pane_color((.93,.95,.97,1))
            axis._axinfo['grid'].update(color=(.70,.76,.82,.40),linewidth=.6)
        axes.append(ax)
    axes[0].plot(*curves[0].T,color='#c15c2e',lw=2.4,solid_capstyle='round')
    axes[0].scatter(*curves[0][0],color='#512b23',s=28,depthshade=False,zorder=5)
    axes[1].scatter(*curves[1].T,s=.8,c='#207f94',alpha=.38,depthshade=False,edgecolors='none')
    for x,title,subtitle,color in [
        (.25,'Certified periodic attractor','8 rows · 5 synchronized clusters','#a54d28'),
        (.75,'Hyperchaos candidate','6 positive finite-time exponents','#176e81')]:
        fig.text(x,.84,title,ha='center',fontsize=15,fontweight='bold',color=color)
        fig.text(x,.805,subtitle,ha='center',fontsize=10,color='#596d80')
    fig.text(.25,.165,r'$\rho=129/55$ · $T\approx25.00168$',ha='center',fontsize=11)
    fig.text(.75,.165,r'$\rho\approx2.93409$ · spectrum over $20\,000\,\tau$',ha='center',fontsize=11)
    fig.text(.5,.092,'Camera rotation only · line: one period · points: the final 1,000 time units',
             ha='center',fontsize=10,color='#4c6276')
    fig.text(.5,.047,'23-dimensional system → 3-dimensional projection · identical axis limits · zero-based row labels',
             ha='center',fontsize=9,color='#6c7d8f')
    frames=96;fps=12
    def update(frame):
        angle=2*np.pi*frame/frames
        for ax in axes:ax.view_init(elev=23+8*np.sin(angle),azim=-55+360*frame/frames)
        return []
    update(0)
    fig.savefig(stem.with_suffix('.png'),dpi=150,facecolor=fig.get_facecolor())
    if args.preview_only:
        plt.close(fig)
        print('Preview saved')
        return
    movie=FuncAnimation(fig,update,frames=frames,interval=1000/fps)
    movie.save(stem.with_suffix('.gif'),writer=PillowWriter(fps=fps),dpi=100)
    plt.close(fig)
    gif=Image.open(stem.with_suffix('.gif'))
    assert gif.n_frames==frames and gif.info.get('loop')==0
    checkpoints=[]
    for i in (0,24,48,72):
        gif.seek(i);checkpoints.append(gif.convert('RGB').copy())
    differences=[ImageChops.difference(checkpoints[0],q).getbbox() for q in checkpoints[1:]]
    assert all(differences)
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_sha256=sha(Path(__file__)),
        inputs={p.relative_to(ROOT).as_posix():sha(p) for p in inputs.values()},
        no_integration_or_proof_changes=True,projection=dict(state_indices=indices,labels=['X0','X1','s0'],
            dimensionality=23,fitted_embedding=False,shared_axis_lower=lower.tolist(),shared_axis_upper=upper.tolist()),
        periodic=dict(samples=len(curves[0]),saved_period=float(cycle['time'][-1]),
            maximum_seed_distance_from_certificate_midpoint=seed_difference,maximum_full_state_closure_error=closure,
            qualification='Floating-point illustration of the separately certified periodic orbit; not an interval-enclosed image.'),
        candidate=dict(samples=len(curves[1]),original_time_range=[float(irregular['time'][0]),float(irregular['time'][-1])],
            spectrum_observed_tau=spectrum['observed_tau'],first_six_finite_time_exponents=spectrum['spectrum'][:6],
            qualification='The final saved 1000-tau block is displayed; finite-time hyperchaos and asymptotic attraction remain distinct claims.'),
        animation=dict(frames=frames,fps_requested=fps,pixel_dimensions=list(gif.size),loop=gif.info['loop'],
            recorded_frame_duration_ms=gif.info.get('duration'),different_quarter_turn_frames=all(differences),
            camera_rotation_degrees=360,trajectory_animation=False),
        outputs={stem.with_suffix(suffix).name:dict(bytes=stem.with_suffix(suffix).stat().st_size,sha256=sha(stem.with_suffix(suffix)))
                 for suffix in ('.png','.gif')})
    stem.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(frames=frames,files=record['outputs'],seed_difference=seed_difference,closure=closure)))


if __name__=='__main__':main()
