"""Two-dimensional numerical covering proposal, retaining all grid points."""
from pathlib import Path
from return_paths import INPUT,STRUCTURE,CYCLES,OUTPUT
import json,time
import numpy as np
import icosa_poincare as flow
from icosa_fast_return import rhs

base=OUTPUT;out=base/'icosa-return-grid';out.mkdir(exist_ok=False);flow.field=rhs
strips=np.array([[-.000008,.000002],[.002105,.00212]]);V=.0002;count=17;values=np.zeros((2,count,count,2));times=np.zeros((2,count,count));jac=np.zeros((2,count,count,2,2));speeds=np.zeros_like(times);starts=np.zeros_like(values);begin=time.time()
for k,(lo,hi) in enumerate(strips):
    for i,a in enumerate(np.linspace(lo,hi,count)):
        for j,b in enumerate(np.linspace(-V,V,count)):
            uv=np.array([a,b]);q,t,J,speed=flow.evaluate(uv);starts[k,i,j]=uv;values[k,i,j]=q;times[k,i,j]=t;jac[k,i,j]=J;speeds[k,i,j]=speed
        if i%4==0:print(json.dumps({'strip':k,'row':i+1,'of':count,'seconds':time.time()-begin}),flush=True)
rows=[]
for k in range(2):
    J=jac[k];expand=np.abs(J[:,:,0,0])-.1*np.abs(J[:,:,0,1]);slope=(np.abs(J[:,:,1,0])+.1*np.abs(J[:,:,1,1]))/expand
    rows.append({'strip':k,'u_interval':strips[k].tolist(),'v_interval':[-V,V],'sample_count':count*count,'minimum_left_edge_u_output':float(values[k,0,:,0].min()),'maximum_right_edge_u_output':float(values[k,-1,:,0].max()),'left_margin_beyond_both_targets':float(values[k,0,:,0].min()-strips[:,1].max()),'right_margin_beyond_both_targets':float(strips[:,0].min()-values[k,-1,:,0].max()),'maximum_abs_v_output':float(abs(values[k,:,:,1]).max()),'stable_margin':float(V-abs(values[k,:,:,1]).max()),'return_time_range':[float(times[k].min()),float(times[k].max())],'minimum_crossing_speed':float(speeds[k].min()),'unstable_cone_slope':.1,'sampled_minimum_unstable_expansion':float(expand.min()),'sampled_maximum_image_cone_slope':float(slope.max())})
np.savez_compressed(out/'data.npz',inputs=starts,outputs=values,times=times,derivatives=jac,crossing_speeds=speeds,strips=strips,stable_half_width=V)
report={'rows':rows,'seconds':time.time()-begin,'all_sampled_margins_positive':all(min(r['left_margin_beyond_both_targets'],r['right_margin_beyond_both_targets'],r['stable_margin'])>0 for r in rows),'limits':'Only 17 by 17 point samples on each strip. No bounds between samples, validated flow tube, all-point crossing, compact invariant set, horseshoe, or entropy lower bound is certified.'};(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
