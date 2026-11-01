"""Rigid-camera portraits and proper/improper symmetry sampling controls."""
from pathlib import Path
import sys
from icosa_paths import ROOT,RESEARCH,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json,sys
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(RESEARCH))
from render_twist import colorize

root=OUTPUT;d=np.load(root/'icosa-v1/input.npz');out=root/'icosa-render';out.mkdir(exist_ok=False)
phi=(1+np.sqrt(5))/2
normals=[np.array(x,dtype=float) for x in [[1,0,0],[1,1,1],[0,1,phi],[1,.37,.19],[1,.7,.9],[.3,1,.5]]]
cameras=[]
for n in normals:
    n/=np.linalg.norm(n);a=np.cross(n,[0,0,1.]);a/=np.linalg.norm(a);b=np.cross(n,a);cameras.append(np.array([a,b,n]))
fig,axs=plt.subplots(3,6,figsize=(18,9),facecolor='#080e19');rows=[]
for row,ident in enumerate([183,43,177]):
    z=np.load(root/f'icosa-long/{ident}-dense.npy');rms=np.sqrt((z*z).sum(1).mean());middle=len(z)//2
    reference=z[:middle:16];query=z[middle::80];tree=cKDTree(reference)
    identity=float(np.sqrt(np.mean(tree.query(query)[0]**2))/rms)
    proper=[float(np.sqrt(np.mean(tree.query(query@g.T)[0]**2))/rms) for g in d['group']]
    improper=[float(np.sqrt(np.mean(tree.query(-query@g.T)[0]**2))/rms) for g in d['group']]
    # Generic orbit sectors, using the original non-axis seed as the anchor.
    anchors=d['group']@d['seed_center'][ident];sector=np.argmax(z[::16]@anchors.T,axis=1);counts=np.bincount(sector,minlength=60)
    info={'id':ident,'observed_states':len(z),'identity_error':identity,'proper_errors':proper,'improper_errors':improper,'max_proper_over_identity':max(proper)/identity,'min_improper_over_identity':min(improper)/identity,'mean_over_rms':float(np.linalg.norm(z.mean(0))/rms),'sector_counts':counts.tolist(),'reference_count':len(reference),'query_count':len(query),'limits':'Finite disjoint-half nearest-neighbor distances; no exact invariant-measure symmetry or chirality proof.'};rows.append(info)
    for j,cam in enumerate(cameras):
        xy=z@cam[:2].T;half=np.abs(xy).max()*1.035
        hist=np.histogram2d(*xy.T,bins=1600,range=[[-half,half]]*2)[0].T[::-1]
        pic,vmax=colorize(hist)
        pic.save(out/f'{ident}-view-{j}.png')
        np.savez_compressed(out/f'{ident}-view-{j}.npz',histogram=hist,camera=cam,half_width=half,vmax=vmax)
        axs[row,j].imshow(pic);axs[row,j].set_axis_off();axs[row,j].set_title(f'{ident}, view {j}',color='white')
    print(json.dumps({k:v for k,v in info.items() if k not in ['proper_errors','improper_errors','sector_counts']}),flush=True)
fig.tight_layout();fig.savefig(out/'gallery.png',dpi=100);plt.close(fig)
(out/'symmetry.json').write_text(json.dumps(rows,indent=2))
