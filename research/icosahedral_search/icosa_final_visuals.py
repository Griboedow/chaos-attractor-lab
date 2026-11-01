"""Publishable orthographic images from the actual long trajectories."""
from pathlib import Path
import sys
from icosa_paths import ROOT,RESEARCH,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import sys,json,hashlib
import numpy as np
from PIL import Image,ImageDraw
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(RESEARCH))
from render_twist import density,colorize

root=OUTPUT;out=root/'icosa-final';out.mkdir(exist_ok=False)
rows=[]
for ident,view in [(177,2),(183,5)]:
    z=np.load(root/f'icosa-long/{ident}-dense.npy');cam=np.load(root/f'icosa-render/{ident}-view-{view}.npz')['camera'];half=float(np.max(np.linalg.norm(z,axis=1))*1.04)
    h=density(z,cam[:2],2400,half);pic,vmax=colorize(h);pic.save(out/f'icosahedral-{ident}.png')
    np.savez_compressed(out/f'portrait-{ident}.npz',histogram=h,camera=cam,half_width=half,vmax=vmax)
    sample=z[::4];frames=[];matrices=[]
    for a in np.linspace(0,2*np.pi,72,endpoint=False):
        rot=Rotation.from_rotvec(a*np.array([.3,.5,.812403840463596])).as_matrix();c=cam@rot;matrices.append(c)
        im,_=colorize(density(sample,c[:2],720,half),vmax=np.log1p(80))
        frame=Image.new('RGB',(720,756),'black');frame.paste(im,(0,18));draw=ImageDraw.Draw(frame)
        draw.text((14,4),f'Icosahedral flow / case {ident}',fill='#cbd7e5');draw.text((14,739),'One fixed trajectory; rotating orthographic camera',fill='#8494a8');frames.append(frame)
    frames[0].save(out/f'icosahedral-{ident}.gif',save_all=True,append_images=frames[1:],duration=90,loop=0,optimize=False)
    rows.append({'id':ident,'states':len(z),'trajectory_sha256_le_float64':hashlib.sha256(z.astype('<f8').tobytes()).hexdigest(),'camera':cam.tolist(),'half_width':half,'histogram_sum':int(h.sum()),'resolution':2400,'vmax':vmax,'intensity_gamma':.72,'palette':'magma','smoothing':False,'group_replication':False,'animation':{'states':len(sample),'subset_stride':4,'frames':72,'resolution':720,'vmax':float(np.log1p(80)),'cameras':np.array(matrices).tolist()}})
    print(json.dumps({'rendered':ident,'states':len(z),'histogram_sum':int(h.sum())}),flush=True)
(out/'render.json').write_text(json.dumps(rows,indent=2))
