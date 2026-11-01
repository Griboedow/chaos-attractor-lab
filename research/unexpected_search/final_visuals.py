from paths import OUTPUT, ROOT, RESEARCH
"""Publishable views, histograms and provenance, all from numerical states."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from PIL import Image,ImageDraw
sys.path.insert(0,str(RESEARCH))
from render_twist import density,colorize
from anderson_extended import survey

root=OUTPUT;out=root/'visuals';out.mkdir(exist_ok=False);metadata=[]
def camera(a,e):return np.array([[np.cos(a),np.sin(a),0],[-np.sin(a)*np.sin(e),np.cos(a)*np.sin(e),np.cos(e)]])
for ident,a,e,label in [(421,0,0,'tetrahedral-flow'),(421,.4,.6,'tetrahedral-flow-oblique'),(565,.4,.6,'tetrahedral-vortex-565')]:
    z=np.load(root/f'tetra-render/dense-{ident}.npy');cam=camera(a,e);half=float(np.abs(z@cam.T).max()*1.04)
    hist=density(z,cam,2400,half);pic,vmax=colorize(hist);pic.save(out/f'{label}.png');np.savez_compressed(out/f'{label}-histogram.npz',histogram=hist)
    metadata.append({'file':label+'.png','id':ident,'points':len(z),'camera':cam.tolist(),'half_width':half,'resolution':2400,'color_log_vmax':vmax,'sha256_points_float64_le':hashlib.sha256(z.astype('<f8').tobytes()).hexdigest()})
z=np.load(root/'tetra-render/dense-421.npy')[::8];half=float(np.linalg.norm(z,axis=1).max()*1.05)
hists=[density(z,camera(2*np.pi*i/72,.35),720,half) for i in range(72)]
vmax=max(float(np.quantile(np.log1p(h[h>0]),.9995)) for h in hists)
frames=[colorize(h,vmax)[0].quantize(256) for h in hists]
frames[0].save(out/'tetrahedral-flow.gif',save_all=True,append_images=frames[1:],duration=84,loop=0,optimize=False,disposal=2)
metadata.append({'file':'tetrahedral-flow.gif','id':421,'points_per_frame':len(z),'source_stride':8,'frames':72,'camera_elevation':.35,'camera_azimuth':'2*pi*frame/72','half_width':half,'resolution':720,'color_log_vmax':vmax,'meaning':'Rigid camera rotation of one fixed cloud. No trajectory evolution or symmetry replication.'})
d=np.load(root/'anderson-v2/data.npz');ids=np.array([54,449]);w,clouds,end,status,steps,bounds=survey(d['parameters'][ids],d['starts'][ids],100000,3200000,3200000)
np.savez_compressed(out/'anderson-render-dynamics.npz',ids=ids,windows=w,end=end,status=status,steps=steps,max_abs=bounds)
combined=Image.new('RGB',(2400,1200))
for lane,ident in enumerate(ids):
    z=clouds[lane];center=z.mean(0);_,basis=np.linalg.eigh(np.cov(z.T));basis=basis[:,::-1]
    if np.linalg.det(basis)<0:basis[:,-1]*=-1
    cam=camera(.4,.6)@basis.T;points=z-center;half=float(np.abs(points@cam.T).max()*1.04)
    hist=density(points,cam,1200,half);pic,vmax=colorize(hist);combined.paste(pic,(lane*1200,0));pic.save(out/f'anderson-{ident}.png')
    np.savez_compressed(out/f'anderson-{ident}-histogram.npz',histogram=hist);np.savez_compressed(out/f'anderson-{ident}-sample.npz',points=z[::32])
    metadata.append({'file':f'anderson-{ident}.png','id':int(ident),'points':len(z),'center':center.tolist(),'camera':cam.tolist(),'half_width':half,'resolution':1200,'color_log_vmax':vmax,'sha256_points_float64_le':hashlib.sha256(z.astype('<f8').tobytes()).hexdigest(),'windows':w[lane].tolist(),'status':int(status[lane])})
combined.save(out/'anderson-memory.png')
(out/'provenance.json').write_text(json.dumps(metadata,indent=2));print(json.dumps({'visuals':len(metadata),'anderson_windows':w.tolist(),'anderson_status':status.tolist()}),flush=True)
