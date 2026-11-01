from pathlib import Path
import sys,json,hashlib
import numpy as np
from PIL import Image,ImageDraw
import matplotlib.pyplot as plt
sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')
from gaussian_vortex import integrate
from render_twist import colorize,density
base=Path('D:/Codex/chaos-attractor-lab');out=base/'results/gaussian-vortex';assets=base/'docs/assets'
p=np.array(json.loads((out/'candidates.json').read_text())['candidates'][0]['parameters'])
hs=np.array([.00625]);w,z,end,means,radius=integrate(p[None,:],np.array([[.1,.2,.3]]),hs,2000.,30000.,1000.,4800000,1)
z=z[0]
np.savez_compressed(out/'render-run.npz',parameters=p,stepsize=.00625,warm_time=2000.,observed_time=30000.,window_time=1000.,windows=w[0],tail=z[-100000:],end=end[0],means=means[0],max_radius=radius[0])
np.savez_compressed(Path(__file__).parent/'render-full.npz',points=z)
basis=np.column_stack([np.array([1.,-1,0])/np.sqrt(2),np.array([1.,1,-2])/np.sqrt(6),np.ones(3)/np.sqrt(3)])
center=z.mean(0);coords=(z-center)@basis
shift=(coords.max(0)+coords.min(0))/2;coords-=shift;center+=basis@shift
halfwidth=float(np.linalg.norm(coords,axis=1).max()*1.06)
def camera(az,el):
 return np.array([[np.cos(az),np.sin(az),0],[-np.sin(az)*np.sin(el),np.cos(az)*np.sin(el),np.cos(el)]])
fig,axes=plt.subplots(2,3,figsize=(12,8),facecolor='black')
for ax,(az,el) in zip(axes.flat,[(0,30),(30,45),(60,60),(0,90),(60,15),(30,75)]):
 h=density(coords,camera(np.deg2rad(az),np.deg2rad(el)),650,halfwidth)
 pic,_=colorize(h);ax.imshow(pic);ax.set_axis_off();ax.set_title(f'az={az} el={el}',color='white')
fig.tight_layout();fig.savefig(Path(__file__).parent/'vortex-views.png',dpi=125);plt.close(fig)
az=np.deg2rad(30);el=np.deg2rad(45)
h=density(coords,camera(az,el),2400,halfwidth);pic,vmax=colorize(h);pic.save(assets/'gaussian-vortex-triskelion.png')
np.savez_compressed(out/'render-histogram.npz',counts=h)
frames=[]
for angle in np.linspace(az,az+2*np.pi,72,endpoint=False):
 hh=density(coords[-600000:],camera(angle,el),600,halfwidth)
 frame,_=colorize(hh,np.log1p(160));canvas=Image.new('RGB',(600,640),'black');canvas.paste(frame,(0,20));draw=ImageDraw.Draw(canvas)
 draw.text((15,6),'Gaussian vortex flow / candidate 77',fill='#cbd7e5')
 draw.text((15,622),'Fixed trajectory; rotating orthographic camera',fill='#8494a8');frames.append(canvas)
frames[0].save(assets/'gaussian-vortex-triskelion.gif',save_all=True,append_images=frames[1:],duration=100,loop=0,optimize=False)
meta={'candidate_id':77,'points':len(z),'trajectory_sha256_le_float64':hashlib.sha256(z.astype('<f8').tobytes()).hexdigest(),
      'basis_columns':basis.tolist(),'center':center.tolist(),'halfwidth':halfwidth,'azimuth':float(az),'elevation':float(el),'projection_rows':camera(az,el).tolist(),
      'resolution':2400,'clipped_points':int(len(z)-h.sum()),'palette':'magma','intensity_gamma':.72,'log_count_saturation':vmax,
      'animation_points':600000,'animation_frames':72,'animation_log_count_saturation':float(np.log1p(160)),'smoothing':False}
(out/'render.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps({'spectrum':w.mean(1).tolist(),'min_leading_window':float(w[:,:,0].min()),'points':len(z)}),flush=True)
