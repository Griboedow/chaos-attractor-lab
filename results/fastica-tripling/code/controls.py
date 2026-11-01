"""Check the polynomial against the literal empirical FastICA expression."""
from pathlib import Path
import json,hashlib
import numpy as np
from model import coefficients,dataset,update
rng=np.random.default_rng(2026110105);rows=[]
for epsilon in [0.,.005,.025,.05,.1,.2,.3,.4]:
 X=dataset(epsilon);errors=[];jacerrors=[];cov=X.T@X/len(X)
 for k in range(128):
  w=rng.normal(size=4);w/=np.linalg.norm(w);v=rng.normal(size=4);v-=w*(v@w);v/=np.linalg.norm(v);f,J,hn=update(w,coefficients(epsilon));u=X@w;h=X.T@(u**3)/len(X)-3*np.mean(u*u)*w;ref=h/np.linalg.norm(h);errors.append(float(abs(f-ref).max()))
  A=3*(X.T*(u*u))@X/len(X)-3*np.mean(u*u)*np.eye(4)-6*np.outer(w,X.T@u/len(X));R=(np.eye(4)-np.outer(ref,ref))@A/np.linalg.norm(h);jacerrors.append(float(np.linalg.norm((J-R)@v)/max(1,np.linalg.norm(J@v))))
 row={'epsilon':epsilon,'rows':len(X),'mean_residual':float(abs(X.mean(0)).max()),'covariance_residual':float(abs(cov-np.eye(4)).max()),'maximum_step_difference':max(errors),'maximum_tangent_relative_difference':max(jacerrors)};rows.append(row)
 assert row['maximum_step_difference']<1e-9 and row['maximum_tangent_relative_difference']<1e-8,row
theta=np.arange(16384)*2*np.pi/16384;maximum=0.
for t in theta:
 w=np.array([np.cos(t),np.sin(t),0.,0.]);f,J,hn=update(w,coefficients(0.));target=np.array([np.cos(3*t),-np.sin(3*t),0.,0.]);maximum=max(maximum,float(abs(f-target).max()))
report={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'empirical_controls':rows,'circle_control_maximum_error':maximum};(Path(__file__).parent/'controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
