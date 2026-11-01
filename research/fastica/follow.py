"""Full tangent spectra and literal weight trajectories on the unit 3-sphere."""
from pathlib import Path
from paths import output_root
import hashlib,json,math,time
import numpy as np
from numba import njit,prange

@njit(cache=True)
def map_jac(w,X,kind,alpha):
 n,d=X.shape;h=np.zeros(d);A=np.zeros((d,d));avgp=0.;b=np.zeros(d)
 for i in range(n):
  u=0.
  for j in range(d):u+=X[i,j]*w[j]
  if kind==0:g=math.tanh(alpha*u);gp=alpha*(1-g*g);gpp=-2*alpha*gp*g
  elif kind==1:e=math.exp(-u*u/2);g=u*e;gp=(1-u*u)*e;gpp=(u*u*u-3*u)*e
  else:g=u*u*u;gp=3*u*u;gpp=6*u
  avgp+=gp/n
  for j in range(d):
   h[j]+=X[i,j]*g/n;b[j]+=X[i,j]*gpp/n
   for k in range(d):A[j,k]+=X[i,j]*X[i,k]*gp/n
 for j in range(d):
  h[j]-=avgp*w[j]
  for k in range(d):A[j,k]-=w[j]*b[k]
  A[j,j]-=avgp
 norm=np.sqrt(np.dot(h,h));f=h/norm;J=np.empty_like(A)
 for k in range(d):
  c=0.
  for j in range(d):c+=f[j]*A[j,k]
  for j in range(d):J[j,k]=(A[j,k]-f[j]*c)/norm
 return f,J,norm

@njit(cache=True)
def integrate(X,w,kind,alpha,warmup,steps,window,stride):
 w=w/np.sqrt(np.dot(w,w));Q=np.eye(4)[:,:3].copy()
 # Householder reflection produces an orthonormal tangent frame at the start.
 a=w.copy();a[3]+=1. if w[3]>=0 else -1.;a/=np.sqrt(np.dot(a,a))
 for j in range(3):
  for i in range(4):Q[i,j]-=2*a[i]*a[j]
 sums=np.zeros((steps//window,3));cloud=np.empty((steps//stride,4));minimum=1e9;min_diag=1e9
 for it in range(warmup+steps):
  w,J,hn=map_jac(w,X,kind,alpha);minimum=min(minimum,hn);V=J@Q
  for j in range(3):
   V[:,j]-=w*np.dot(w,V[:,j])
   for k in range(j):V[:,j]-=Q[:,k]*np.dot(Q[:,k],V[:,j])
   g=np.sqrt(np.dot(V[:,j],V[:,j]));min_diag=min(min_diag,g);Q[:,j]=V[:,j]/g
   if it>=warmup:sums[(it-warmup)//window,j]+=math.log(g)/window
  if it>=warmup and (it-warmup)%stride==0:cloud[(it-warmup)//stride]=w
 return sums,cloud,minimum,min_diag,w

if __name__=='__main__':
 root=output_root();inp=np.load(root/'pilot-01/inputs.npz');report=json.loads((root/'pilot-01/report.json').read_text());out=root/'follow-01';out.mkdir(exist_ok=False)
 protocol={'cases':report['retained'],'warmup':5000,'observed':200000,'window':20000,'stride':2,'dimension':4,'state':'Unmodified unit-norm weight vector. No sign folding.','tangent':'Analytic normalized-update Jacobian; QR on the three-dimensional tangent space.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
 rows=[]
 for i in protocol['cases']:
  start=time.time();X=inp['X'][i,:inp['counts'][i]];w=inp['starts'][i];kind=inp['kinds'][i];alpha=inp['alphas'][i]
  rate,cloud,mn,mg,last=integrate(X,w,kind,alpha,5000,200000,20000,2)
  np.savez_compressed(out/f'case-{i}.npz',cloud=cloud,rates=rate,endpoint=last)
  row={**report['cases'][i],'full_spectrum_mean':rate.mean(axis=0).tolist(),'spectrum_windows':rate.tolist(),'minimum_update_norm':float(mn),'minimum_qr_diagonal':float(mg),'sphere_residual':float(abs(np.sum(cloud**2,axis=1)-1).max()),'wall_seconds':time.time()-start};rows.append(row)
  (out/f'case-{i}-input.json').write_text(json.dumps({'case':i,'X':X.tolist(),'start':w.tolist(),'kind':int(kind),'alpha':float(alpha)},indent=2)+'\n');print(json.dumps(row),flush=True)
 (out/'report.json').write_text(json.dumps(rows,indent=2)+'\n')
