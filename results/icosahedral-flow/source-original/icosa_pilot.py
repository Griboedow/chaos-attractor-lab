"""A further polyhedral screen, not yet a validated discovery."""
from pathlib import Path
import json,time
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import qmc
from tetra_vortex import rhs,survey

root=Path(__file__).parent;out=root/'icosa-v1';out.mkdir(exist_ok=False)
group=Rotation.create_group('I').as_matrix();u=qmc.Sobol(9,scramble=True,seed=2026110143).random_base2(8)
def unit(phi,cost):return np.column_stack([np.sqrt(1-cost**2)*np.cos(phi),np.sqrt(1-cost**2)*np.sin(phi),cost])
seed_c=unit(2*np.pi*u[:,0],2*u[:,1]-1);seed_axis=unit(2*np.pi*u[:,2],2*u[:,3]-1)
centers=np.einsum('gij,nj->ngi',group,seed_c);axes=np.einsum('gij,nj->ngi',group,seed_axis)
sigma=.12+.43*u[:,4];gamma=10**(-3+2.5*u[:,5]);starts=(2*u[:,6:9]-1)*1.2
directions=unit(np.arange(48)*np.pi*(3-np.sqrt(5)),1-2*(np.arange(48)+.5)/48)
probes=np.vstack([r*directions for r in [.75,1.,1.25]])
amplitude=np.zeros(256);error=0.
for i in range(256):
    f=np.array([rhs(z,np.zeros(3),centers[i],axes[i],sigma[i],0.)[0] for z in probes]);scale=float(np.sqrt((f*f).sum(1).mean()))
    assert scale>1e-8
    amplitude[i]=1/scale;axes[i]*=amplitude[i]
    if i<8:
        z=starts[i];base=rhs(z,np.zeros(3),centers[i],axes[i],sigma[i],gamma[i])[0]
        for r in group:error=max(error,float(np.max(np.abs(rhs(r@z,np.zeros(3),centers[i],axes[i],sigma[i],gamma[i])[0]-r@base))))
assert error<1e-9
np.savez_compressed(out/'input.npz',group=group,centers=centers,axes=axes,seed_center=seed_c,seed_axis=seed_axis,sigma=sigma,gamma=gamma,amplitude=amplitude,starts=starts,normalization_probes=probes)
(out/'design.json').write_text(json.dumps({'seed':2026110143,'rows':256,'normalization':'Common vortex strength chosen to make undamped field RMS speed one at 144 fixed probes. This changes the time unit relative to gamma, not the shape of a cloud after rendering.','max_equivariance_error':error,'warm_time':100,'observed_time':400,'h':.025,'keep':4000,'status':'exploratory only'},indent=2))
begin=time.time();w,z,end,speed,radius=survey(centers,axes,sigma,gamma,starts,4000,16000,.025,4000)
np.savez_compressed(out/'data.npz',windows=w,tail=z,end=end,speed=speed,radius=radius)
ids=np.where((w.min(1)>.01)&(speed>1e-4))[0]
(out/'summary.json').write_text(json.dumps({'seconds':time.time()-begin,'positive_ids':ids.tolist(),'min_amplitude':float(amplitude.min()),'max_amplitude':float(amplitude.max()),'max_equivariance_error':error},indent=2));print((out/'summary.json').read_text(),flush=True)
