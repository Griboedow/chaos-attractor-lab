"""Resolve the regular outcome at an OCR-transcribed 1988 figure parameter pair."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from numba import set_num_threads
from twist_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    u=qmc.Sobol(3,scramble=True,seed=2026110121).random_base2(5)
    r=u[:,0]**(1/3);c=2*u[:,1]-1;a=2*np.pi*u[:,2]
    starts=np.r_[[[.1,.2,.3]],np.c_[r*np.sqrt(1-c*c)*np.cos(a),r*np.sqrt(1-c*c)*np.sin(a),r*c]]
    p=np.array([.5,0.,-np.pi/2,0.,0.,0.,-5.,0.,0.,0.]);parameters=np.tile(p,(len(starts),1))
    windows,tails,ends,means,radii=spectra(parameters,starts,100000,1000000,100000,2000)
    residual=np.array([[np.linalg.norm(t[k:]-t[:-k],axis=1).max() for k in range(1,65)] for t in tails])
    file=out/'histories.npz';np.savez_compressed(file,parameters=parameters,starts=starts,windows=windows,tails=tails,ends=ends,means=means,radii=radii,return_residual=residual)
    periods=[next((k+1 for k,e in enumerate(row) if e<1e-8),None) for row in residual]
    result=dict(source='1988 public OCR figure1 caption reads A=-5,b=.5; typeset figure not inspected. Eq.(3) transcribed separately.',
        seed=2026110121,starts='Original declared start plus32Sobol unit-ball-volume starts in the affine-conjugate coordinates.',
        warmup=100000,observed=1000000,window=100000,tail=2000,parameters=p.tolist(),
        periods=periods,positive_all_windows=int(np.count_nonzero(windows[:,:,0].min(1)>0)),
        spectrum_range=[windows.mean(1).min(0).tolist(),windows.mean(1).max(0).tolist()],
        sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        scope='A bounded reproduction check, not a correction of the source paper. OCR, source initial conditions, coexistence and finite-time effects are not resolved by this experiment.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in ['check_twist_prior_replay.py','twist_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
