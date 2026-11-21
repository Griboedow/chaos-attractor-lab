"""Regenerate one extension cloud from its portable record and print the actual SHA-256."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from events import reconstruct
from adaptive import full_rhs
from field import integrate,ring_tangent
from relative_cores import sector_cloud


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--manifest',required=True,type=Path)
    ap.add_argument('--id',required=True);ap.add_argument('--out',required=True,type=Path)
    a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    a.out.parent.mkdir(exist_ok=True,parents=True)
    row=next(r for r in json.loads(a.manifest.read_text())['cases'] if r['identifier']==a.id)
    p=np.array(row['parameters']);n=row['shape'][0]
    if row['recipe']=='rotated_nodes':
        z,_=reconstruct(row,p,n)
    elif row['recipe']=='rk4':
        h=row['step'];_,z=integrate(p,np.array(row['initial']),h,round(row['burn']/h),round(row['window']/h),
                                      row['blocks'],row['stride'],ring_tangent)
    elif row['recipe']=='single_relative_sector':
        z,_=sector_cloud(np.array(row['initial']),row['sector_period'],p,row['rotation_index'],n)
    else:
        T=row['period'];sol=solve_ivp(full_rhs,(0,T),np.r_[np.array(row['initial']),np.eye(3).ravel(),0],
                        args=(p,1),method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        z=sol.sol(np.linspace(0,T,n,endpoint=False))[:3].T
    np.savez_compressed(a.out,points=z)
    actual=hashlib.sha256(a.out.read_bytes()).hexdigest()
    print(json.dumps(dict(identifier=a.id,sha256=actual,original=row['sha256'],shape=z.shape)))


if __name__=='__main__':
    main()
