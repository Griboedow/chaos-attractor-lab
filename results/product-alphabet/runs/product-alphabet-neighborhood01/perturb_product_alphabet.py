"""Declared parameter neighborhoods of the two selected product maps."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from numba import set_num_threads
from search_product_alphabet import follow


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4);begin=time.time()
    for name in [Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    original=json.loads(args.protocol.read_text(encoding='utf-8'))
    u=2*qmc.Sobol(8,scramble=True,seed=2026110147).random_base2(6)-1
    rows=[];params=[];levels=[];orders=[]
    for ident in [988,489]:
        source=original['cases'][ident];a=np.array(source['halfwidths']);c=np.array(source['center'])
        for radius in [.01,.05]:
            for k,v in enumerate(u):
                # Beta is additive relative to its distance from one. Noise and
                # halfwidths change multiplicatively; shifts use halfwidth units.
                beta=source['beta']+radius*(1-source['beta'])*v[0]
                tau=source['tau']*np.exp(radius*v[1]);aa=a*np.exp(radius*v[2:5]);cc=c+radius*a*v[5:8]
                n=source['order'];ll=np.array([cc[j]+aa[j]*np.linspace(-1,1,n) for j in range(3)])
                params.append([beta,tau]);levels.append(np.pad(ll,((0,0),(0,4-n))));orders.append(n)
                rows.append(dict(id=len(rows),base=ident,radius=radius,sobol_row=k,beta=beta,tau=tau,halfwidths=aa.tolist(),center=cc.tolist()))
    params=np.array(params);levels=np.array(levels);orders=np.array(orders)
    for offset in range(0,len(rows),8):
        sl=slice(offset,offset+8);rates,states,status,minima=follow(params[sl],levels[sl],orders[sl],100000,500000)
        np.savez_compressed(args.out/f'block-{offset:03}.npz',params=params[sl],levels=levels[sl],orders=orders[sl],rates=rates,states=states,status=status,minimum_radius=minima)
        for j in range(len(rates)):
            row=rows[offset+j];row.update(status=status[j].tolist(),minimum_radius=float(minima[j]))
            assert status[j].tolist()==[0,500000]
            z=states[j,-2048:];scale=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
            residuals=[float(np.max(np.linalg.norm(z[k:]-z[:-k],axis=1))) for k in range(1,513)]
            periods=[k+1 for k,r in enumerate(residuals) if r<1e-8*scale]
            row.update(spectrum=rates[j].mean(0).tolist(),minimum_leading=float(rates[j,:,0].min()),period=min(periods) if periods else None,return_residuals=residuals)
        print(json.dumps(dict(completed=offset+8,total=len(rows),seconds=time.time()-begin)),flush=True)
    report=dict(source_sha256=hashlib.sha256(args.protocol.read_bytes()).hexdigest(),seed=2026110147,initial=[.317,-.211,.419],warm=100000,observed=500000,windows=20,design='64 Sobol directions per base, repeated at radius 0.01 and 0.05. See executed source for exact parameter scaling.',cases=rows,scope='Finite parameter sensitivity from the declared original start; no open-set robustness or basin-probability proof.',seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--protocol',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
