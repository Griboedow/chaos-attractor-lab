"""Track the 16T C7-relative orbit through its next multiplier crossing."""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from branches import rotation,polar_ring
from relative_cores import solve_sector,sector_cloud
from field import integrate,ring_tangent
from periodic_catalog import shifts
from scipy.integrate import solve_ivp


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--old',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--parent32',type=Path,help='Instead continue the archived 32T parent and test for a 64T daughter')
    a=ap.parse_args()
    if (a.out/'report.json').exists():return
    a.out.mkdir(exist_ok=False,parents=True)
    old=json.loads((a.old/'long-doubling.json').read_text());q=old['accepted_cycles'][0]
    if a.parent32:
        q=next(rr for rr in json.loads(a.parent32.read_text())['postflip'] if rr.get('cloud_relative'))
        p=np.array(q['parameters']);v=np.array(q['initial']);t=q['period']/7;parent_k=1
        grid=np.arange(q['f'],.9532,.000025);offsets=[6e-6,1.3e-5,2.1e-5]
    else:
        p=np.array(old['parameters']);v=np.array(q['refined_start']);t=q['period']/7;parent_k=4
        grid=np.arange(.9498,.95505,.0001);offsets=[3e-5,7e-5,.00012,.00018,.00026,.00032,.0004]
    child_k=2*parent_k%7;target_period=64 if a.parent32 else 32
    e0=-.4297716347835525;inc=.0759375;R=rotation(2*np.pi*parent_k/7)
    log=[];candidates=[];prev=None
    for f in grid:
        p[8]=e0+inc*f;got=solve_sector(v,t,p,R)
        if got is None:continue
        v,t=got;_,stats=sector_cloud(v,t,p,parent_k,56);mu=stats['dominant_multiplier'][0]
        log.append(dict(f=f,e=p[8],period=t*7,multiplier=mu,seed=v.tolist(),closure=stats['sector_return_error_into_rotated_center']))
        if prev is not None and prev>-1>=mu:
            pair=[log[-2],log[-1]];break
        prev=mu
    else:pair=[]
    result=dict(target_period_order=target_period,tracked=log,first_flip_bracket=pair,postflip=[])
    save=lambda:(a.out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    save();print('16T flip',[(c['f'],c['multiplier']) for c in pair],flush=True)
    if not pair:return
    low,hi=pair;left=np.array(low['seed']);tleft=low['period']/7
    for _ in range(16):
        f=(low['f']+hi['f'])/2;p[8]=e0+inc*f
        v,t=solve_sector(left,tleft,p,R);_,st=sector_cloud(v,t,p,parent_k,56)
        row=dict(f=f,e=p[8],period=t*7,multiplier=st['dominant_multiplier'][0]);s=st['dominant_multiplier'][0]
        if s>-1:low=row;left=v;tleft=t
        else:hi=row
    result['refined_flip_bracket']=[low,hi];save()
    fflip=(low['f']+hi['f'])/2
    for delta in offsets:
        f=fflip+delta;p[8]=e0+inc*f
        h=.008;start=left+np.array([.0002,-.00014,.00021])
        rates,z=integrate(p,start,h,round(29000/h),round(750/h),4,4,ring_tangent);T0=2*tleft
        guesses=[]
        for i in np.linspace(0,len(z)-int(T0/(h*4))-180,70).astype(int):
            center=i+round(T0/h/4);near=z[center-30:center+31];err=np.linalg.norm(near-rotation(2*np.pi*child_k/7)@z[i],axis=1)
            j=np.argmin(err);guesses.append((float(err[j]),i,(center-30+j-i)*h*4))
        row=dict(f=f,e=p[8],offset_from_flip=delta,rates=rates.tolist());result['postflip'].append(row)
        for err,i,tguess in sorted(guesses)[:6]:
            if err>.003:break
            got=solve_sector(z[i],tguess,p,rotation(2*np.pi*child_k/7))
            if got is None:continue
            y,tt=got;pts,st=sector_cloud(y,tt,p,child_k,65536);floor,div=shifts(pts)
            if st['sector_return_error_into_rotated_center']>2e-7 or abs(st['dominant_multiplier'][0])>=.995 or floor<.0001:continue
            check=solve_ivp(lambda _,q:polar_ring(q,p),(0,tt),y,method='RK45',rtol=1e-12,atol=1e-14)
            row.update(period=7*tt,multiplier=st['dominant_multiplier'][0],initial=y.tolist(),sector_period=tt,
                       rotation_index=child_k,primitive_floor=floor,worst_divisor=div,RK45_closure=float(np.linalg.norm(check.y[:,-1]-rotation(2*np.pi*child_k/7)@y)),**{'DOP853_closure':st['sector_return_error_into_rotated_center']})
            if not any(x.get('cloud_relative') for x in result['postflip']):
                path=a.out/f'{target_period}T.npz';np.savez_compressed(path,points=pts)
                row.update(cloud_relative=path.resolve().relative_to(a.out.resolve().parents[2]).as_posix(),path=str(path.resolve()),shape=list(pts.shape),point_key='points',parameters=p.tolist(),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            break
        save();print(delta,row.get('period'),row.get('multiplier'),flush=True)


if __name__=='__main__':
    main()
