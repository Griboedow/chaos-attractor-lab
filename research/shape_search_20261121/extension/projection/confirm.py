"""Five-million-step QR, multistart controls, circular counterfactual and point clouds."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from affine_model import step,orbit,independent
from survey import BASE
from select_shapes import descriptors


def period(z):
    last=z[-384:];span=max(1.,np.linalg.norm(np.ptp(z,axis=0)))
    for n in range(1,257):
        if np.max(np.linalg.norm(last-z[-384-n:-n],axis=1))/span<1e-8:return n
    return None


def convert(raw):
    n=int(raw[0,3]);return dict(beta=float(raw[0,0]),inv_tau=float(raw[0,1]),tube_radius=float(raw[0,2]),scales=raw[2:5,3].tolist(),warp=float(raw[5,3]),sites=raw[2:2+n,:3].tolist())


def from_record(rec):
    p=rec['parameters'];raw=np.zeros((10,4));raw[0]=[p['beta'],p['inv_tau'],p['tube_radius'],len(p['sites'])]
    raw[1,:3]=rec['initial'];raw[2:2+len(p['sites']),:3]=p['sites'];raw[2:5,3]=p['scales'];raw[5,3]=p.get('warp',0)
    return raw


def run(args):
    with np.load(args.shortlist) as f:
        ids=f['ids'];raws=f['raw'];parents=f['parents'];previous=f['refs'];spread=f['spread']
    args.output.mkdir(parents=True,exist_ok=True);args.clouds.mkdir(parents=True,exist_ok=True)
    for j in range(args.rank,len(ids),args.workers):
        case=int(ids[j]);raw=raws[:,j,:].astype(float);start=raw[1,:3].copy()
        record=args.output/f'case-{case}.json'
        if record.exists():continue
        quick,rs=orbit(raw,start,65536,32768,4,3072)
        p=period(quick)
        if min(rs[:,2])<.002 or p is not None:
            record.write_text(json.dumps(dict(case=case,shortlist_position=j,retained=False,phase='cpu rescreen',parent=int(parents[j]),parameters=convert(raw),initial=start.tolist(),qr=rs.tolist(),period=p),indent=2)+'\n',encoding='utf8')
            continue
        z,rs=orbit(raw,start,100000,500000,10,250000)
        delta=[];values=[]
        for x in z[::32700]:
            actual,jac=step(x,raw);h=2e-6
            approx=np.column_stack([(independent(x+h*d,raw)-independent(x-h*d,raw))/(2*h) for d in np.eye(3)])
            delta.append(np.linalg.norm(jac-approx)/max(1.,np.linalg.norm(jac)))
            values.append(np.max(abs(actual-independent(x,raw))))
        cohort=[];rng=np.random.default_rng(508922+case)
        for x in np.vstack((start,rng.uniform(-2,2,(16,3)))):
            tail,sr=orbit(raw,x,32768,16384,4,1024);pr=period(tail)
            cohort.append(dict(initial=x.tolist(),windows=sr.tolist(),period=pr,three_positive=bool(sr.min()>.001 and pr is None)))
        hits=sum(t['three_positive'] for t in cohort);pr=period(z)
        circle=raw.copy();circle[2:5,3]=1;circle[5,3]=0
        cz,cr=orbit(circle,start,65536,32768,4,3072)
        samples=np.stack([descriptors(z[offset:offset+1536]) for offset in (20000,65000,155000,204000)])
        old=previous.reshape(31,4,-1).mean(1)
        target=samples.mean(0)
        norm=lambda v:np.sqrt(np.mean((v/spread)**2,axis=-1))
        matches=norm(old-target);ix=int(np.argmin(matches))
        ref_spread=norm(previous.reshape(31,4,-1)[ix]-old[ix])
        same=norm(samples-target)
        mat=np.linalg.eigvalsh(np.cov(z.T));mat/=mat.sum()
        pass_all=bool(rs.min()>.001 and hits>=9 and pr is None and max(delta)<1e-5)
        extras={}
        if pass_all:
            dest=args.clouds/f'affine-{case}.npz'
            np.savez_compressed(dest,points=z.astype('f4'),case=case,initial=start)
            extras=dict(cloud_relative='extension/projection/clouds/'+dest.name,cloud_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),points_saved=len(z))
        doc=dict(case=case,shortlist_position=j,retained=pass_all,phase='long validation',parent=int(parents[j]),family='Quadratically sheared torus' if raw[5,3] else 'Affine-torus metric projection',parameters=convert(raw),initial=start.tolist(),burn=100000,windows=10,window_length=500000,rates=rs.mean(0).tolist(),minimum=rs.min(0).tolist(),qr_windows=rs.tolist(),period=pr,cohort_three=hits,cohort=cohort,cohort_burn=32768,cohort_window=16384,jacobian_finite_difference_max_relative=max(delta),nested_map_max_error=max(values),pca_fractions=mat.tolist(),normalized_descriptor_distance_to_nearest_old_mean=float(matches[ix]),nearest_old_index=ix,within_candidate_max=float(same.max()),within_old_max=float(ref_spread.max()),circular_counterfactual=dict(rates=cr.mean(0).tolist(),minimum=cr.min(0).tolist(),period=period(cz)),**extras)
        tmp=record.with_suffix('.tmp.json');tmp.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8');os.replace(tmp,record)
        print(case,pass_all,hits,np.round(doc['rates'],4),round(matches[ix],3),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--shortlist',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--clouds',type=Path,required=True)
    p.add_argument('--rank',type=int,default=0);p.add_argument('--workers',type=int,default=1)
    run(p.parse_args())
