"""Exploratory long-horizon excursion audit with exact checkpoint replay.

The reference box comes from the published two-million-state portrait. Its
enlargement is a descriptive threshold, not an invariant neighborhood.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from search_product_alphabet import step


@njit(cache=True)
def inspect(x, beta, tau, levels, order, count, lower, upper, stride):
    sampled=np.empty((count//stride,3));indices=[np.int64(0) for _ in range(0)]
    events=[np.zeros(3) for _ in range(0)]; radii=[0. for _ in range(0)]
    minimum=1e100;minimum_at=-1;minimum_state=x.copy();peak=0.;peak_at=-1;peak_state=x.copy()
    for t in range(count):
        y,_,r,_,_,flag=step(x,beta,tau,levels,order)
        if flag: return x,sampled,np.array(indices),np.empty((0,3)),np.array(radii),minimum,minimum_at,minimum_state,peak,peak_at,peak_state,flag
        excess=max(np.max(lower-y),np.max(y-upper),0.)
        if excess>0: indices.append(t+1);events.append(y.copy());radii.append(r)
        if r<minimum: minimum=r;minimum_at=t;minimum_state=x.copy()
        if excess>peak: peak=excess;peak_at=t+1;peak_state=y.copy()
        x=y
        if (t+1)%stride==0:sampled[(t+1)//stride-1]=x
    cloud=np.empty((len(events),3))
    for i in range(len(events)):cloud[i]=events[i]
    return x,sampled,np.array(indices),cloud,np.array(radii),minimum,minimum_at,minimum_state,peak,peak_at,peak_state,0


@njit(cache=True)
def segment(x,beta,tau,levels,order,skip,count):
    for _ in range(skip):x=step(x,beta,tau,levels,order)[0]
    states=np.empty((count+1,3));states[0]=x
    for t in range(count):x=step(x,beta,tau,levels,order)[0];states[t+1]=x
    return states


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);start=time.time()
    names=[Path(__file__).name,'search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']
    for name in names:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    protocol=json.loads((args.horizon/'protocol.json').read_text(encoding='utf-8'));k=protocol['ids'].index(args.case)
    params=np.array(protocol['params'][k]);levels=np.array(protocol['levels'][k]);order=protocol['orders'][k]
    with np.load(args.reference) as f:z=f['states'];reference_lower=z.min(0);reference_upper=z.max(0)
    lower=reference_lower-args.pad;upper=reference_upper+args.pad
    blocks=sorted(args.horizon.glob('block-*.json'));rows=[];clouds=[];all_events=[];all_indices=[];event_segments=[]
    for file in blocks:
        b=json.loads(file.read_text());initial=np.array(b['incoming_state'][k]);count=b['count']
        x,samples,indices,events,radii,minr,at,state,peak,peak_at,peak_state,flag=inspect(initial,*params,levels,order,count,lower,upper,args.stride)
        assert not flag and np.array_equal(x,np.array(b['outgoing_state'][k])),file
        assert minr==b['minimum_radius'][k]
        row=dict(block=b['block'],start=b['start'],count=count,exceedance_count=len(indices),minimum_radius=float(minr),minimum_radius_at=b['start']+int(at),minimum_radius_state=state.tolist(),maximum_excess=float(peak),maximum_excess_at=b['start']+int(peak_at),maximum_excess_state=peak_state.tolist(),checkpoint_exact=True)
        rows.append(row);clouds.append(samples);all_events.append(events);all_indices.append(indices+b['start'])
        # Capture the block's normalization minimum with 32 preceding updates.
        # Starting from this exact binary64 state permits a short independent
        # arbitrary-precision replay without claiming a 100-million-step shadow.
        left=max(0,int(at)-32);take=min(count-left,288)
        seq=segment(initial,*params,levels,order,left,take)
        np.savez_compressed(args.out/f'block-{b["block"]:03}-minimum.npz',states=seq,start=b['start']+left,minimum_at=int(at)-left,params=params,levels=levels,order=order)
        event_segments.append(dict(block=b['block'],start=b['start']+left,minimum_offset=int(at)-left,source=f'block-{b["block"]:03}-minimum.npz'))
        (args.out/'progress.json').write_text(json.dumps(dict(completed=b['start']+count,exceedances=sum(r['exceedance_count'] for r in rows),minimum_radius=min(r['minimum_radius'] for r in rows),seconds=time.time()-start),indent=2)+'\n',encoding='utf-8',newline='\n')
        print(json.dumps(row),flush=True)
    events=np.concatenate(all_events);indices=np.concatenate(all_indices)
    groups=1+np.sum(np.diff(indices)>args.event_gap) if len(indices) else 0
    np.savez_compressed(args.out/'uniform-sample.npz',states=np.concatenate(clouds),params=params,levels=levels,order=order,stride=args.stride,first_observed_update=args.stride)
    np.savez_compressed(args.out/'all-threshold-exceedances.npz',states=events,indices=indices,params=params,levels=levels,order=order)
    report=dict(case=args.case,horizon_protocol_sha256=hashlib.sha256((args.horizon/'protocol.json').read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),reference_lower=reference_lower.tolist(),reference_upper=reference_upper.tolist(),padding=args.pad,threshold_lower=lower.tolist(),threshold_upper=upper.tolist(),observed=sum(r['count'] for r in rows),exceedance_count=len(indices),fraction=len(indices)/sum(r['count'] for r in rows),group_gap=args.event_gap,groups=int(groups),uniform_stride=args.stride,rows=rows,minimum_segments=event_segments,seconds=time.time()-start,scope='Exploratory threshold selected after inspecting extended-run bounds. Every block endpoint and normalization minimum replay exactly. Excursion counts are finite-orbit observations, not probabilities, a crisis classification or an asymptotic invariant-set result.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['horizon','reference','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--case',type=int,default=988);p.add_argument('--pad',type=float,default=.02)
    p.add_argument('--stride',type=int,default=50);p.add_argument('--event-gap',type=int,default=128)
    main(p.parse_args())
