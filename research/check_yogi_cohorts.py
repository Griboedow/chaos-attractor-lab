"""Cold starts, nearby parameters and full-state boxes for four Yogi candidates."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.spatial import cKDTree
from scipy.stats import qmc
from yogi_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    source=np.load(args.source);reference=np.load(args.reference);ids=[44,7590,4214,5174]
    amplitudes=np.r_[-np.logspace(-6,2,10),np.logspace(-6,2,10),-.73,.73]
    factors=np.exp(.04*(qmc.Sobol(3,scramble=True,seed=2026110142).random_base2(5)-.5))
    volume=qmc.Sobol(3,scramble=True,seed=2026110143).random_base2(6)
    parameters=[];starts=[];labels=[];boxes=[];references={}
    for ident in ids:
        p=source['parameters'][ident];k=int(np.flatnonzero(np.all(reference['parameters']==p,axis=1))[0]);cloud=reference['tails'][k]
        extent=np.ptp(cloud,axis=0);lower=cloud.min(0)-.1*extent;upper=cloud.max(0)+.1*extent;lower[2]=max(lower[2],1e-8)
        boxes.append(dict(id=ident,lower=lower.tolist(),upper=upper.tolist(),reference_id=int(reference['ids'][k])))
        references[ident]=(cKDTree(cloud),np.linalg.norm(extent))
        for amplitude in amplitudes:
            parameters.append(p.copy());starts.append([amplitude,0.,0.]);labels.append(dict(id=ident,kind='cold_start',amplitude=float(amplitude)))
        for index,factor in enumerate(factors):
            near=p.copy();near[:2]=1-(1-p[:2])*factor[:2];near[2]*=factor[2]
            assert np.all((near[:2]>0)&(near[:2]<1))
            parameters.append(near);starts.append([.73,0.,0.]);labels.append(dict(id=ident,kind='parameter_neighbor',index=index))
        for index,u in enumerate(volume):
            parameters.append(p.copy());starts.append(lower+u*(upper-lower));labels.append(dict(id=ident,kind='state_box',index=index))
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        reference=str(args.reference),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),ids=ids,
        cold_amplitudes=amplitudes.tolist(),parameter_design='32 Sobol triples, seed2026110142; multiply(1-b,1-c,e) by exp(.04*(u-.5)), approximately+/-2%.',
        state_boxes=boxes,box_design='64 Sobol triples, seed2026110143, in full3D reference bounding boxes expanded10% on both sides; lower v clipped to1e-8.',
        warm=100000,observed=200000,window=50000,tail=4096,bias_correction=False,
        interpretation='Every declared run retained. Positive rule: all leading windows>.005,tail x span>1e-4,no lag<=128 closure. Second-direction positivity reported separately.',
        support_distance='Directed Euclidean nearest-neighbor distance to100000 reference states divided by reference extent norm; finite sampled support diagnostic.',
        scope='Finite cohorts, not basin probabilities, an open chaotic region or proof. Cold +/- pairs are exactly symmetry-related, not independent attractor discoveries.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for n in ['check_yogi_cohorts.py','yogi_map.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    before=time.perf_counter();parameters=np.array(parameters);starts=np.array(starts)
    windows,tails,ends,margins=spectra(parameters,starts,100000,200000,50000,4096)
    assert np.isfinite(windows).all() and np.isfinite(tails).all()
    np.savez_compressed(out/'histories.npz',parameters=parameters,starts=starts,windows=windows,tails=tails,ends=ends,minimum_margins=margins,neighbor_factors=factors,volume_points=volume)
    rows=[]
    for i,label in enumerate(labels):
        tail=tails[i];tol=1e-9*(1+np.max(abs(tail)))
        period=next((lag for lag in range(1,129) if np.max(abs(tail[lag:]-tail[:-lag]))<tol),None)
        unique,inverse=np.unique(tail,axis=0,return_inverse=True);tree,scale=references[label['id']];distance=tree.query(unique)[0][inverse]/scale
        passed=bool(windows[i,:,0].min()>.005 and np.ptp(tail[:,0])>1e-4 and period is None)
        row=dict(**label,parameters=parameters[i].tolist(),spectrum=windows[i].mean(0).tolist(),
            window_minima=windows[i].min(0).tolist(),period_up_to128=period,tail_extents=np.ptp(tail,axis=0).tolist(),
            support_distances_q50_q95_q99_max=np.quantile(distance,[.5,.95,.99,1]).tolist(),passed_leading=passed,
            passed_two=bool(passed and windows[i,:,1].min()>.001),minimum_scaled_switch_margin=float(margins[i]))
        rows.append(row)
    summary=[]
    for ident in ids:
        for kind in ['cold_start','parameter_neighbor','state_box']:
            selected=[r for r in rows if r['id']==ident and r['kind']==kind]
            record=dict(id=ident,kind=kind,total=len(selected),passed_leading=sum(r['passed_leading'] for r in selected),
                passed_two=sum(r['passed_two'] for r in selected),periodic=sum(r['period_up_to128'] is not None for r in selected))
            summary.append(record);print(record,flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows,summary=summary,completed=len(rows),declared=len(labels),wall_seconds=time.perf_counter()-before),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
