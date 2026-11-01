"""Predeclared cold-start and local parameter cohorts for selected Adamax cases."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.spatial import cKDTree
from scipy.stats import qmc
from adamax_map import spectra,trajectory


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    source=np.load(args.source);long=np.load(args.reference);ids=[847,239,4099,4096]
    amplitudes=np.r_[-np.logspace(-6,2,10),np.logspace(-6,2,10),-.73,.73]
    # A fixed relative perturbation of both memory loss rates stays inside (0,1).
    offsets=qmc.Sobol(2,scramble=True,seed=2026110123).random_base2(5)
    factors=np.exp(.04*(offsets-.5))
    parameters=[];starts=[];labels=[]
    for ident in ids:
        p=source['parameters'][ident]
        for amplitude in amplitudes:
            parameters.append(p.copy());starts.append([amplitude,0.,0.]);labels.append(dict(id=ident,kind='cold_start',amplitude=float(amplitude)))
        for index,factor in enumerate(factors):
            nearby=p.copy();nearby[:2]=1-(1-p[:2])*factor
            assert np.all((nearby[:2]>0)&(nearby[:2]<1))
            parameters.append(nearby);starts.append([.73,0.,0.]);labels.append(dict(id=ident,kind='parameter_neighbor',index=index))
    protocol=dict(source=str(args.source.resolve()),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        reference=str(args.reference.resolve()),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        ids=ids,cold_amplitudes=amplitudes.tolist(),neighbor_design='32 scrambled Sobol points, seed 2026110123, multiply each of (1-beta1,1-beta2) by exp(.04*(u-.5)). Same offsets for all four centers.',
        warmup=100000,observed=200000,window=50000,tail=4096,
        status_rule='Report every run, including nonpositive windows and period closures; no posthoc selection.',
        support_distance='Each candidate tail to 100000 reference tail states, Euclidean distance divided by reference extent norm; directed finite-sample diagnostic, not a basin proof.',
        zero_start='Exact absorbing origin is checked separately without tangent normalization.',
        scope='Finite cohorts, not an open parameter neighborhood, basin probability or asymptotic-chaos theorem.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_adamax_cohorts.py','adamax_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    before=time.perf_counter();parameters=np.array(parameters);starts=np.array(starts)
    windows,tails,ends,resets,margins=spectra(parameters,starts,100000,200000,50000,4096,True)
    np.savez_compressed(out/'histories.npz',parameters=parameters,starts=starts,windows=windows,tails=tails,ends=ends,reset_counts=resets,minimum_margins=margins,neighbor_factors=factors)
    references={}
    for ident in ids:
        k=list(long['ids']).index(ident);reference=long['tails'][k]
        references[ident]=(cKDTree(reference),np.linalg.norm(np.ptp(reference,axis=0)))
    rows=[]
    for i,label in enumerate(labels):
        valid=bool(np.isfinite(windows[i]).all() and np.isfinite(tails[i]).all())
        row=dict(**label,parameters=parameters[i].tolist(),completed_finite=valid)
        if valid:
            tail=tails[i,-2048:];returns=[float(np.max(np.linalg.norm(tail[lag:]-tail[:-lag],axis=1))) for lag in range(1,65)]
            period=next((k+1 for k,r in enumerate(returns) if r<1e-10*(1+np.linalg.norm(np.ptp(tail,axis=0)))),None)
            tree,scale=references[label['id']];dist=tree.query(tails[i])[0]/scale
            row.update(spectrum_two=windows[i].mean(0).tolist(),minimum_leading_window=float(windows[i,:,0].min()),
                leading_windows=windows[i,:,0].tolist(),tail_period_up_to64=period,tail_extent=np.ptp(tails[i],axis=0).tolist(),
                support_distances_q50_q95_q99_max=np.quantile(dist,[.5,.95,.99,1]).tolist(),
                reset_fraction=float(resets[i]/200000),minimum_scaled_switch_margin=float(margins[i]))
        rows.append(row)
    zero=[]
    for ident in ids:
        p=source['parameters'][ident];z=trajectory(p,np.zeros(3),0,10000)
        assert np.all(z[:,:2]==0) and np.all(z[:,2]==p[2])
        zero.append(dict(id=ident,steps=len(z),final_state=z[-1].tolist(),classification='Exact absorbing fixed state, not chaotic.'))
    result=dict(runs=rows,zero_controls=zero,completed=len(rows),declared=len(labels),seconds=time.perf_counter()-before)
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    for ident in ids:
        for kind in ['cold_start','parameter_neighbor']:
            selected=[r for r in rows if r['id']==ident and r['kind']==kind]
            passed=[r for r in selected if r.get('minimum_leading_window',-np.inf)>.005 and r['tail_period_up_to64'] is None]
            print(ident,kind,len(passed),'/',len(selected),'positive nonperiodic finite checks',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
