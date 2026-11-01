"""Fresh starts and nearby parameters for selected pre-existing tanh Yogi maps."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from scipy.stats import qmc
from yogi_tanh_map import spectra


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    with np.load(args.source) as f:s={k:f[k] for k in f.files}
    with np.load(args.reference) as f:r={k:f[k] for k in f.files}
    amplitudes=np.r_[-np.logspace(-6,2,10),np.logspace(-6,2,10),-.73,.73]
    factors=np.exp(.04*(qmc.Sobol(4,scramble=True,seed=2026110163).random_base2(5)-.5))
    volume=qmc.Sobol(3,scramble=True,seed=2026110164).random_base2(6)
    parameters=[];starts=[];bias=[];labels=[];boxes=[]
    for ident in args.ids:
        p=s['parameters'][ident];j=list(r['ids']).index(ident);z=r['tails'][j]
        extent=np.ptp(z,axis=0);low=z.min(0)-.1*extent;high=z.max(0)+.1*extent;low[2]=max(low[2],1e-12)
        boxes.append(dict(id=ident,lower=low.tolist(),upper=high.tolist()))
        for kind,biased in [('autonomous_cold',False),('library_style_cold',True)]:
            for amplitude in amplitudes:
                parameters.append(p.copy());starts.append([amplitude,0.,1e-5/p[3] if biased else 0.]);bias.append(biased);labels.append(dict(id=ident,kind=kind,amplitude=float(amplitude)))
        for index,factor in enumerate(factors):
            near=p.copy();near[:2]=1-(1-p[:2])*factor[:2];near[2:]*=factor[2:];assert np.all((near[:2]>0)&(near[:2]<1))
            biased=bool(s['bias'][ident]);parameters.append(near);starts.append([.73,0.,1e-5/near[3] if biased else 0.]);bias.append(biased);labels.append(dict(id=ident,kind='parameter_neighbor',index=index))
        for index,u in enumerate(volume):
            parameters.append(p.copy());starts.append(low+u*(high-low));bias.append(False);labels.append(dict(id=ident,kind='state_box',index=index))
    p=np.array(parameters);z=np.array(starts);bias=np.array(bias);Q=np.linalg.qr(np.array([[1.,2.,3.],[4.,5.,7.],[2.,6.,8.]]))[0]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),reference=str(args.reference),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),ids=args.ids,
        declared=len(labels),amplitudes=amplitudes.tolist(),neighbor_design='32Sobol4tuples seed2026110163; multiply(1-b,1-c,e,k) by exp(.04*(u-.5)), approximately+/-2%. Retain original source-row bias convention.',
        boxes=boxes,box_design='64Sobol3tuples seed2026110164 in reference-tail box expanded10%; lower v clipped to1e-12. Autonomous histories.',
        warmup=100000,observed=500000,window=50000,tail=4096,
        scope='Finite robustness and attraction evidence only. Library-style runs are independent transcriptions at scout coefficients, not actual TensorFlow executions. Reflection-related plus/minus starts are not independent discoveries.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for n in ['check_yogi_tanh_cohorts.py','yogi_tanh_map.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    before=time.perf_counter();w,det,tail,end,status,negative=spectra(p,z,bias,100000,500000,50000,4096,Q)
    np.savez_compressed(out/'histories.npz',parameters=p,starts=z,bias=bias,windows=w,determinant_windows=det,tails=tail,ends=end,status=status,negative_determinant_counts=negative,neighbor_factors=factors,box_points=volume)
    rows=[]
    for i,label in enumerate(labels):
        if status[i]:rows.append(dict(**label,status=int(status[i])));continue
        tol=1e-9*(1+np.max(abs(tail[i])));period=next((lag for lag in range(1,129) if np.max(abs(tail[i,lag:]-tail[i,:-lag]))<tol),None)
        lead=bool(w[i,:,0].min()>.005 and np.ptp(tail[i,:,0])>1e-4 and period is None)
        residual=float(np.max(abs(w[i].sum(1)-det[i])));assert residual<1e-7
        rows.append(dict(**label,status=0,spectrum=w[i].mean(0).tolist(),window_minima=w[i].min(0).tolist(),period_up_to128=period,leading=lead,two=bool(lead and w[i,:,1].min()>.001),determinant_sum_max_residual=residual))
    summary=[]
    for ident in args.ids:
        for kind in ['autonomous_cold','library_style_cold','parameter_neighbor','state_box']:
            chosen=[r for r in rows if r['id']==ident and r['kind']==kind]
            summary.append(dict(id=ident,kind=kind,total=len(chosen),leading=sum(r.get('leading',False) for r in chosen),two=sum(r.get('two',False) for r in chosen),periodic=sum(r.get('period_up_to128') is not None for r in chosen)))
            print(summary[-1],flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(completed=int(np.sum(status==0)),declared=len(labels),runs=rows,summary=summary,seconds=time.perf_counter()-before),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['source','reference','output']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',default=[39,7239,2850]);main(p.parse_args())
