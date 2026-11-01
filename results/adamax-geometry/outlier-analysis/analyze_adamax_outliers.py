"""Analyze frozen outlier trajectories; deduplicate distance queries exactly."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.spatial import cKDTree
from adamax_map import step,step_tangent


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    d=np.load(args.source/'histories.npz');tails=d['tails'];w=d['windows'];protocol=json.loads((args.source/'protocol.json').read_text())
    original=Path(protocol['source']);approach=np.load(original/'case-4099.npz');volume=np.load(original/'volume-checks.npz')
    far=np.flatnonzero(approach['endpoint_distances']>.01);weak=np.flatnonzero(np.any(volume['windows'][:,:,0]<=.005,axis=1))
    labels=[dict(kind='far_endpoint',index=int(i)) for i in far]+[dict(kind='weak_volume_check',index=int(i)) for i in weak]
    reference=np.load(args.reference,mmap_mode='r');tree=cKDTree(reference);scale=np.linalg.norm(np.ptp(reference,axis=0))
    unique,inverse=np.unique(tails.reshape(-1,3),axis=0,return_inverse=True)
    print('Unique distance queries',len(unique),'of',len(inverse),flush=True)
    distances=(tree.query(unique)[0]/scale)[inverse].reshape(tails.shape[:2]);rows=[]
    for i,label in enumerate(labels):
        tail=tails[i,-1024:];returns=[float(np.max(np.linalg.norm(tail[k:]-tail[:-k],axis=1))) for k in range(1,257)]
        period=next((k+1 for k,r in enumerate(returns) if r<1e-9*(1+np.linalg.norm(np.ptp(tail,axis=0)))),None)
        rows.append(dict(**label,spectrum_two=w[i].mean(0).tolist(),minimum_leading_window=float(w[i,:,0].min()),tail_period_up_to256=period,
            minimum_return_residual=min(returns),minimum_return_lag=int(np.argmin(returns)+1),
            support_distance_q50_q99_max=np.quantile(distances[i],[.5,.99,1]).tolist()))
    # Close the first return satisfying the declared tolerance, not the first
    # bitwise-exact floating-point repetition (which may be a multiple).
    p=d['parameters'][0];period=rows[0]['tail_period_up_to256'];assert period is not None
    def return_map(z):
        for _ in range(period):z=step(z,p)
        return z
    sol=root(lambda z:return_map(z)-z,d['ends'][0],tol=1e-11);z=sol.x.copy();orbit=[];product=np.eye(3);margins=[]
    for _ in range(period):
        orbit.append(z.copy());margins.append(abs(abs(z[0])+p[2]-p[1]*z[2]));z,J,_=step_tangent(z,np.eye(3),p);product=J@product
    closure=float(np.max(abs(z-sol.x)));eig=np.linalg.eigvals(product)
    assert closure<1e-10 and max(abs(eig))<1 and min(margins)>1e-8
    orbit=np.array(orbit);distinct=float(min(np.linalg.norm(orbit[k]-orbit[0]) for k in range(1,period)));assert distinct>1e-5
    cycle_dist=np.max(cKDTree(orbit).query(tails.reshape(-1,3))[0])
    cycle=dict(period=period,parameters=p.tolist(),closure=closure,solver_success=bool(sol.success),minimum_other_phase_distance=distinct,
        minimum_switch_margin=float(min(margins)),multipliers=[[float(v.real),float(v.imag)] for v in eig],spectral_radius=float(max(abs(eig))),
        maximum_distance_all437_tails_to_closed_cycle=float(cycle_dist),scope='Floating-point closure and stability only, not an interval existence certificate.')
    np.savez_compressed(out/'cycle.npz',orbit=orbit,monodromy=product,parameters=p,unique_distance_states=unique,distances=distances)
    result=dict(runs=rows,completed=len(rows),declared=len(labels),positive_all_windows=sum(r['minimum_leading_window']>.005 for r in rows),
        closed_period_counts={str(k):sum(r['tail_period_up_to256']==k for r in rows) for k in sorted({r['tail_period_up_to256'] for r in rows if r['tail_period_up_to256'] is not None})},
        cycle=cycle,unique_query_count=len(unique),total_query_count=len(inverse),
        source_sha256=hashlib.sha256((args.source/'histories.npz').read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest())
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in ['analyze_adamax_outliers.py','adamax_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(result['closed_period_counts'],cycle,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
