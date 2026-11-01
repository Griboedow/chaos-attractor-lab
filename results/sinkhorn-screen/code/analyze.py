"""Classify retained final windows without equating recurrence with a proof."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import model
from spectrum import state_jacobian
import precision
import mpmath as mp

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();rows=[]
    for n in (3,4):
        data=json.loads((args.input/f'n{n}.json').read_text());a=np.load(args.input/f'n{n}-final-25000-steps.npz')
        for r,z in zip(data['runs'],a['states']):
            assert int(a['case'][len([x for x in rows if x['shape'][0]==n])])==r['case']
            radius=float(np.max(abs(z)));tolerance=1e-9*(1+radius);period=None;error=None
            # The entire final 4096 samples must match, not only the endpoint.
            for p in range(1,2049):
                e=float(np.max(abs(z[-4096:]-z[-4096-p:-p])))
                if e<tolerance:period=p;error=e;break
            centered=z-z.mean(0);sv=np.linalg.svd(centered,compute_uv=False)
            scaled=sv/sv[0] if sv[0]>0 else np.zeros_like(sv)
            variance=(scaled*scaled/np.sum(scaled*scaled)).tolist() if sv[0]>0 else [0.]*len(sv)
            item=dict(case=r['case'],shape=[n,n],stop_step=r['blocks'][-1]['stop_step'],last_block_rate=r['blocks'][-1]['rate'],
                final_window_max_abs=radius,recurrence_threshold=tolerance,recurrence_test_samples=4096,
                smallest_tested_recurrence=period,recurrence_max_error=error,pca_variance_fraction=variance)
            if period is not None and period<=256:
                L=np.array(r['log_kernel']);la,lb=model.targets(L);w=r['omega'];point=z[-1].copy();J=np.eye(len(point))
                for _ in range(period):point,A=state_jacobian(point,L,la,lb,w);J=A@J
                item['float_cycle_residual']=float(np.max(abs(point-z[-1])))
                item['float_multiplier_moduli']=np.abs(np.linalg.eigvals(J)).tolist()
                mp.mp.dps=100;K,aa,bb,ww=precision.prepare(r);point=[precision.exact(v) for v in z[-1]];v=[mp.mpf(0)]*len(point)
                for _ in range(period):point,v=precision.step(point,v,K,aa,bb,ww)
                item['mp100_cycle_residual']=mp.nstr(max(abs(u-precision.exact(vv)) for u,vv in zip(point,z[-1])),30)
                item['interpretation']='Numerical fixed-state or periodic recurrence; not an interval certificate.'
            else:item['interpretation']='No recurrence up to period 2048 at the declared tolerance; neither chaos nor attraction follows.'
            rows.append(item)
    result=dict(runs=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([dict(case=r['case'],shape=r['shape'],period=r['smallest_tested_recurrence'],rate=r['last_block_rate']) for r in rows]))
