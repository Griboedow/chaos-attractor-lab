"""Independent two-direction finite-separation Yogi checks at 70/90 digits."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,json,time
from pathlib import Path
import mpmath as mp
import numpy as np


def run_one(job):
    ident,raw_parameters,raw_start,digits,output=job;mp.mp.dps=digits
    b,c,e=[mp.mpf(float(v)) for v in raw_parameters];z=tuple(mp.mpf(float(v)) for v in raw_start)
    delta=mp.mpf('1e-30' if digits==70 else '1e-40');before=time.perf_counter()
    stem=Path(output)/f'case-{ident}-{digits}digits'
    def advance(w):
        x,m,v=w;mn=b*m+(1-b)*x;vn=v-(1-c)*mp.sign(v-x*x)*x*x
        return (x-mn/(mp.sqrt(vn)+e),mn,vn)
    for _ in range(100000):z=advance(z)
    nearby=[(z[0]+delta,z[1],z[2]),(z[0],z[1]+delta,z[2])]
    windows=[];total=[mp.mpf(0),mp.mpf(0)];tail=[];minimum=mp.inf
    for t in range(105000):
        margin=abs(z[2]-z[0]*z[0])/(1+z[2]);z=advance(z);next_nearby=[advance(w) for w in nearby]
        basis=[];growth=[]
        for w in next_nearby:
            difference=[a-bb for a,bb in zip(w,z)]
            for _ in range(2):
                for q in basis:
                    dot=sum(a*bb for a,bb in zip(difference,q))
                    difference=[a-dot*bb for a,bb in zip(difference,q)]
            length=mp.sqrt(sum(v*v for v in difference))
            if length==0:raise RuntimeError('Finite separation collapsed')
            basis.append([v/length for v in difference]);growth.append(mp.log(length/delta))
        nearby=[tuple(a+delta*bb for a,bb in zip(z,q)) for q in basis]
        if t>=5000:
            minimum=min(minimum,margin);total=[a+bb for a,bb in zip(total,growth)];n=t-5000+1
            if n%10000==0:
                windows.append([float(v/10000) for v in total]);total=[mp.mpf(0),mp.mpf(0)]
                stem.with_suffix('.progress.json').write_text(json.dumps(dict(observed=n,windows=windows,wall_seconds=time.perf_counter()-before),indent=2)+'\n')
            if n>99000:tail.append([float(v) for v in z])
    row=dict(id=int(ident),digits=digits,parameters_binary64=raw_parameters,start_binary64=raw_start,
        parameters_exact_promoted=[mp.nstr(v,digits) for v in (b,c,e)],separation=mp.nstr(delta,digits),
        spectrum_two=np.mean(windows,axis=0).tolist(),windows=windows,minimum_windows=np.min(windows,axis=0).tolist(),
        minimum_scaled_switch_margin=float(minimum),final_state=[mp.nstr(v,digits) for v in z],seconds=time.perf_counter()-before)
    np.savez_compressed(stem.with_suffix('.npz'),tail=tail,windows=windows)
    row['tail_sha256']=hashlib.sha256(stem.with_suffix('.npz').read_bytes()).hexdigest()
    stem.with_suffix('.json').write_text(json.dumps(row,indent=2)+'\n');return row


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);source=np.load(args.source)
    jobs=[(i,source['parameters'][i].tolist(),source['starts'][i].tolist(),d,str(out)) for i in args.ids for d in (70,90)]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=args.ids,digits=[70,90],
        warmup=100000,alignment=5000,observed=100000,window=10000,
        method='Independent scalar recurrence and two finite-separation directions; repeated Euclidean Gram-Schmidt and renormalization each step. No analytic Jacobian or production-map import.',
        parameter_semantics='Exact promotion of the selected binary64 coefficients and initial states; same real values at both precisions.',
        scope='Finite numerical evidence for leading and second growth, including discontinuity checks; not asymptotic chaos or historical uniqueness.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(run_one,job) for job in jobs]):
            row=future.result();rows.append(row);print(row['id'],row['digits'],row['spectrum_two'],row['minimum_windows'],flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=sorted(rows,key=lambda r:(r['id'],r['digits'])),completed=len(rows),declared=len(jobs)),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())
