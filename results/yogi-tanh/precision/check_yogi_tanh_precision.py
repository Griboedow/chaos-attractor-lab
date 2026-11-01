"""Independent finite-separation checks of smooth Yogi at70and90digits."""
import argparse,hashlib,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import mpmath as mp
import numpy as np


def run(job):
    ident,raw_p,raw_z,bias,digits,output=job;mp.mp.dps=digits
    b,c,e,k=map(lambda x:mp.mpf(float(x)),raw_p);z=tuple(mp.mpf(float(x)) for x in raw_z)
    initial=[mp.nstr(x,digits) for x in z];delta=mp.mpf('1e-30' if digits==70 else '1e-40');before=time.perf_counter()
    stem=Path(output)/f'case-{ident}-{digits}digits'
    def advance(w,factor=1):
        x,m,v=w;M=b*m+(1-b)*x;V=v+(1-c)*mp.tanh(k*(x*x-v))*x*x
        return (x-factor*M/(mp.sqrt(V)+e),M,V)
    pb=mp.mpf(1);pc=mp.mpf(1);cutoff=mp.mpf(10)**(-digits-10)
    for t in range(100000):
        if bias:
            pb=pb*b if pb>cutoff else mp.mpf(0);pc=pc*c if pc>cutoff else mp.mpf(0)
            factor=mp.sqrt(1-pc)/(1-pb)
        else:factor=1
        z=advance(z,factor)
    assert not bias or pb==pc==0
    nearby=[(z[0]+delta,z[1],z[2]),(z[0],z[1]+delta,z[2])];windows=[];acc=[mp.mpf(0),mp.mpf(0)];tail=[]
    for t in range(105000):
        z=advance(z);other=[advance(w) for w in nearby];basis=[];logs=[]
        for w in other:
            diff=[a-bb for a,bb in zip(w,z)]
            for _ in range(2):
                for q in basis:
                    dot=sum(a*bb for a,bb in zip(diff,q));diff=[a-dot*bb for a,bb in zip(diff,q)]
            length=mp.sqrt(sum(x*x for x in diff))
            if not length:raise RuntimeError('Finite separation collapsed')
            basis.append([x/length for x in diff]);logs.append(mp.log(length/delta))
        nearby=[tuple(a+delta*bb for a,bb in zip(z,q)) for q in basis]
        if t>=5000:
            acc=[a+bb for a,bb in zip(acc,logs)];n=t-5000+1
            if n%10000==0:
                windows.append([float(x/10000) for x in acc]);acc=[mp.mpf(0),mp.mpf(0)]
                stem.with_suffix('.progress.json').write_text(json.dumps(dict(observed=n,windows=windows,seconds=time.perf_counter()-before),indent=2)+'\n')
            if n>99000:tail.append([float(x) for x in z])
    np.savez_compressed(stem.with_suffix('.npz'),tail=tail,windows=windows)
    r=dict(id=ident,digits=digits,bias_warmup=bias,parameters_binary64=raw_p,parameters_exact_promoted=[mp.nstr(x,digits) for x in (b,c,e,k)],
        initial_exact_promoted=initial,separation=mp.nstr(delta,digits),spectrum_two=np.mean(windows,axis=0).tolist(),minimum_windows=np.min(windows,axis=0).tolist(),windows=windows,
        final_state=[mp.nstr(x,digits) for x in z],seconds=time.perf_counter()-before,tail_sha256=hashlib.sha256(stem.with_suffix('.npz').read_bytes()).hexdigest())
    stem.with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n');return r


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    with np.load(args.source) as f:d={k:f[k] for k in f.files}
    jobs=[(i,d['parameters'][i].tolist(),d['starts'][i].tolist(),bool(d['bias'][i]),digits,str(out)) for i in args.ids for digits in [70,90]]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=args.ids,digits=[70,90],warmup=100000,alignment=5000,observed=100000,window=10000,
        method='Independent scalar mpmath recurrence and two finite-separation directions, reorthogonalized twice each step. No production map or analytic derivative is imported.',
        initialization='Exact promotion of stored binary64 coefficients and coordinates; retain original declared bias warmup. Bias powers below10^(-digits-10) are set to0, after their subtraction has rounded to1 at working precision.',
        scope='Finite numerical reproducibility, not asymptotic Lyapunov positivity or historical uniqueness.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());rows=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            row=future.result();rows.append(row);print(row['id'],row['digits'],row['spectrum_two'],row['minimum_windows'],flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=sorted(rows,key=lambda x:(x['id'],x['digits'])),completed=len(rows),declared=len(jobs)),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())
