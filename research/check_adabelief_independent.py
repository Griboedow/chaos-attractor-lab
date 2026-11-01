"""Unmodified upstream PyTorch and arbitrary-precision secant controls."""
import argparse
import hashlib
import importlib.util
import json
import time
from pathlib import Path
import numpy as np
from adabelief_map import step, jacobian


def author_check(source, output):
    import torch
    spec=importlib.util.spec_from_file_location('upstream_adabelief',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    reports=[]
    for b,c in [(.52,.15),(.9,.99),(.9,.999)]:
        p=torch.tensor([.0073],dtype=torch.float64,requires_grad=True)
        opt=module.AdaBelief([p],lr=.01,betas=(b,c),eps=1e-16,weight_decay=0,
                            weight_decouple=False,amsgrad=False,rectify=False,print_change_log=False)
        z=np.array([.73,0.,0.]);short_error=0.;acc=0.;windows=[];tail=[]
        v=np.array([1.,.3,.2]);v/=np.linalg.norm(v)
        minimum=np.full(3,np.inf);maximum=-minimum
        for t in range(120000):
            old_state=opt.state[p]
            if t>=20000:
                actual=np.array([p.item()/.01,old_state['exp_avg'].item()/.01,
                                 old_state['exp_avg_var'].item()/1e-4])
                v=jacobian(actual,b,c)@v
                norm=np.linalg.norm(v);v/=norm;acc+=np.log(norm)
            p.grad=p.detach().clone()
            opt.step()
            state=opt.state[p]
            actual=np.array([p.item()/.01,state['exp_avg'].item()/.01,state['exp_avg_var'].item()/1e-4])
            if t<25:
                z=step(z,b,c,iteration=t+1)
                short_error=max(short_error,float(abs(actual-z).max()))
            if t>=20000:
                minimum=np.minimum(minimum,actual);maximum=np.maximum(maximum,actual)
                if (t+1)%10000==0: windows.append(acc/10000);acc=0.
                if t>=110000:tail.append(actual)
        np.savez_compressed(output/f'author-{b}-{c}.npz',tail=tail,windows=windows,end=actual)
        assert short_error<1e-10
        reports.append(dict(b=b,c=c,leading_exponent=float(np.mean(windows)),windows=windows,
                            short_error=short_error,minimum=minimum.tolist(),maximum=maximum.tolist()))
        print(json.dumps(reports[-1]),flush=True)
    return dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),torch=torch.__version__,
                warmup=20000,observed=100000,runs=reports,
                note='Actual upstream bias-corrected state updates; tangent uses independently checked limiting Jacobian.')


def precision_check(output, digits, separation):
    import mpmath as mp
    mp.mp.dps=digits
    alpha=mp.mpf('0.01');eps=mp.mpf('1e-16');scale=[alpha,alpha,alpha**2]
    delta=mp.mpf(separation)
    reports=[]
    for bs,cs in [('0.52','0.15'),('0.9','0.99')]:
        b,c=mp.mpf(bs),mp.mpf(cs)
        # Original units, separately transcribed source recurrence, no Jacobian.
        def advance(z):
            x,m,s=z
            mm=b*m+(1-b)*x
            ss=c*s+(1-c)*(x-mm)**2+eps
            return [x-alpha*mm/(mp.sqrt(ss)+eps),mm,ss]
        z=[mp.mpf('.0073'),mp.mpf(0),mp.mpf(0)]
        v=[mp.mpf(1),mp.mpf('.3'),mp.mpf('.2')]
        norm=mp.sqrt(sum(w*w for w in v));v=[w/norm for w in v]
        sums=mp.mpf(0);windows=[];tail=[]
        for t in range(120000):
            shifted=[z[i]+delta*v[i]*scale[i] for i in range(3)]
            zn=advance(z);sn=advance(shifted)
            diff=[(sn[i]-zn[i])/scale[i] for i in range(3)]
            norm=mp.sqrt(sum(w*w for w in diff))
            v=[w/norm for w in diff];z=zn
            if t>=20000:
                sums+=mp.log(norm/delta)
                if (t+1)%10000==0:
                    windows.append(float(sums/10000));sums=mp.mpf(0)
                if t>=118000: tail.append([float(z[i]/scale[i]) for i in range(3)])
        np.savez_compressed(output/f'precision-{bs}-{cs}.npz',tail=tail,windows=windows)
        reports.append(dict(b=bs,c=cs,leading_exponent=float(np.mean(windows)),windows=windows,
                            terminal_original_units=[mp.nstr(w,digits) for w in z]))
        print(json.dumps(reports[-1]),flush=True)
    return dict(digits=digits,separation=separation,warmup=20000,observed=100000,
                mpmath=mp.__version__,runs=reports,
                note='Renormalized finite secants in original units, no analytic Jacobian or QR.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--author-source',type=Path);p.add_argument('--digits',type=int,default=70)
    p.add_argument('--separation',default='1e-25');args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);t=time.time()
    if args.author_source:report=author_check(args.author_source,args.output)
    else:report=precision_check(args.output,args.digits,args.separation)
    report['seconds']=time.time()-t
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__': main()
