"""Independent high-precision containment tests for the soft-map extension."""
import argparse,hashlib,json,math,sys,time
from fractions import Fraction
from pathlib import Path
import numpy as np
from soft_interval import exp_point,evaluate,remainder
from search_soft_constellations import step


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    if args.mpmath_path:sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps=100;begin=time.time();rng=np.random.default_rng(2026110122)
    bound=Fraction(8,7)*Fraction(1,8)**13/math.factorial(13)
    assert bound<Fraction(1,2**71)
    def contained(v,iv):return mp.mpf(float(iv[0]))<=v<=mp.mpf(float(iv[1]))
    arguments=np.r_[0.,.125,-.125,700.,-700.,np.nextafter(.125,np.inf),np.nextafter(-.125,-np.inf),rng.uniform(-700,700,1024)]
    largest_relative=0.
    for x in arguments:
        interval=exp_point(float(x));actual=mp.exp(mp.mpf(float(x)));assert contained(actual,interval),(x,interval,actual)
        largest_relative=max(largest_relative,float((mp.mpf(interval[1])-mp.mpf(interval[0]))/actual))
    checks=[]
    for ident in [29,156]:
        source=args.source/f'case-{ident}.npz'
        with np.load(source) as f:beta,tau=f['params'];p=f['points'];states=f['states']
        b=mp.mpf(float(beta));temperature=mp.mpf(float(tau));sites=[[mp.mpf(float(a)) for a in q] for q in p]
        def independent(x):
            logits=[(2*sum(q[j]*x[j] for j in range(3))-sum(a*a for a in q))/temperature for q in sites];offset=max(logits)
            weights=[mp.exp(a-offset) for a in logits];total=sum(weights);weights=[a/total for a in weights]
            mean=[sum(weights[k]*sites[k][j] for k in range(len(p))) for j in range(3)]
            rb=[2*mean[j]-x[j] for j in range(3)];r=mp.sqrt(rb[0]**2+rb[1]**2);a=[rb[0]/r,rb[1]/r,mp.mpf(0)]
            return [b/2*(2*a[j]-rb[j]+x[j])+(1-b)*mean[j] for j in range(3)]
        count=0
        for k in range(16):
            mid=states[int(rng.integers(len(states)))];width=2.**(-8-k%5);box=np.c_[mid-width/2,mid+width/2]
            fb,jb,axis,flag=evaluate(box,beta,tau,p);assert not flag
            center,J,*_=step(mid,beta,p,tau);radii,errors,axis2,flag=remainder(box,mid,center,J,beta,tau,p,2.**-35);assert not flag
            samples=[mid]+[box[:,0]+np.array([(c>>j)&1 for j in range(3)])*(box[:,1]-box[:,0]) for c in range(8)]
            samples+=list(rng.uniform(box[:,0],box[:,1],(3,3)))
            for sample in samples:
                x=[mp.mpf(float(a)) for a in sample];value=independent(x)
                for i in range(3):
                    assert contained(value[i],fb[i]),(ident,k,'value',i)
                    linear=mp.mpf(float(center[i]))+sum(mp.mpf(float(J[i,j]))*(x[j]-mp.mpf(float(mid[j]))) for j in range(3))
                    assert abs(value[i]-linear)<=mp.mpf(float(errors[i])),(ident,k,'remainder',i)
                # Differentiate the independent nested-reflector expression.
                for j in range(3):
                    def column(t):return mp.matrix(independent([t if h==j else x[h] for h in range(3)]))
                    derivative=mp.diff(column,x[j])
                    for i in range(3):assert contained(derivative[i],jb[i,j]),(ident,k,'Jacobian',i,j)
                count+=1
        checks.append(dict(id=ident,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),boxes=16,points=count,derivatives=9*count))
    files={}
    for name in [Path(__file__).name,'soft_interval.py','raar_interval_scalar.py']:
        data=Path(__file__).with_name(name).read_bytes();(args.out/name).write_bytes(data);files[name]=hashlib.sha256(data).hexdigest()
    result=dict(status='passed',precision_digits=100,taylor_remainder_bound='(8/7)*(1/8)^13/13! < 2^-71, checked by exact rational arithmetic',exponential_points=len(arguments),maximum_exponential_relative_width=largest_relative,cases=checks,sources=files,seconds=time.time()-begin,scope='Containment regression tests and exact Taylor remainder arithmetic. The tests support implementation validation, not exhaustive correctness or a trapping theorem.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--mpmath-path',type=Path);main(p.parse_args())
