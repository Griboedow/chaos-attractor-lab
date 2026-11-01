"""Independent arbitrary-precision RAAR reflector replay, without analytic Jacobian."""
import argparse,concurrent.futures,hashlib,json,time
from pathlib import Path
import mpmath as mp

def projection(x,center,normal,radius):
    d=[x[j]-center[j] for j in range(3)]
    dot=sum(d[j]*normal[j] for j in range(3))
    v=[d[j]-dot*normal[j] for j in range(3)]
    length=mp.sqrt(sum(a*a for a in v))
    if length<mp.mpf('1e-12'):raise RuntimeError('Projection axis encountered')
    return [center[j]+radius*v[j]/length for j in range(3)]


def step(x,params):
    beta,radius,normal,center=params
    b=projection(x,center,normal,radius)
    rb=[2*b[j]-x[j] for j in range(3)]
    a=projection(rb,[mp.mpf(0)]*3,[mp.mpf(0),mp.mpf(0),mp.mpf(1)],mp.mpf(1))
    ra=[2*a[j]-rb[j] for j in range(3)]
    return [beta*(ra[j]+x[j])/2+(1-beta)*b[j] for j in range(3)]


def norm(v):return mp.sqrt(sum(x*x for x in v))

def run(task):
    case,digits,outdir=task;mp.mp.dps=digits
    beta,radius,tilt,*center=[mp.mpf(float(x)) for x in case['params']]
    params=(beta,radius,[mp.mpf(0),-mp.sin(tilt),mp.cos(tilt)],center)
    s=[mp.mpf(float(x)) for x in [.317,-.211,.419]];begin=time.time()
    for _ in range(100000):s=step(s,params)
    delta=mp.mpf(10)**(-30 if digits==70 else -40)
    q=[[mp.mpf(1),mp.mpf(0),mp.mpf(0)],[mp.mpf(0),mp.mpf(1),mp.mpf(0)]]
    rates=[[mp.mpf(0),mp.mpf(0)] for _ in range(10)];low=list(s);high=list(s);samples=[]
    for n in range(105000):
        shadows=[step([s[j]+delta*v[j] for j in range(3)],params) for v in q]
        new=step(s,params)
        v=[[ (shadow[j]-new[j])/delta for j in range(3)] for shadow in shadows]
        r0=norm(v[0]);q0=[x/r0 for x in v[0]]
        dot=sum(q0[j]*v[1][j] for j in range(3));v1=[v[1][j]-dot*q0[j] for j in range(3)]
        r1=norm(v1);q=[q0,[x/r1 for x in v1]];s=new
        if n>=5000:
            index=(n-5000)//10000;rates[index][0]+=mp.log(r0);rates[index][1]+=mp.log(r1)
            low=[min(a,b) for a,b in zip(low,s)];high=[max(a,b) for a,b in zip(high,s)]
            if (n-5000)%100==0:samples.append([float(x) for x in s])
    rates=[[v/10000 for v in row] for row in rates]
    report=dict(id=case['id'],digits=digits,parameters=case['params'],parameter_hex=[float(x).hex() for x in case['params']],
                warm=100000,alignment=5000,observed=100000,windows=10,delta=str(delta),
                rates_decimal=[[mp.nstr(v,40) for v in row] for row in rates],rates=[[float(v) for v in row] for row in rates],
                means=[float(sum(row[j] for row in rates)/10) for j in range(2)],
                minima=[float(min(row[j] for row in rates)) for j in range(2)],bounds=[[float(v) for v in low],[float(v) for v in high]],
                samples_stride=100,samples=samples,seconds=time.time()-begin,mpmath_version=mp.__version__,
                method='Independent nested reflection formula and two finite-separation directions in actual x,y,z; no production map or analytic Jacobian imported. Parameters are exact binary64 lifts; trigonometry is recomputed at the stated precision.')
    path=Path(outdir)/f'case-{case["id"]}-{digits}.json';path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    return {k:report[k] for k in ['id','digits','means','minima','seconds']}

def main(args):
    args.out.mkdir(parents=True,exist_ok=False);cases=json.loads(args.cases.read_text())['cases']
    (args.out/'input.json').write_bytes(args.cases.read_bytes());(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    tasks=[(case,digits,str(args.out)) for case in cases for digits in [70,90]];results=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        for future in concurrent.futures.as_completed([pool.submit(run,task) for task in tasks]):
            result=future.result();results.append(result);print(json.dumps(result),flush=True)
            (args.out/'progress.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    (args.out/'report.json').write_text(json.dumps(dict(runs=results,input_sha256=hashlib.sha256(args.cases.read_bytes()).hexdigest()),indent=2),encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cases',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
