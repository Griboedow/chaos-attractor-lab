"""Independent MP replay using the four RK stages, not the polynomial/Jacobian.

Two finite-separation directions are reorthogonalized in actual (q,p,h) states.
Parameters and initial state are exact lifts of recorded binary64 values.
"""
import argparse,concurrent.futures,hashlib,json,time
from pathlib import Path
import mpmath as mp

def field(y,zeta):return [y[1],-y[0]-2*zeta*y[1]]

def step(s,zeta,rtol,safety,minimum,maximum):
    y=s[:2];h=s[2];reject=False
    for _ in range(100):
        k1=field(y,zeta)
        k2=field([y[j]+h*k1[j]/2 for j in range(2)],zeta)
        k3=field([y[j]+3*h*k2[j]/4 for j in range(2)],zeta)
        out=[y[j]+h*(2*k1[j]/9+k2[j]/3+4*k3[j]/9) for j in range(2)]
        k4=field(out,zeta)
        error=[h*(5*k1[j]/72-k2[j]/12-k3[j]/9+k4[j]/8) for j in range(2)]
        norm=mp.sqrt(sum((error[j]/(1+rtol*max(abs(y[j]),abs(out[j]))))**2 for j in range(2))/2)
        factor=maximum if norm==0 else safety*mp.power(norm,-mp.mpf(1)/3)
        if norm<1:
            factor=min(maximum,factor)
            if reject:factor=min(mp.mpf(1),factor)
            return out+[h*factor]
        h*=max(minimum,factor);reject=True
    raise RuntimeError('Rejected 100 trial steps')

def norm(v):return mp.sqrt(sum(x*x for x in v))

def run(task):
    case,digits,outdir=task;mp.mp.dps=digits
    zeta,rtol=[mp.mpf(float(x)) for x in case['params']]
    safety=mp.mpf(float(.9));minimum=mp.mpf(float(.2));maximum=mp.mpf(10)
    s=[mp.mpf(1),mp.mpf(0),mp.mpf(float(.1))];begin=time.time()
    for _ in range(100000):s=step(s,zeta,rtol,safety,minimum,maximum)
    delta=mp.mpf(10)**(-30 if digits==70 else -40)
    q=[[mp.mpf(1),mp.mpf(0),mp.mpf(0)],[mp.mpf(0),mp.mpf(1),mp.mpf(0)]]
    rates=[[mp.mpf(0),mp.mpf(0)] for _ in range(10)];low=list(s);high=list(s);samples=[]
    for n in range(105000):
        shadows=[step([s[j]+delta*v[j] for j in range(3)],zeta,rtol,safety,minimum,maximum) for v in q]
        new=step(s,zeta,rtol,safety,minimum,maximum)
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
                method='Independent four-stage Bogacki-Shampine formulas and two finite-separation directions in actual q,p,h; no production map or analytic Jacobian imported.')
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
