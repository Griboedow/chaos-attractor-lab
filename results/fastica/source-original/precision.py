"""Independent mpmath tangent growth for exact binary64 input lifts."""
from pathlib import Path
import argparse, json, time, hashlib
import mpmath as mp

def calculate(data,digits,warmup,steps,window):
    mp.mp.dps=digits
    X=[[mp.mpf(float(v)) for v in row] for row in data['X']]
    w=[mp.mpf(float(v)) for v in data['start']]
    dot=lambda a,b:mp.fsum(x*y for x,y in zip(a,b))
    norm=lambda a:mp.sqrt(dot(a,a))
    z=norm(w);w=[x/z for x in w]
    v=[mp.mpf(1)]*4;c=dot(w,v);v=[a-c*b for a,b in zip(v,w)]
    z=norm(v);v=[x/z for x in v]
    alpha=mp.mpf(float(data['alpha']));kind=data['kind'];N=len(X)
    sums=[mp.mpf(0)]*(steps//window);minimum=mp.inf;residual=mp.mpf(0)
    for it in range(warmup+steps):
        h=[mp.mpf(0)]*4;jv=[mp.mpf(0)]*4;p=mp.mpf(0);b=mp.mpf(0)
        for x in X:
            u=dot(x,w);z=dot(x,v)
            if kind==0:g=mp.tanh(alpha*u);gp=alpha*(1-g*g);gpp=-2*alpha*gp*g
            elif kind==1:e=mp.exp(-u*u/2);g=u*e;gp=(1-u*u)*e;gpp=(u*u*u-3*u)*e
            else:g=u**3;gp=3*u*u;gpp=6*u
            p+=gp;b+=gpp*z
            for k in range(4):h[k]+=x[k]*g;jv[k]+=x[k]*gp*z
        h=[(a-p*c)/N for a,c in zip(h,w)]
        jv=[(a-p*c-b*d)/N for a,c,d in zip(jv,v,w)]
        hn=norm(h);minimum=min(minimum,hn);w=[a/hn for a in h]
        c=dot(w,jv);v=[(a-c*b)/hn for a,b in zip(jv,w)]
        gain=norm(v);v=[a/gain for a in v]
        residual=max(residual,abs(dot(w,w)-1))
        if it>=warmup:sums[(it-warmup)//window]+=mp.log(gain)/window
    return {'case':data['case'],'digits':digits,'warmup':warmup,'observed':steps,'window':window,'rates':[str(v) for v in sums],'minimum_update_norm':str(minimum),'sphere_residual':str(residual),'input_semantics':'Exact lifts of parsed binary64 X and initial state; normalization and all subsequent dynamics in mpmath.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('D:/Codex/unexpected-attractors-20261101/fastica'));p.add_argument('--digits',type=int,default=70);p.add_argument('--cases',type=int,nargs='+',default=[15,204,246,351]);a=p.parse_args();rows=[]
    for case in a.cases:
        path=a.root/f'follow-01/case-{case}-input.json';t=time.time();row=calculate(json.loads(path.read_text()),a.digits,1000,8000,2000);row.update(wall_seconds=time.time()-t,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),mpmath_version=mp.__version__);rows.append(row);print(json.dumps(row),flush=True);(a.root/f'precision-{a.digits}.json').write_text(json.dumps(rows,indent=2)+'\n')
