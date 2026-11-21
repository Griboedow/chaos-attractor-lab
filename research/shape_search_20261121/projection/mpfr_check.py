"""Independent multiprecision two-shadow expansion with periodic orthogonalization."""
import argparse
import json
import os
from pathlib import Path
import gmpy2 as g


def add(x,y):
    return tuple(a+b for a,b in zip(x,y))


def norm(x):
    return g.sqrt(sum(v*v for v in x))


def build_map(params):
    beta=g.mpfr(params['beta']); inv=g.mpfr(params['inv_tau'])
    tube=g.mpfr(params['tube_radius'])
    points=[tuple(g.mpfr(v) for v in b) for b in params['sites']]
    linked=params['second_ring']
    temperature=g.mpfr(params.get('mix_temperature',0))
    def project(v,second):
        x,y,z=v
        if second: x,y,z=x-1,z,y
        rho=g.sqrt(x*x+y*y)
        ca=x/rho; sa=y/rho
        length=g.sqrt((rho-1)**2+z*z)
        radial=1+tube*(rho-1)/length
        px,py,pz=radial*ca,radial*sa,tube*z/length
        p=(px+1,pz,py) if second else (px,py,pz)
        return p,abs(length-tube)
    def step(x):
        logs=[-inv*sum((a-b)**2 for a,b in zip(x,p)) for p in points]
        peak=max(logs)
        weights=[g.exp(w-peak) for w in logs]
        den=sum(weights)
        m=tuple(sum(w*p[j] for w,p in zip(weights,points))/den for j in range(3))
        v=tuple(2*a-b for a,b in zip(m,x))
        p,d=project(v,False)
        if linked:
            alt,dist=project(v,True)
            if temperature:
                # Subtract the larger log weight to avoid large exponentials.
                l1=-d*d/temperature; l2=-dist*dist/temperature
                mx=max(l1,l2)
                w1=g.exp(l1-mx); w2=g.exp(l2-mx)
                p=tuple((w1*a+w2*b)/(w1+w2) for a,b in zip(p,alt))
            elif dist<d:
                p=alt
        return tuple(beta*(2*a-b+c)/2+(1-beta)*e for a,b,c,e in zip(p,v,x,m))
    return step


def run(case,bits,burn,window,windows,stride):
    with g.context(precision=bits):
        # Displacement shrinks at higher precision. All normals use MPFR.
        shift=bits//3
        half=g.mpfr(2)**-shift
        step=build_map(case['parameters'])
        x=tuple(g.mpfr(a) for a in case['initial'])
        q1=(half,0,0); q2=(0,half,0)
        rows=[[g.mpfr(0),g.mpfr(0)] for _ in range(windows)]
        for begin in range(0,burn+windows*window,stride):
            x1=add(x,q1);x2=add(x,q2)
            for _ in range(stride):
                x,x1,x2=step(x),step(x1),step(x2)
            a=tuple(v-u for v,u in zip(x1,x))
            b=tuple(v-u for v,u in zip(x2,x))
            n1=norm(a)
            unit=tuple(v/n1 for v in a)
            scalar=sum(v*w for v,w in zip(unit,b))
            orth=tuple(v-scalar*w for v,w in zip(b,unit))
            n2=norm(orth)
            q1=tuple(half*v for v in unit)
            q2=tuple(half*v/n2 for v in orth)
            if begin>=burn:
                row=rows[(begin-burn)//window]
                row[0]+=g.log(n1/half)/window
                row[1]+=g.log(n2/half)/window
        return dict(bits=bits,shift_power_of_two=-shift,burn=burn,window_length=window,windows=[[float(a),float(b)] for a,b in rows],qr_stride=stride,final_state=[str(a) for a in x])


def main(args):
    assert args.burn%args.stride==args.window%args.stride==0
    files=[args.source/f'case-{i}.json' for i in args.cases]
    args.out.mkdir(parents=True,exist_ok=True)
    for file in files:
        case=json.loads(file.read_text())
        if not case['passes']: continue
        result=args.out/f'mpfr-{case["case"]}-{args.bits}.json'
        if result.exists(): continue
        run_data=run(case,args.bits,args.burn,args.window,args.windows,args.stride)
        run_data.update(case=case['case'],gmpy2=g.version(),mpfr=g.mpfr_version(),method='Three independently advanced MPFR trajectories. Renormalized displacement vectors and orthogonalization every stride steps; no analytic Jacobian.')
        tmp=result.with_suffix('.tmp.json')
        tmp.write_text(json.dumps(run_data,indent=2)+'\n',encoding='utf8')
        os.replace(tmp,result)
        print(case['case'],args.bits,run_data['windows'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cases',type=int,nargs='+',required=True)
    p.add_argument('--bits',type=int,default=192)
    p.add_argument('--burn',type=int,default=4096)
    p.add_argument('--window',type=int,default=4096)
    p.add_argument('--windows',type=int,default=4)
    p.add_argument('--stride',type=int,default=4)
    main(p.parse_args())
