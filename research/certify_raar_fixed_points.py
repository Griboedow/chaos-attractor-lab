"""Interval Sturm certificate for fixed points of selected spatial-circle RAAR maps.

The mathematical parameters are exact binary64 lifts; sin/cos are interval
enclosures at 80 decimal digits. Float polynomial roots only propose rational
isolating intervals. Completeness and signs are checked independently here.
"""
import argparse,hashlib,json
from pathlib import Path
import mpmath as mp

iv=mp.iv


def exact(x):
    n,d=float(x).as_integer_ratio();return iv.mpf(n)/d


def add(a,b):
    c=[iv.mpf(0) for _ in range(max(len(a),len(b)))]
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return c


def scale(a,s):return [x*s for x in a]


def mul(a,b):
    c=[iv.mpf(0) for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def sign(x):
    if x.a>0:return 1
    if x.b<0:return -1
    raise ArithmeticError('Sign not certified: '+str(x))


def contains_zero(x):return x.a<=0 and x.b>=0


def evaluate(p,x):
    y=iv.mpf(0)
    for c in p[::-1]:y=y*x+c
    return y


def remainder(a,b):
    sign(b[-1]);a=list(a)
    while len(a)>=len(b):
        q=a[-1]/b[-1];offset=len(a)-len(b)
        for j in range(len(b)-1):a[offset+j]-=q*b[j]
        # The exact leading coefficient cancels identically in polynomial division.
        a.pop()
    return a


def sturm(p):
    sequence=[p,[i*p[i] for i in range(1,len(p))]]
    while len(sequence[-1])>1:
        row=scale(remainder(sequence[-2],sequence[-1]),iv.mpf(-1));sign(row[-1]);sequence.append(row)
    sign(sequence[-1][0]);return sequence


def variations(values):return sum(a!=b for a,b in zip(values,values[1:]))


def at(sequence,x):return [sign(evaluate(p,x)) for p in sequence]


def interval_string(x):return str(x)


def one(case):
    beta,r,theta,cx,cy,cz=map(exact,case['params']);ct=iv.cos(theta);st=iv.sin(theta)
    D=list(map(iv.mpf,[1,0,1]));S=list(map(iv.mpf,[0,2]));C=list(map(iv.mpf,[1,0,-1]))
    X=add(scale(D,cx),scale(C,r));Y=add(scale(D,cy),scale(S,r*ct))
    fn=add(scale(S,-cx),scale(C,cy*ct+cz*st))
    en=add(add(scale(mul(S,D),cx),scale(mul(C,D),-ct*cy)),scale(mul(S,C),r*st*st))
    L2=add(mul(X,X),mul(Y,Y));polynomial=add(mul(mul(fn,fn),L2),scale(mul(en,en),iv.mpf(-1)))
    assert len(polynomial)==9;sign(polynomial[-1]);sequence=sturm(polynomial)
    left=[sign(p[-1])*(-1)**(len(p)-1) for p in sequence];right=[sign(p[-1]) for p in sequence]
    total=variations(left)-variations(right)
    assert total==len(case['root_intervals'])
    # B never crosses the z-axis, so each stationary A point has one of two signs.
    axis_test=(cx/r)**2+(cy/(r*ct))**2-1;sign(axis_test)
    # gamma=pi is outside the tan(gamma/2) chart; check both stationary branches.
    bpi=[cx-r,cy,cz];Lpi=iv.sqrt(bpi[0]**2+bpi[1]**2)
    pi_stationarity=[-(cy*ct+cz*st)+s*ct*cy/Lpi for s in [-1,1]]
    assert all(sign(x) in [-1,1] for x in pi_stationarity)
    rows=[];fixed=[]
    previous=None
    for a,b in case['root_intervals']:
        if previous is not None:assert a>previous
        previous=b;aa=exact(a);bb=exact(b);sa=at(sequence,aa);sb=at(sequence,bb)
        assert variations(sa)-variations(sb)==1
        u=iv.mpf([aa.a,bb.b]);den=1+u*u;sn=2*u/den;cs=(1-u*u)/den
        pointB=[cx+r*cs,cy+r*ct*sn,cz+r*st*sn];L=iv.sqrt(pointB[0]**2+pointB[1]**2);sign(L)
        f=-cx*sn+(cy*ct+cz*st)*cs;e=cx*sn-ct*cy*cs+r*st*st*sn*cs
        branch=[]
        for sigma in [-1,1]:
            residual=f+sigma*e/L
            if not contains_zero(residual):
                branch.append(dict(sigma=sigma,rejected='stationarity',residual=interval_string(residual)));continue
            A=[sigma*pointB[0]/L,sigma*pointB[1]/L,iv.mpf(0)]
            k=beta/(1-beta);x=[pointB[j]+k*(A[j]-pointB[j]) for j in range(3)]
            radialA=sum((2*pointB[j]-x[j])*A[j] for j in range(3))
            unitB=[cs,ct*sn,st*sn];center=[cx,cy,cz]
            radialB=sum((x[j]-center[j])*unitB[j] for j in range(3))
            signA=sign(radialA);signB=sign(radialB);valid=signA>0 and signB>0
            row=dict(sigma=sigma,stationarity_interval=interval_string(residual),radial_A=interval_string(radialA),radial_B=interval_string(radialB),valid=valid)
            if valid:
                # J has an exact eigenvalue beta. Bound the remaining quadratic.
                ta=[-A[1],A[0],iv.mpf(0)];tb=[-sn,ct*cs,st*cs];cosine=sum(ta[j]*tb[j] for j in range(3))
                al=1/radialA;be=r/radialB
                trace=2*beta-beta*al+(1-2*beta)*be+2*beta*al*be*cosine**2
                determinant=beta*((1-al)*(beta+(1-2*beta)*be)+al*be*cosine**2)
                q1=1-trace+determinant;qm1=1+trace+determinant
                unstable=(determinant.a>1 or determinant.b < -1 or q1.b<0 or qm1.b<0)
                row.update(state_box=[interval_string(v) for v in x],quadratic_trace=interval_string(trace),quadratic_determinant=interval_string(determinant),quadratic_at_plus_one=interval_string(q1),quadratic_at_minus_one=interval_string(qm1),unstable_certified=bool(unstable))
                fixed.append(row)
            branch.append(row)
        # Squaring introduces two branches; one and only one must remain stationary.
        assert sum('valid' in row for row in branch)==1
        rows.append(dict(interval=[a,b],left_signs=sa,right_signs=sb,branches=branch))
    # Infeasibility of the two circle constraints: A's intersections with B's plane.
    plane_y=cy-cz*ct/st
    if plane_y.a>1 or plane_y.b < -1:
        intersection=dict(disjoint=True,reason='B plane misses the closed A disk',plane_y=interval_string(plane_y))
    else:
        root=iv.sqrt(1-plane_y**2)
        residuals=[(s*root-cx)**2+(plane_y-cy)**2+cz**2-r**2 for s in [-1,1]]
        assert all(sign(v) in [-1,1] for v in residuals)
        intersection=dict(disjoint=True,reason='Both plane-line candidates fail the B circle equation',squared_radius_residuals=[interval_string(v) for v in residuals])
    return dict(id=case['id'],parameters=case['params'],parameter_hex=[float(x).hex() for x in case['params']],
        polynomial_coefficients=[interval_string(v) for v in polynomial],sturm_degrees=[len(p)-1 for p in sequence],
        sturm_minus_infinity_signs=left,sturm_plus_infinity_signs=right,total_real_roots=total,
        axis_exclusion=interval_string(axis_test),pi_stationarity_exclusion=[interval_string(v) for v in pi_stationarity],
        roots=rows,number_of_fixed_points=len(fixed),all_fixed_points_unstable=all(row['unstable_certified'] for row in fixed),intersection=intersection)


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);iv.dps=80
    cases=json.loads(args.input.read_text())['cases'];rows=[]
    (args.out/'input.json').write_bytes(args.input.read_bytes());(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    for case in cases:
        row=one(case);rows.append(row);print(json.dumps(dict(id=row['id'],real_roots=row['total_real_roots'],fixed_points=row['number_of_fixed_points'],all_fixed_points_unstable=row['all_fixed_points_unstable'])),flush=True)
        (args.out/f'case-{case["id"]}.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
    (args.out/'report.json').write_text(json.dumps(dict(status='passed',precision=80,mpmath_version=mp.__version__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),cases=rows,scope='Complete fixed-point and circle-disjointness classification at the six exact parameter settings. Not a chaos, basin, priority or structural-stability certificate.'),indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
