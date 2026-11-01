"""Complete interval-Sturm preimage counts for selected spatial-circle RAAR targets."""
import argparse,hashlib,json,time
from pathlib import Path
import mpmath as mp
from certify_raar_fixed_points import exact,add,scale,mul,sign,contains_zero,evaluate,sturm,variations,at


def one(case):
    iv=mp.iv;beta,r,theta,cx,cy,cz=map(exact,case['params']);yx,yy,yz=map(exact,case['target']);ct=iv.cos(theta);st=iv.sin(theta)
    D=list(map(iv.mpf,[1,0,1]));S=list(map(iv.mpf,[0,2]));C=list(map(iv.mpf,[1,0,-1]))
    X=add(scale(D,cx-yx),scale(C,r));Y=add(scale(D,cy-yy),scale(S,r*ct))
    f=add(scale(S,-(yx-(1-beta)*cx)),scale(C,ct*(yy-(1-beta)*cy)+st*(yz-(1-beta)*cz)))
    e=add(add(scale(mul(S,D),-(cx-yx)),scale(mul(C,D),ct*(cy-yy))),scale(mul(S,C),-r*st*st))
    polynomial=add(mul(mul(f,f),add(mul(X,X),mul(Y,Y))),scale(mul(e,e),-beta*beta));assert len(polynomial)==9;sign(polynomial[-1]);sequence=sturm(polynomial)
    sl=[sign(p[-1])*(-1)**(len(p)-1) for p in sequence];sr=[sign(p[-1]) for p in sequence];total=variations(sl)-variations(sr);assert total==len(case['root_intervals'])
    ellipse=((yx-cx)/r)**2+((yy-cy)/(r*ct))**2-1;sign(ellipse)
    lpi=iv.sqrt((cx-r-yx)**2+(cy-yy)**2);fpi=-(ct*(yy-(1-beta)*cy)+st*(yz-(1-beta)*cz));epi=-ct*(cy-yy)
    pi_values=[fpi*lpi-beta*s*epi for s in [-1,1]];assert all(sign(v) in [-1,1] for v in pi_values)
    rows=[];valid=[];previous=None
    for j,(lo,hi) in enumerate(case['root_intervals']):
        assert lo<hi and (previous is None or previous<lo);previous=hi
        lo,hi=exact(lo),exact(hi);ls=at(sequence,lo);rs=at(sequence,hi);assert variations(ls)-variations(rs)==1
        t=iv.mpf([lo.a,hi.b]);den=1+t**2;sn=2*t/den;cs=(1-t**2)/den;b=[cx+r*cs,cy+r*ct*sn,cz+r*st*sn];dx=b[0]-yx;dy=b[1]-yy;length=iv.sqrt(dx**2+dy**2);sign(length)
        ff=evaluate(f,t)/den;ee=evaluate(e,t)/den**2;sigma=sign(ff)*sign(ee);assert contains_zero(ff*length-beta*sigma*ee)
        a=[sigma*dx/length,sigma*dy/length,iv.mpf(0)];x=[(y-beta*aa-(1-2*beta)*bb)/beta for y,aa,bb in zip([yx,yy,yz],a,b)]
        radialA=1+sigma*length/beta;radialB=(x[0]-cx)*cs+(x[1]-cy)*ct*sn+(x[2]-cz)*st*sn;sa=sign(radialA);sb=sign(radialB);good=sa>0 and sb>0
        row=dict(root=j,half_angle_interval=[str(lo),str(hi)],sturm_left_signs=ls,sturm_right_signs=rs,sigma=sigma,radial_A=str(radialA),radial_B=str(radialB),valid_preimage=good,preimage_box=[str(v) for v in x],point_B_box=[str(v) for v in b]);rows.append(row)
        if good:valid.append(row)
    assert len(valid)>=2
    # Distinct image roots must also have disjoint reported coordinate boxes.
    for j,a in enumerate(valid):
        for b in valid[j+1:]:
            def parse(s):return iv.mpf([v.strip() for v in s.strip('[]').split(',')])
            aa=list(map(parse,a['preimage_box']));bb=list(map(parse,b['preimage_box']));assert any(x.b<y.a or y.b<x.a for x,y in zip(aa,bb))
    return dict(label=case['label'],parameters=case['params'],target=case['target'],degree=8,sturm_degrees=[len(p)-1 for p in sequence],sturm_negative_infinity_signs=sl,sturm_positive_infinity_signs=sr,total_real_roots=total,number_of_preimages=len(valid),projected_circle_target_exclusion=str(ellipse),pi_chart_exclusions=[str(v) for v in pi_values],roots=rows)


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);mp.iv.dps=100;begin=time.time();cases=json.loads(args.input.read_text())['cases'];rows=[]
    for name in ['certify_raar_inverse.py','certify_raar_fixed_points.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    (args.out/'input.json').write_bytes(args.input.read_bytes())
    for case in cases:
        row=one(case);rows.append(row);(args.out/(row['label']+'.json')).write_text(json.dumps(row,indent=2),encoding='utf-8');print(json.dumps({k:row[k] for k in ['label','total_real_roots','number_of_preimages']}),flush=True)
    report=dict(status='passed',interval_digits=100,mpmath=mp.__version__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),cases=rows,seconds=time.time()-begin,scope='Complete preimage count for each specified exact binary64 target in the regular domain, with the two omitted angular/axis cases explicitly excluded. This proves noninjectivity of the full map. It does not prove that multiple preimages lie on the asymptotic attracting set, or inequivalence of that set to all prior attractors.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
