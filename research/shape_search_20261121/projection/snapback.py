"""Validate two inverse branches near a periodic full-dimensional repeller."""
import argparse
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from certify_cycles import iterate,verify,eye,matvec,matmul,det,exactfloat

iv=mp.iv
enclosure_model=iterate


def block(center,radius):
    return [iv.mpf([(v-radius).a,(v+radius).b]) for v in center]


def infnorm(matrix):
    return max(sum(abs(a).b for a in row).b for row in matrix)


def frac_iv(v):
    return iv.mpf(v.numerator)/v.denominator


def exact_basis(raw):
    b=[[Fraction(float(v)) for v in row] for row in raw]
    determinant=det(b)
    # Cyclic minor ordering already includes the cofactor sign.
    inv=[[(b[(j+1)%3][(i+1)%3]*b[(j+2)%3][(i+2)%3]-b[(j+1)%3][(i+2)%3]*b[(j+2)%3][(i+1)%3])/determinant for j in range(3)] for i in range(3)]
    assert matmul(inv,b)==[[int(i==j) for j in range(3)] for i in range(3)]
    return [[frac_iv(v) for v in row] for row in b],[[frac_iv(v) for v in row] for row in inv]


def decimal(v): return mp.nstr(v,145)


def floats(x): return [mp.mpf(v) for v in x]


def refined(seed,target,p,rec,periodic=False):
    x=floats(seed)
    for _ in range(20):
        y,j,_=iterate(x,p,rec,mp.mp)
        dest=x if periodic else target
        f=mp.matrix([v-t for v,t in zip(y,dest)])
        m=mp.matrix(j)
        if periodic:m-=mp.eye(3)
        if mp.norm(f,mp.inf)<mp.mpf('1e-125'): break
        s=mp.lu_solve(m,f)
        x=[x[i]-s[i] for i in range(3)]
    else: raise ArithmeticError('Newton refinement failed')
    inv=m**-1
    return dict(period=p,center=[decimal(v) for v in x],preconditioner=[[decimal(inv[i,j]) for j in range(3)] for i in range(3)])


def residual_matrix(cm,j):
    prod=matmul(cm,j)
    return [[(i==k)-prod[i][k] for k in range(3)] for i in range(3)]


def enclosed_iterate(space,p,rec):
    """Accumulate one-step Jacobians and recenter interval states after each step."""
    matrix=eye(iv)
    for _ in range(p):
        middle=[v.mid if v.mid.a>=v.a and v.mid.b<=v.b else v.a for v in space]
        value,_,_=enclosure_model(middle,1,rec,iv)
        _,j,_=enclosure_model(space,1,rec,iv)
        delta=matvec(j,[v-m for v,m in zip(space,middle)])
        space=[v+d for v,d in zip(value,delta)]
        matrix=matmul(j,matrix)
    return space,matrix


def inverse_block(cert,target_box,rec):
    x=[iv.mpf(v) for v in cert['center']]
    radius=iv.mpf(cert['radius'])
    domain=block(x,radius)
    value,_,_=enclosure_model(x,cert['period'],rec,iv)
    _,j=enclosed_iterate(domain,cert['period'],rec)
    c=[[iv.mpf(v) for v in row] for row in cert['preconditioner']]
    rem=residual_matrix(c,j)
    q=infnorm(rem)
    correction=matvec(c,[v-t for v,t in zip(value,target_box)])
    estimate=max(abs(v).b for v in correction)+q*radius
    if not q<1 or not estimate<radius: raise ArithmeticError('Parameterized inverse block failed')
    return domain,float((estimate/radius).b)


def inverse_chart_block(cert,target_center,target_radius,b,bi,rec,derivative_bounds=None):
    x=[iv.mpf(v) for v in cert['center']]
    radius=iv.mpf(cert['radius'])
    delta=matvec(b,[iv.mpf([-radius.b,radius.b])]*3)
    domain=[v+d for v,d in zip(x,delta)]
    value,_,_=enclosure_model(x,cert['period'],rec,iv)
    _,j=enclosed_iterate(domain,cert['period'],rec)
    c=[[iv.mpf(v) for v in row] for row in cert['preconditioner']]
    q=infnorm(residual_matrix(c,matmul(bi,matmul(j,b))))
    correction=matvec(c,matvec(bi,[v-t for v,t in zip(value,target_center)]))
    estimate=max(abs(v).b for v in correction)+infnorm(c)*target_radius+q*radius
    if not q<1 or not estimate<radius:raise ArithmeticError(f'Adapted inverse block failed: q={float(q.b):g}, box fraction={float((estimate/radius).b):g}')
    if derivative_bounds is not None:
        derivative_bounds.append(infnorm(c)/(1-q))
    return x,radius,float((estimate/radius).b)


def local_chart(base,raw_basis,rec,existing=None,power=1):
    b,bi=exact_basis(raw_basis)
    root=[iv.mpf(v) for v in base['center']]
    root_box=block(root,iv.mpf(base['root_radius']))
    p=base['period']*power
    if existing:
        pre=existing['preconditioner']
    else:
        _,m,_=iterate(floats(base['center']),p,rec,mp.mp)
        bm=mp.matrix(raw_basis)
        cm=(bm**-1*mp.matrix(m)*bm)**-1
        pre=[[decimal(cm[i,j]) for j in range(3)] for i in range(3)]
    c=[[iv.mpf(v) for v in row] for row in pre]
    cn=infnorm(c)
    choices=[existing['radius']] if existing else [f'1e-{j}' for j in range(3,35)]
    for radius in choices:
        rr=iv.mpf(radius)
        spread=matvec(b,[iv.mpf([-rr.b,rr.b])]*3)
        space=[x+y for x,y in zip(root_box,spread)]
        try:
            _,m=enclosed_iterate(space,p,rec)
            q=infnorm(residual_matrix(c,matmul(bi,matmul(m,b))))
            if not q+cn<1: continue
            k=cn/(1-q)
            return dict(radius=radius,preconditioner=pre,q_display=float(q.b),inverse_lipschitz_display=float(k.b)),k,b,bi,root_box
        except (ArithmeticError,ZeroDivisionError,ValueError): continue
    raise ArithmeticError('Uniform local inverse contraction failed')


def close_branches(doc,rec):
    base=doc['fixed'];power=doc.get('local_power',1);p=base['period']*power
    vr=verify(base,rec,base['root_radius'],enclosure_model)
    if vr['unstable_dimension']!=3:raise ArithmeticError('Return is not a repeller')
    chart,k,b,bi,root_box=local_chart(base,doc['basis'],rec,doc['local'],power)
    root=[iv.mpf(v) for v in base['center']]
    target_radius=iv.mpf(doc['target_radius'])
    target=block(root,target_radius)
    adapted=doc.get('adapted_chain',False)
    small=iv.mpf(doc['small_chart_radius'])
    # The common target chart is centered on the exact periodic point.
    needed=small+infnorm(bi)*iv.mpf(base['root_radius']) if adapted else infnorm(b)*small+iv.mpf(base['root_radius'])
    if not needed<target_radius:
        raise ArithmeticError('Common target exceeds the parameterized inverse domain')
    ratios=[];total=0
    derivative_bounds=[]
    target_center=root
    for cert in doc['blocks']:
        if adapted:
            target_center,target_radius,ratio=inverse_chart_block(cert,target_center,target_radius,b,bi,rec,derivative_bounds)
        else:
            target,ratio=inverse_block(cert,target,rec)
        ratios.append(ratio);total+=cert['period']
    if total%p:raise ArithmeticError('Compound return is not an iterate of the local return')
    if adapted:
        diff=[v+iv.mpf([-target_radius.b,target_radius.b]) for v in matvec(bi,[u-v for u,v in zip(target_center,root_box)])]
    else:
        diff=matvec(bi,[u-v for u,v in zip(target,root_box)])
    rr=iv.mpf(chart['radius'])
    if not small<rr:raise ArithmeticError('Common target exceeds the local chart')
    if not max(abs(v).b for v in diff)<rr:raise ArithmeticError('Return preimage leaves the local chart')
    n=0
    if doc.get('direct',False):
        if not max(abs(v).b for v in diff)<small:raise ArithmeticError('Second inverse branch exceeds its common target')
        if not any(abs(v).a>k**(total//p)*small for v in diff):raise ArithmeticError('Return branches overlap')
    else:
        if not any(abs(v).a>small for v in diff):raise ArithmeticError('Return branches overlap')
        while not k**n*rr<small:
            n+=1
            if n>50000:raise ArithmeticError('Unresolved local inverse contraction')
    steps=total+n*p
    contraction={}
    if adapted and doc.get('direct',False):
        excursion=iv.mpf(1)
        for bound in derivative_bounds:excursion*=bound
        local=k**(total//p)
        if not excursion<1 or not local<1:raise ArithmeticError('Both branches must contract in chart coordinates')
        contraction=dict(excursion_lipschitz_upper_display=float(excursion.b),local_branch_lipschitz_upper_display=float(local.b))
    return dict(**contraction,fixed_period=base['period'],local_return_period=p,excursion_steps=total,local_inverse_repetitions=n,uniform_return_steps=steps,entropy_lower_numerator='log(2)',entropy_lower_denominator=steps,local_contraction_upper_display=float(k.b),root_excluded_from_other_branch=True,largest_inverse_box_fraction_display=max(ratios),preimage_chart_distance_upper_display=float(max(abs(v).b for v in diff)),periodic_root_check=vr)


def generate(proposal,rec,target_power,adapted=False,chart_override=None,direct=False,local_power=1):
    base=refined(proposal['base'],None,proposal['period'],rec,True)
    base['root_radius']='1e-115'
    vr=verify(base,rec,base['root_radius'],enclosure_model)
    assert vr['unstable_dimension']==3
    local,k,b,bi,root_box=local_chart(base,proposal['basis'],rec,power=local_power)
    if chart_override:
        local['radius']=chart_override
        local,k,b,bi,root_box=local_chart(base,proposal['basis'],rec,local,local_power)
    print('local',local['radius'],float(k.b),flush=True)
    root=[iv.mpf(v) for v in base['center']]
    target_radius=f'1e-{target_power}'
    target=block(root,iv.mpf(target_radius));width=mp.mpf(target_radius)
    current_center=root;current_radius=iv.mpf(target_radius)
    bm=mp.matrix(proposal['basis'])
    centers=floats(base['center'])
    blocks=[]
    small=decimal(min(mp.mpf(target_radius)*mp.mpf('.5'),mp.mpf(local['radius'])*mp.mpf('.5'))) if direct else f'1e-{target_power+2}'
    local_period=base['period']*local_power
    itinerary=list(proposal['blocks'])
    if local_power>1:
        # Once close to the repeller, invert multiple primitive periods jointly.
        # Grouped derivatives avoid loss from independent axis-aligned bounds.
        length=0
        for i,step in enumerate(itinerary):
            length+=step['period']
            if length%local_period==0 and max(abs(float(v)-float(w)) for v,w in zip(step['seed'],base['center']))<.0003:
                itinerary=itinerary[:i+1]
                break
    extension=[dict(seed=base['center'],period=local_period) for _ in range(320)]
    for index,row in enumerate(itinerary+extension):
        new=refined(row['seed'],centers,row['period'],rec)
        c=mp.matrix([[mp.mpf(v) for v in line] for line in new['preconditioner']])
        if adapted:
            c=bm**-1*c*bm
            new['preconditioner']=[[decimal(c[i,j]) for j in range(3)] for i in range(3)]
        # A radius proposal only. The interval Krawczyk test below proves it.
        amp=max(sum(abs(c[i,j]) for j in range(3)) for i in range(3))
        radius=max(mp.mpf('1e-112'),width*amp*mp.mpf('1.03'))
        for factor in ('1','1.16','1.6','2.5'):
            trial_radius=radius*mp.mpf(factor)
            new['radius']=decimal(trial_radius)
            try:
                if adapted:
                    nc,nr,_=inverse_chart_block(new,current_center,current_radius,b,bi,rec)
                else:
                    nt,_=inverse_block(new,target,rec)
            except (ArithmeticError,ValueError,ZeroDivisionError) as exc:
                last_error=exc
                continue
            radius=trial_radius
            if adapted:current_center,current_radius=nc,nr
            else:target=nt
            break
        else:raise ArithmeticError(f'Inverse block {index+1} failed') from last_error
        blocks.append(new);centers=floats(new['center']);width=radius
        if adapted:
            dist=[v+iv.mpf([-current_radius.b,current_radius.b]) for v in matvec(bi,[a-bb for a,bb in zip(current_center,root_box)])]
        else:
            dist=matvec(bi,[a-bb for a,bb in zip(target,root_box)])
        if index%10==0:
            print('block',index+1,'radius',float(width),'distance',float(max(abs(v).b for v in dist)),flush=True)
        lo=k**(sum(c['period'] for c in blocks)//local_period)*iv.mpf(small) if direct else iv.mpf(small)
        hi=iv.mpf(small) if direct else iv.mpf(local['radius'])
        if sum(c['period'] for c in blocks)%local_period==0 and max(abs(v).b for v in dist)<hi and any(abs(v).a>lo for v in dist):
            doc=dict(case=rec['case'],precision_decimal=mp.mp.dps,basis=proposal['basis'],fixed=base,local=local,local_power=local_power,target_radius=target_radius,small_chart_radius=small,adapted_chain=adapted,direct=direct,blocks=blocks)
            doc['proof']=close_branches(doc,rec)
            return doc
    raise ArithmeticError('No nontrivial preimage entered the local chart')


def main(args):
    global enclosure_model
    if args.engine=='analytic':
        from analytic_bounds import iterate as enclosure_model
    mp.mp.dps=mp.iv.dps=160
    rec=json.loads(args.record.read_text())
    if args.verify:
        doc=json.loads(args.out.read_text())
        assert doc['case']==rec['case']
        check=close_branches(doc,rec)
        assert check['uniform_return_steps']<=doc['proof']['uniform_return_steps']
        print('verified',rec['case'],check,flush=True)
    else:
        prop=json.loads(args.proposal.read_text())
        doc=generate(prop,rec,args.target_power,args.adapted_chain,args.chart_radius,args.direct,args.local_power)
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf8')
        print('certified',rec['case'],doc['proof'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--record',type=Path,required=True)
    p.add_argument('--proposal',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    p.add_argument('--target-power',type=int,default=65)
    p.add_argument('--engine',choices=('jet','analytic'),default='jet')
    p.add_argument('--adapted-chain',action='store_true')
    p.add_argument('--chart-radius')
    p.add_argument('--local-power',type=int,default=1)
    p.add_argument('--direct',action='store_true')
    main(p.parse_args())
