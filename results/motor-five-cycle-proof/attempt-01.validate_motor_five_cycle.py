"""Interval sixth-period proof for the weighted five-cluster attracting cycle."""
from __future__ import annotations
import argparse
import datetime
import hashlib
import json
import platform
import time
from pathlib import Path
from types import SimpleNamespace
import flint
from flint import arb, arb_mat, ctx
import validate_motor_cycle as base
import validate_motor_paired_cycle as paired
import certify_motor_saddle_spectrum as spectral
from validate_motor_cycle import (ROOT, ONE, ZERO, PolynomialODE, integrate,
    ball_record, matrix_record, midpoint_matrix, read_ball)
from validate_motor_paired_cycle import simplify, identity, extract, schur_certificate

K=5
DIM=14
MULTIPLICITIES=(2,2,1,2,1)
PERMUTATION=(3,0,4,1,2)
FREE=tuple(j for j in range(DIM) if j != 2*K)

def weights():
    return [arb(m)/8 for m in MULTIPLICITIES]

def positions(state):
    w=weights()
    return list(state[2*K:DIM])+[-sum((w[j]*state[2*K+j] for j in range(K-1)),ZERO)/w[-1]]

def permutation_matrix():
    p=arb_mat(DIM,DIM)
    w=weights()
    for j,k in enumerate(PERMUTATION):
        p[j,k]=1
        p[K+j,K+k]=1
        if j < K-1:
            if k < K-1:
                p[2*K+j,2*K+k]=1
            else:
                for l in range(K-1):
                    p[2*K+j,2*K+l]=-w[l]/w[-1]
    assert p**6 == identity(DIM)
    assert all(MULTIPLICITIES[j] == MULTIPLICITIES[k] for j,k in enumerate(PERMUTATION))
    return p

class FiveODE(PolynomialODE):
    def __init__(self,rho,r,with_core=True,with_transverse=False,order=48,generic_jets=False):
        self.rho,self.r,self.order,self.generic_jets=rho,r,order,generic_jets
        w=weights()
        velocities=[]
        for j in range(K):
            v=[(r,K+j)]+[(-r*w[k],K+k) for k in range(K)]
            v += [(-rho,2*K+j)] if j < K-1 else [(rho*w[k]/w[-1],2*K+k) for k in range(K-1)]
            velocities.append(v)
        terms=[]
        for j in range(K):
            terms.append(simplify([(-ONE,j,None),(ONE,None,None)]+[(-c,k,K+j) for c,k in velocities[j]]))
        for j in range(K):
            terms.append(simplify([(-ONE,K+j,None),(ONE,None,None)]+[(c,k,j) for c,k in velocities[j]]))
        for j in range(K-1):
            terms.append(simplify([(c,k,None) for c,k in velocities[j]]))
        self.field_terms=terms
        self.terms=list(terms)
        self.offsets={}
        self.matrix_jac={}
        if with_core:
            jac=[[[] for _ in range(DIM)] for _ in range(DIM)]
            for row,ts in enumerate(terms):
                for c,i,j in ts:
                    if i is not None: jac[row][i].append((c,j))
                    if j is not None: jac[row][j].append((c,i))
            self.add_matrix("core",jac,DIM)
        if with_transverse:
            for j,m in enumerate(MULTIPLICITIES):
                if m == 1: continue
                v=velocities[j]
                jac=[[[(-ONE,None)],[(-c,k) for c,k in v]+[(-r,K+j)],[(rho,K+j)]],
                     [v,[(-ONE,None),(r,j)],[(-rho,j)]],
                     [[],[(r,None)],[(-rho,None)]]]
                self.add_matrix(str(j),jac,3)
        self.row_disk_upper=ZERO

    add_matrix=paired.PairedODE.add_matrix
    initial=paired.PairedODE.initial

    def lipschitz(self,x,duration,unused_weights):
        matrix_weight=max(ONE,max((v.abs_upper() for v in x[DIM:]),default=ONE))
        scales=[ONE]*DIM+[matrix_weight]*(len(x)-DIM)
        rows=[]
        for row,ts in enumerate(self.terms):
            jac={}
            for c,i,j in ts:
                if i is not None: jac[i]=jac.get(i,ZERO)+c*(ONE if j is None else x[j])
                if j is not None: jac[j]=jac.get(j,ZERO)+c*(ONE if i is None else x[i])
            rows.append((sum((v.abs_upper()*scales[i] for i,v in jac.items()),ZERO)/scales[row]).upper())
        return duration.abs_upper()*max(rows)

    def jets(self,x,duration,order):
        if order == self.order+1:
            self.row_disk_upper=max(self.row_disk_upper,max((x[j]**2+x[K+j]**2).upper() for j in range(K)))
        if self.generic_jets:
            return super().jets(x,duration,order)
        core=PolynomialODE.jets(SimpleNamespace(terms=self.field_terms),x[:DIM],duration,order)
        coefficients=list(core)
        for name,(off,size) in self.offsets.items():
            jac=self.matrix_jac[name]
            jac_jets=[arb_mat([[sum((c*(ONE if k==0 else ZERO) if q is None else c*core[q][k]
                for c,q in jac[i][j]),ZERO) for j in range(size)] for i in range(size)]) for k in range(order)]
            matrices=[arb_mat([[x[off+size*i+j] for j in range(size)] for i in range(size)])]
            for k in range(order):
                product=arb_mat(size,size)
                for l in range(k+1): product+=jac_jets[l]*matrices[k-l]
                matrices.append((duration/(k+1))*product)
            coefficients += [[m[i,j] for m in matrices] for i in range(size) for j in range(size)]
        return coefficients

def flow(ode,state,duration,steps,order,audit):
    x=integrate(ode,ode.initial(state),duration,steps,order,audit)
    if audit is not None:
        audit[-1].pop("core_disk_upper",None)
        audit[-1]["five_cluster_row_disk_tube_upper"]=ball_record(ode.row_disk_upper)
        audit[-1]["jet_implementation"]="generic sparse" if ode.generic_jets else "Arb matrix convolution"
    return x

def initial_state(q):
    state=[ZERO]*DIM
    for i,j in enumerate(FREE): state[j]=q[i]
    return state

def segment_states(q,segments):
    states=[initial_state(q)]
    for s in range(1,segments):
        start=DIM-1+(s-1)*DIM
        states.append(list(q[start:start+DIM]))
    return states

def shooting(q,rho,r,steps,order,audit=None,generic=False,segments=1):
    states=segment_states(q,segments)
    size=DIM*segments
    if len(q)!=size: raise ValueError("Wrong multiple-shooting dimension")
    p=permutation_matrix()
    residual,jac=arb_mat(size,1),arb_mat(size,size)
    for s,state in enumerate(states):
        ode=FiveODE(rho,r,order=order,generic_jets=generic)
        x=flow(ode,state,q[-1],steps,order,audit)
        m=extract(x,ode.offsets["core"])
        f=FiveODE(rho,r,with_core=False).value(x[:DIM],ONE)
        last=s==segments-1
        target=p*arb_mat([[v] for v in states[0]]) if last else arb_mat([[v] for v in states[s+1]])
        for i in range(DIM):
            row=DIM*s+i
            residual[row,0]=x[i]-target[i,0]
            jac[row,size-1]=f[i]
            for j in (FREE if s==0 else range(DIM)):
                col=FREE.index(j) if s==0 else DIM-1+(s-1)*DIM+j
                jac[row,col]+=m[i,j]
            if last:
                for col,j in enumerate(FREE): jac[row,col]-=p[i,j]
            else:
                jac[row,DIM-1+s*DIM+i]-=1
    return residual,jac


# Sparse field and multiple-shooting layout adapted from the preserved saddle source
# SHA-256 01128b550bbae6b22f4c8aaf7c828556ed7a7e06ebb37d899e627f2ccd4c6364. No original module is modified.
# The relative flow has order-six row closure. No published module is mutated.

def spectrum(box,rho,r,args,audit):
    relative_flow=identity(DIM)
    pair_flow={j:identity(3) for j in (0,1,3)}
    for state in segment_states(box,args.segments):
        ode=FiveODE(rho,r,with_transverse=True,order=args.order)
        x=flow(ode,state,box[-1],args.steps,args.order,audit)
        relative_flow=extract(x,ode.offsets['core'])*relative_flow
        for j in pair_flow: pair_flow[j]=extract(x,ode.offsets[str(j)])*pair_flow[j]
    relative=permutation_matrix()**5*relative_flow
    core=spectral.similarity_schur(relative,args.seed,arb('0.9'))
    sharp_core=spectral.similarity_schur(relative,args.seed,arb('0.7'))
    contrasts={}
    for row in (0,1,3):
        product=identity(3);cursor=row;path=[]
        for _ in range(3):
            path.append(cursor);product=pair_flow[cursor]*product;cursor=PERMUTATION[cursor]
        assert cursor==row
        full=product**2;poly=full.charpoly()
        ordinary=schur_certificate([poly[k] for k in range(4)],arb('0.5'))
        sharp=schur_certificate([poly[k] for k in range(4)],arb('0.13'))
        contrasts[str(row)]=dict(three_relative_intervals=path,half_return=matrix_record(product),
            full_return=matrix_record(full),stability=ordinary,sharper_radius_test=sharp)
    passed=bool(core['passed'] and all(v['stability']['passed'] for v in contrasts.values()))
    return dict(passed=passed,relative_flow=matrix_record(relative_flow),relative_return=matrix_record(relative),
        quotient_relative_stability=core,quotient_sharper_relative_radius=sharp_core,transverse=contrasts,
        quotient_full_return_identity='M=(P^-1 Dphi_L)^6; exact flow eigenvalue 1 is deflated before Schur testing.',
        transverse_identity='For pair0: full M=(B1 B3 B0)^2. The other two full matrices are cyclic conjugates.',
        floquet_count=dict(unstable=0,neutral=1,stable=22) if passed else None,
        all_nonflow_moduli_below_point13=bool(sharp_core['passed'] and all(v['sharper_radius_test']['passed'] for v in contrasts.values())))


def prove(args,record):
    ctx.prec=args.bits
    rho=arb(129)/55
    r=1+arb('19/20')/arb('32001/500000').sqrt()
    seed=json.loads(args.seed.read_text())
    q=[arb(str(seed['states'][0][j])).mid() for j in FREE]
    q += [arb(str(value)).mid() for state in seed['states'][1:] for value in state]
    q += [arb(str(seed['segment_duration'])).mid()]
    size=DIM*args.segments
    if len(q)!=size:raise ValueError('Seed does not match shooting segments')
    if args.verify_input:
        old=json.loads(args.verify_input.read_text())
        q=[read_ball(v).mid() for v in old['shooting_center']]
        record['center_input_sha256']=hashlib.sha256(args.verify_input.read_bytes()).hexdigest()
    audit=[]
    record.update(exact_rho='129/55',parameter_rho=ball_record(rho),parameter_r=ball_record(r),
        seed_path=str(args.seed),seed_sha256=hashlib.sha256(args.seed.read_bytes()).hexdigest(),
        multiplicities=list(MULTIPLICITIES),relative_permutation=list(PERMUTATION),relative_permutation_order=6,
        shooting_dimension=size,shooting_segments=args.segments,integrations=audit,newton_history=[])
    for iteration in range(args.newton):
        f,j=shooting(q,rho,r,args.steps,args.order,audit,segments=args.segments)
        correction=midpoint_matrix(j).inv()*midpoint_matrix(f)
        record['newton_history'].append(dict(iteration=iteration,residual=matrix_record(f)))
        q=[(v-correction[k,0]).mid() for k,v in enumerate(q)]
        record['shooting_center']=[ball_record(v) for v in q]
        print('Newton',iteration,'residual',max(f[k,0].abs_upper() for k in range(size)),flush=True)
    center_f,center_j=shooting(q,rho,r,args.steps,args.order,audit,segments=args.segments)
    radius=arb(args.radius).upper();box=[arb(v,radius) for v in q]
    _,box_j=shooting(box,rho,r,args.steps,args.order,audit,segments=args.segments)
    c=midpoint_matrix(midpoint_matrix(center_j).inv());defect=identity(size)-c*box_j
    image=-c*center_f+defect*arb_mat([[arb(0,radius)] for _ in range(size)])
    ratios=[image[k,0].abs_upper()/radius for k in range(size)]
    contraction=max(sum((defect[i,j].abs_upper() for j in range(size)),ZERO).upper() for i in range(size))
    state=initial_state(box);xs=positions(state)
    velocity=FiveODE(rho,r,with_core=False).value(state,ONE)[2*K]
    distinct=[]
    for a in range(K):
        for b in range(a+1,K):
            differences=[state[a]-state[b],state[K+a]-state[K+b],xs[a]-xs[b]]
            distinct.append(dict(rows=[a,b],differences=[ball_record(v) for v in differences],
                                 passed=bool(any(not v.contains(0) for v in differences))))
    passed=bool(all(v<1 for v in ratios) and contraction<1 and not c.det().contains(0)
                and box[-1]>0 and velocity>0 and all(v['passed'] for v in distinct))
    record.update(shooting_center=[ball_record(v) for v in q],shooting_box=[ball_record(v) for v in box],
        segment_initial_boxes=[[ball_record(v) for v in state] for state in segment_states(box,args.segments)],
        center_residual=matrix_record(center_f),interval_jacobian=matrix_record(box_j),preconditioner=matrix_record(c),
        krawczyk_ratios=[ball_record(v) for v in ratios],krawczyk_contraction=ball_record(contraction),
        initial_state_box=[ball_record(v) for v in state],cluster_distinctions=distinct,
        section_velocity=ball_record(velocity),positive_period=ball_record(6*args.segments*box[-1]),
        krawczyk_passed=passed,status='periodic_five_cluster_cycle_certified' if passed else 'not_certified')
    if not passed:return
    print('Multiple-shooting inclusion certified',flush=True)
    record['spectrum']=spectrum(box,rho,r,args,audit)
    if record['spectrum']['passed']:record['status']='attracting_five_cluster_cycle_certified'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=Path,default=ROOT/'results/motor-five-cycle-proof/seed/seed.json')
    parser.add_argument('--bits',type=int,default=448)
    parser.add_argument('--steps',type=int,default=144)
    parser.add_argument('--order',type=int,default=48)
    parser.add_argument('--segments',type=int,choices=(3,),default=3)
    parser.add_argument('--radius',default='1e-38')
    parser.add_argument('--newton',type=int,default=3)
    parser.add_argument('--verify-input',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose an unused output')
    if args.bits<64 or args.steps<1 or args.order<2 or args.newton<0 or not arb(args.radius)>0:
        parser.error('Positive settings required')
    paths=[Path(__file__),Path(base.__file__),Path(paired.__file__),Path(spectral.__file__)]
    sources={p.name:p.read_bytes() for p in paths}
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        arguments={k:str(v) for k,v in vars(args).items()},python=platform.python_version(),python_flint=flint.__version__,
        source_sha256={name:hashlib.sha256(data).hexdigest() for name,data in sources.items()},
        limitations=['No least-period theorem','No validated parameter continuation or bifurcation theorem',
                     'No basin-size or historical-priority claim'])
    start=time.monotonic()
    try:prove(args,record)
    except Exception as error:
        record.update(status='error',error=repr(error));raise
    finally:
        record['wall_seconds']=time.monotonic()-start
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        for name,data in sources.items():args.output.with_name(args.output.stem+'.'+name).write_bytes(data)
        print(record.get('status'),args.output,flush=True)
    if record.get('status')!='attracting_five_cluster_cycle_certified':raise SystemExit(2)


if __name__=='__main__':main()
