"""Validate the two singleton contrast blocks needed for N=8k replication."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time
from flint import arb,arb_mat,ctx
import validate_motor_five_cycle as five
from validate_motor_cycle import ZERO,ONE,ROOT,ball_record,matrix_record,read_ball


class SingletonODE(five.FiveODE):
    def __init__(self,rho,r,order):
        super().__init__(rho,r,with_core=False,with_transverse=False,order=order)
        w=five.weights()
        for j in (2,4):
            velocity=[(r,5+j)]+[(-r*w[k],5+k) for k in range(5)]
            velocity += [(-rho,10+j)] if j<4 else [(rho*w[k]/w[4],10+k) for k in range(4)]
            jac=[[[(-ONE,None)],[(-c,k) for c,k in velocity]+[(-r,5+j)],[(rho,5+j)]],
                 [velocity,[(-ONE,None),(r,j)],[(-rho,j)]],
                 [[],[(r,None)],[(-rho,None)]]]
            self.add_matrix(str(j),jac,3)


def full_jacobian(state,rho,r):
    n=16;c=state[:n];s=state[n:2*n];x=list(state[2*n:])+[-sum(state[2*n:],ZERO)]
    velocity=[-rho*x[i]+r*(s[i]-sum(s,ZERO)/n) for i in range(n)]
    jac=arb_mat(47,47)
    for i in range(n):
        jac[i,i]=-1
        for j in range(n):
            projection=arb(int(i==j))-arb(1)/n
            jac[i,n+j]=-int(i==j)*velocity[i]-r*s[i]*projection
            jac[n+i,j]=int(i==j)*velocity[i]
            jac[n+i,n+j]=-int(i==j)+r*c[i]*projection
            if i<15:jac[32+i,n+j]=r*projection
        for j in range(15):
            embedding=arb(int(i==j)) if i<15 else -ONE
            jac[i,32+j]=rho*s[i]*embedding
            jac[n+i,32+j]=-rho*c[i]*embedding
            if i<15:jac[32+i,32+j]=-rho*int(i==j)
    return jac


def controls():
    errors=[]
    group_indices=[0]*4+[1]*4+[2]*2+[3]*4+[4]*2
    for trial in range(1,7):
        state=[arb((i+2)*trial-13)/31 for i in range(14)]
        rho=arb(3)+arb(trial)/17;r=arb(4)+arb(trial)/19
        positions=five.positions(state)
        full=[state[j] for j in group_indices]+[state[5+j] for j in group_indices]+[positions[j] for j in group_indices[:-1]]
        jac=full_jacobian(full,rho,r);ode=SingletonODE(rho,r,8)
        values=ode.value(ode.initial(state),ONE)
        for name,(a,b) in {'2':(8,9),'4':(14,15)}.items():
            embedding=arb_mat(47,3)
            embedding[a,0]=1;embedding[b,0]=-1
            embedding[16+a,1]=1;embedding[16+b,1]=-1
            if a<15:embedding[32+a,2]=1
            if b<15:embedding[32+b,2]=-1
            block=five.extract(values,ode.offsets[name])
            error=jac*embedding-embedding*block
            errors += [error[i,j].abs_upper() for i in range(47) for j in range(3)]
    return dict(passed=bool(max(errors)<arb('1e-95')),components=len(errors),maximum_error=ball_record(max(errors)),
        scope='Independent full47D N16 Jacobian and both singleton internal contrast embeddings, including the eliminated last X coordinate.')


def outside_real_root(poly):
    attempts=[]
    for power in range(1,21):
        radius=2**power
        for left,right in ((-radius,-1),(1,radius)):
            a,b=poly(arb(left)),poly(arb(right))
            passed=bool((a<0 and b>0) or (a>0 and b<0))
            attempts.append(dict(interval=[left,right],left=ball_record(a),right=ball_record(b),passed=passed))
            if passed:return dict(certified=True,root_bracket=attempts[-1],tests=attempts)
    return dict(certified=False,tests=attempts)


def run(args,record):
    ctx.prec=args.bits
    source=json.loads(args.input.read_text())
    if source.get('status')!='attracting_five_cluster_cycle_certified' or not source['krawczyk_passed']:
        raise ValueError('Requires the complete certified five-cluster input')
    if source['exact_rho']!='129/55' or source['multiplicities']!=[2,2,1,2,1] or source['relative_permutation']!=[3,0,4,1,2]:
        raise ValueError('Unexpected model or symmetry')
    record['independent_N16_embedding_control']=controls()
    if not record['independent_N16_embedding_control']['passed']:raise ArithmeticError('N16 embedding control failed')
    box=[read_ball(v) for v in source['shooting_box']]
    rho=arb(129)/55;r=1+arb('19/20')/arb('32001/500000').sqrt()
    blocks={j:five.identity(3) for j in (2,4)};audit=[]
    for state in five.segment_states(box,3):
        ode=SingletonODE(rho,r,args.order)
        end=five.flow(ode,state,box[-1],args.steps,args.order,audit)
        for j in blocks:blocks[j]=five.extract(end,ode.offsets[str(j)])*blocks[j]
    result={}
    for j,k in ((2,4),(4,2)):
        two=blocks[k]*blocks[j];full=two**3;poly=full.charpoly()
        stable=five.schur_certificate([poly[n] for n in range(4)])
        root=dict(certified=False) if stable['passed'] else outside_real_root(poly)
        result[str(j)]=dict(two_interval_path=[j,k],two_interval_return=matrix_record(two),
            full_return=matrix_record(full),characteristic_coefficients=[ball_record(poly[n]) for n in range(4)],
            unit_disk_schur=stable,separate_real_root_outside_unit_disk=root)
    stable=all(v['unit_disk_schur']['passed'] for v in result.values())
    unstable=any(v['separate_real_root_outside_unit_disk']['certified'] for v in result.values())
    record.update(status='replication_attraction_certified' if stable else 'replication_transverse_instability_certified' if unstable else 'singleton_stability_inconclusive',
        singleton_blocks=result,integrations=audit,parameter_rho=ball_record(rho),parameter_r=ball_record(r),
        inherited_positive_period=source['positive_period'],
        composition='M2=(B4 B2)^3; M4=(B2 B4)^3 over the inherited positive period T=18 segment durations.',
        scope='Exact replication has cluster sizes(2k,2k,k,2k,k), N=8k. Original three pair blocks inherit stability from the input. New singleton contrasts exist only for k>=2.',
        limitations=['No repeated periodic shooting; this extension depends on the immutable certified input.',
                     'Failure of a Schur test alone is inconclusive. Instability is asserted only after a separate strict real-root bracket outside[-1,1].',
                     'The fixed N8 theorem is unchanged; no least-period or new global-attractor claim.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'results/motor-five-cycle-proof/recheck-01.json')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--bits',type=int,default=512)
    parser.add_argument('--steps',type=int,default=160)
    parser.add_argument('--order',type=int,default=50)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Choose an unused output')
    paths=[Path(__file__),Path(five.__file__),Path(five.base.__file__),Path(five.paired.__file__),Path(five.spectral.__file__)]
    sources={p.name:p.read_bytes() for p in paths}
    record=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        arguments={k:str(v) for k,v in vars(args).items()},input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        source_sha256={name:hashlib.sha256(data).hexdigest() for name,data in sources.items()})
    start=time.monotonic()
    try:run(args,record)
    except Exception as error:
        record.update(status='error',error=repr(error));raise
    finally:
        record['wall_seconds']=time.monotonic()-start
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        for name,data in sources.items():args.output.with_name(args.output.stem+'.'+name).write_bytes(data)
        print(record.get('status'),args.output,flush=True)


if __name__=='__main__':main()
