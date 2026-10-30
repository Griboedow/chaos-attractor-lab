"""Exact polynomial checks for the weighted motor-volume identity."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import flint
from flint import arb, ctx, fmpq_mpoly_ctx
from validate_motor_cycle import ball_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output; preserve prior evidence')
    checks = {}

    def zero(name, expression):
        checks[name] = bool(expression == 0)
        if not checks[name]:
            raise ArithmeticError(name + ': ' + str(expression))

    for n in (1, 2, 3, 4, 8):
        names = ['rho','r','b'] + [f'c{j}' for j in range(n)] + [f's{j}' for j in range(n)] + [f'X{j}' for j in range(n-1)]
        pc = fmpq_mpoly_ctx.get(tuple(names))
        variables = pc.gens()
        rho, r, b = variables[:3]
        cs, ss = variables[3:3+n], variables[3+n:3+2*n]
        independent_x = variables[3+2*n:]
        xs = list(independent_x) + [-sum(independent_x, pc.constant(0))]
        theta = pc.constant(n-1)/n
        mean = sum(ss)/n
        vs = [-rho*x+r*(s-mean) for x,s in zip(xs,ss)]
        cdot = [1-c-v*s for c,s,v in zip(cs,ss,vs)]
        sdot = [b-s+v*c for c,s,v in zip(cs,ss,vs)]
        field = cdot+sdot+vs[:-1]
        state_names = names[3:]
        trace = sum((f.derivative(name) for f,name in zip(field,state_names)), pc.constant(0))
        expected_trace = -2*n-rho*(n-1)+r*theta*sum(cs)
        zero(f'N{n}:zero_sum_velocity', sum(vs))
        zero(f'N{n}:full_reduced_trace', trace-expected_trace)
        speed2 = sum(v*v for v in vs)
        zero(f'N{n}:mechanical_cancellation', r*sum(v*s for v,s in zip(vs,ss))-speed2-rho*sum(x*v for x,v in zip(xs,vs)))
        phi = theta*(r*sum(cs)+rho*sum(x*x for x in xs)/2)
        phidot = sum((phi.derivative(name)*f for f,name in zip(field,state_names)), pc.constant(0))
        zero(f'N{n}:scalar_derivative', phidot-r*theta*(n-sum(cs))+theta*speed2)
        constant = -2*n+(n-1)*(r-rho)
        zero(f'N{n}:weighted_divergence', trace+phidot-constant+theta*speed2)

    ctx.prec = 256
    r = 1+(arb(19)/20)/(arb(32001)/500000).sqrt()
    parameters = {}
    for label, rho in [('stable_four_pairs', arb('3.3652183995790423')), ('six_cluster_saddle',arb(817)/250), ('activity_path_candidate',arb(1291)/440)]:
        constant = -16+7*(r-rho)
        parameters[label] = dict(rho=ball_record(rho), r=ball_record(r), exponent_sum_upper_bound=ball_record(constant), strictly_negative=bool(constant<0))
    here = Path(__file__).resolve()
    report = here.parents[1]/'docs/motor-volume-contraction.md'
    record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed',
        identities=checks, parameters=parameters, python_flint=flint.__version__, precision_bits=ctx.prec,
        source_sha256=hashlib.sha256(here.read_bytes()).hexdigest(), report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
        scope='Exact polynomial and interval-constant controls. The report supplies the analytic weighted-measure and invariant-set argument; no trajectories or chaos certification.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    snapshot = args.output.with_suffix('.source.py')
    if snapshot.exists():
        raise FileExistsError('Source snapshot already exists')
    snapshot.write_bytes(here.read_bytes())
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('passed',len(checks),'identities',{k:v['exponent_sum_upper_bound']['display'] for k,v in parameters.items()})


if __name__ == '__main__':
    main()
