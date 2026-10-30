"""Exact rest-state linearization and interval source-parameter controls.

No trajectories are computed. The unstable-manifold and connectedness
arguments are analytic statements in docs/motor-global-geometry.md.
"""
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import permutations
import json
from pathlib import Path

import flint
from flint import arb, ctx, fmpq_mat, fmpq_mpoly_ctx


def ball_record(x):
    m,e = x.mid().man_exp()
    r,s = x.rad().man_exp()
    return dict(display=x.str(25), midpoint=[str(m),int(e)], radius=[str(r),int(s)])


def determinant3(matrix):
    total = 0
    for p in permutations(range(3)):
        parity = sum(p[i] > p[j] for i in range(3) for j in range(i+1,3))
        total += (-1)**parity * matrix[0][p[0]] * matrix[1][p[1]] * matrix[2][p[2]]
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve prior evidence: choose a new output')
    checks = {}

    def zero(name, expression):
        checks[name] = bool(expression == 0)
        if not checks[name]:
            raise ArithmeticError(name+': '+str(expression))

    pc = fmpq_mpoly_ctx.get(('lam','rho','r','b'))
    lam,rho,r,b = pc.gens()
    block = [[-1,-b*r,b*rho],[0,r-1,-rho],[0,r,-rho]]
    charpoly = determinant3([[lam*int(i==j)-block[i][j] for j in range(3)] for i in range(3)])
    quadratic = lam*lam+(rho+1-r)*lam+rho
    zero('difference_block_characteristic_polynomial',charpoly-(lam+1)*quadratic)
    zero('characteristic_polynomial_independent_of_b',charpoly.derivative('b'))
    zero('minus_one_not_a_quadratic_root_when_r_positive',quadratic.subs({'lam':-1})-r)

    decompositions = []
    for n in (1,2,3,8):
        names = ['rho','r','b']+[f'u{j}' for j in range(n)]+[f'v{j}' for j in range(n)]+[f'X{j}' for j in range(n-1)]
        local = fmpq_mpoly_ctx.get(tuple(names))
        variables = local.gens()
        rho_n,r_n,b_n = variables[:3]
        us,vs = variables[3:3+n],variables[3+n:3+2*n]
        xs = list(variables[3+2*n:])
        xs += [-sum(xs,local.constant(0))]
        speeds = [-rho_n*x+r_n*(v-sum(vs)/n) for x,v in zip(xs,vs)]
        # c=1+u and s=b+v: exact shifted nonlinear equations.
        field = [-u-speed*(b_n+v) for u,v,speed in zip(us,vs,speeds)]
        field += [-v+speed*(1+u) for u,v,speed in zip(us,vs,speeds)]
        field += speeds[:-1]
        origin = {name:0 for name in names[3:]}
        d = 3*n-1
        for i,f in enumerate(field):
            zero(f'N{n}:rest_field_{i}',f.subs(origin))
        jacobian = [[f.derivative(name).subs(origin) for name in names[3:]] for f in field]
        basis = [[0]*d for _ in range(d)]
        for row in range(n):
            basis[row][0] = 1
            basis[n+row][1] = 1
        expected = [[local.constant(0) for _ in range(d)] for _ in range(d)]
        expected[0][0] = expected[1][1] = local.constant(-1)
        difference = [[-1,-b_n*r_n,b_n*rho_n],[0,r_n-1,-rho_n],[0,r_n,-rho_n]]
        for j in range(n-1):
            col = 2+3*j
            basis[j][col],basis[n-1][col] = 1,-1
            basis[n+j][col+1],basis[2*n-1][col+1] = 1,-1
            basis[2*n+j][col+2] = 1
            for a in range(3):
                for b_index in range(3):
                    expected[col+a][col+b_index] = difference[a][b_index]
        det = fmpq_mat(basis).det()
        assert abs(det) == n*n
        for i in range(d):
            for j in range(d):
                left = sum((jacobian[i][k]*basis[k][j] for k in range(d)),local.constant(0))
                right = sum((basis[i][k]*expected[k][j] for k in range(d)),local.constant(0))
                zero(f'N{n}:full_J_basis_entry_{i}_{j}',left-right)
        decompositions.append(dict(N=n,dimension=d,basis_determinant=str(det),
            common_modes=2,difference_blocks=n-1,coefficient_identities=d*d,passed=True))

    ctx.prec = 256
    r_source = 1+arb('19/20')/arb('32001/500000').sqrt()
    parameters = {}
    for label,rho_source in (('six_cluster_saddle',arb('817/250')),
                              ('activity_path_candidate',arb('1291/440'))):
        gap = r_source-rho_source-1
        discriminant = gap*gap-4*rho_source
        volume_constant = -16+7*(r_source-rho_source)
        passed = bool(rho_source>0 and gap>0 and discriminant<0 and volume_constant<0)
        assert passed
        parameters[label] = dict(N=8,b='1',r=ball_record(r_source),rho=ball_record(rho_source),
            r_minus_rho_minus_one=ball_record(gap),quadratic_discriminant=ball_record(discriminant),
            unstable_real_part=ball_record(gap/2),
            unstable_imaginary_magnitude=ball_record((-discriminant).sqrt()/2),
            weighted_volume_constant=ball_record(volume_constant),
            unstable_real_dimension=14,stable_real_dimension=9,
            global_attractor_dimension_lower_bound=14,phase_space_dimension=23,
            positive_gap_and_zero_volume_condition_certified=passed)
    here = Path(__file__).resolve()
    root = here.parents[2]
    notes = [root/'docs/motor-global-geometry.md',root/'docs/motor-global-attractor.md',
             root/'docs/motor-volume-contraction.md']
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),status='passed',
        python_flint=flint.__version__,precision_bits=ctx.prec,
        identities=checks,identities_count=len(checks),decompositions=decompositions,
        parameters=parameters,source_sha256=hashlib.sha256(here.read_bytes()).hexdigest(),
        note_sha256={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in notes},
        scope='Exact polynomial rest linearization and block decomposition controls, plus Arb inequalities at exact source constants. General-N decomposition, unstable-manifold inclusion and connectedness are proved analytically in the note.',
        limitations=['No new trajectories','No construction or interval parameterization of an unstable manifold',
            'No exact attractor dimension or noninteger dimension claim','No chaos or scientific-priority assertion'])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=result['status'],identities=len(checks),
        source_gaps={name:record['r_minus_rho_minus_one']['display'] for name,record in parameters.items()})))


if __name__ == '__main__':
    main()
