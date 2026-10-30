"""Exact polynomial controls for the global-attractor boundedness argument.

Inequalities, compactness and Sobolev-space claims are proved in the note;
this script checks only the algebra used in those arguments.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import flint
from flint import fmpq_mpoly_ctx


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve prior evidence: choose a new output')
    ctx = fmpq_mpoly_ctx.get(('C','S','V','A','B','F','Z0','E','rho','g','U0','H'))
    C,S,V,A,B,F,Z0,E,rho,g,U0,H = ctx.gens()
    checks = {}

    def zero(name, expression):
        checks[name] = expression == 0
        if not checks[name]:
            raise AssertionError(f'Nonzero identity {name}: {expression}')

    Cdot, Sdot = -C-V*S+A, -S+V*C+B
    zero('transport_cancels_from_amplitude_balance',
         C*Cdot+S*Sdot+C*C+S*S-A*C-B*S)
    zero('forcing_cauchy_schwarz_sum_of_squares',
         (A*A+B*B)*(C*C+S*S)-(A*C+B*S)**2-(A*S-B*C)**2)
    R = F+1
    zero('harmonic_disk_strict_inward_bound', -R*R+F*R+R)
    radial_comparison = F+(Z0-F)*E
    radial_derivative = -(Z0-F)*E  # E'= -E.
    zero('radial_scalar_comparison_equation', radial_derivative+radial_comparison-F)
    G = 2*g*R  # g=|gamma|, R=F+1 in the analytic inequalities.
    scaled_filter = G+(rho*U0-G)*H
    scaled_filter_derivative = -rho*(rho*U0-G)*H  # H'= -rho H.
    zero('mechanical_filter_comparison_equation_times_rho',
         scaled_filter_derivative+rho*scaled_filter-rho*G)
    scaled_boundary = G+rho  # rho*M, where M=G/rho+1.
    zero('positive_displacement_face_outward_bound_times_rho',
         -rho*scaled_boundary+rho*G+rho*rho)
    zero('negative_displacement_face_outward_bound_times_rho',
         -(rho*scaled_boundary-rho*G)+rho*rho)
    for n in (1,2,8):
        local = fmpq_mpoly_ctx.get(tuple(['rho','gamma']+
            [f'X{j}' for j in range(n)]+[f'S{j}' for j in range(n)]))
        values = local.gens()
        rho_n, gamma_n = values[:2]
        xs, ss = values[2:2+n], values[2+n:]
        mean = sum(ss)/n
        velocities = [-rho_n*x-gamma_n*(s-mean) for x,s in zip(xs,ss)]
        zero(f'displacement_sum_constraint_N{n}', sum(velocities)+rho_n*sum(xs))
    here = Path(__file__).resolve()
    root = here.parents[2]
    prerequisite = root/'docs/motor-pde-lifting.md'
    report = root/'docs/motor-global-attractor.md'
    result = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        status='passed', all_zero=all(checks.values()), identities=checks,
        method='Exact rational multivariate polynomial identities in Python-FLINT',
        python_flint=flint.__version__, script_sha256=hashlib.sha256(here.read_bytes()).hexdigest(),
        prerequisite_note_sha256=hashlib.sha256(prerequisite.read_bytes()).hexdigest(),
        report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
        assumptions=['Finite fixed N>=1','rho>0','F=sqrt(A^2+B^2)>=0','g=abs(gamma)>=0',
            'Constant total switching rate and first-harmonic force/forcing'],
        scope='Polynomial equalities only. The analytic note supplies inequalities, uniform absorption, compactness, strong continuity and the complete-trajectory argument.',
        not_computed=['No trajectories','No PDE truncation','No chaos or finite-attractor exponential-rate assertion'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],identities=len(checks),all_zero=result['all_zero'])))


if __name__ == '__main__':
    main()
