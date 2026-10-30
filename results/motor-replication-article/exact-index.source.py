"""Certify the complete Floquet index of the replicated five-cluster cycle.

Uses preserved interval matrices, without replacing or narrowing their radii.
The negative root is bracketed by strict polynomial signs. Synthetic division
then encloses the exact remaining quadratic; real quadratic Jury inequalities
certify its two roots inside the unit disk. No floating eigenvalues are used.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from flint import arb, arb_mat, ctx
from validate_motor_cycle import ROOT, ball_record, read_ball


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def quadratic_jury(linear, constant, radius):
    a, b = linear / radius, constant / radius**2
    margins = [1 - b.abs_upper(), 1 + a + b, 1 - a + b]
    return dict(radius=ball_record(radius),
                margins=[ball_record(v) for v in margins],
                passed=bool(all(v > 0 for v in margins)))


def certify(matrix):
    p = matrix.charpoly()
    assert p.degree() == 3 and p[3] == 1
    # Independent explicit 3x3 coefficient expressions check the input layout.
    trace = sum((matrix[j, j] for j in range(3)), arb(0))
    minors = sum((matrix[i, i]*matrix[j, j]-matrix[i, j]*matrix[j, i]
                  for i in range(3) for j in range(i+1, 3)), arb(0))
    determinant = (matrix[0, 0]*(matrix[1, 1]*matrix[2, 2]-matrix[1, 2]*matrix[2, 1])
                   - matrix[0, 1]*(matrix[1, 0]*matrix[2, 2]-matrix[1, 2]*matrix[2, 0])
                   + matrix[0, 2]*(matrix[1, 0]*matrix[2, 1]-matrix[1, 1]*matrix[2, 0]))
    coefficient_control = all(x.contains(0) for x in
                              [p[2]+trace, p[1]-minors, p[0]+determinant])
    left, right = arb(-32), arb(-1)
    assert p(left) < 0 and p(right) > 0
    steps = 0
    # Stop well before input uncertainty could dominate the sign margins.
    for _ in range(80):
        if right-left < arb('1e-8'):
            break
        middle = (left+right)/2
        value = p(middle)
        if value < 0:
            left = middle
        elif value > 0:
            right = middle
        else:
            break
        steps += 1
    root = arb((left+right)/2, (right-left)/2)
    linear = p[2]+root
    constant = p[1]+root*linear
    stable = quadratic_jury(linear, constant, arb(1))
    sharp = quadratic_jury(linear, constant, arb('0.13'))
    passed = bool(coefficient_control and p(left)<0 and p(right)>0
                  and right < -1 and stable['passed'])
    return dict(passed=passed, characteristic_coefficients=[ball_record(p[j]) for j in range(4)],
                independent_coefficient_control=coefficient_control,
                root_left=ball_record(left), root_right=ball_record(right),
                polynomial_left=ball_record(p(left)), polynomial_right=ball_record(p(right)),
                root_enclosure=ball_record(root), bisection_steps=steps,
                quotient_linear=ball_record(linear), quotient_constant=ball_record(constant),
                unit_disk_jury=stable, sharper_jury=sharp,
                count=dict(unstable=1, neutral=0, stable=2) if passed else None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT/'results/motor-five-replication-proof/singleton-blocks.json')
    parser.add_argument('--cycle', type=Path, default=ROOT/'results/motor-five-cycle-proof/recheck-01.json')
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new output to preserve evidence')
    ctx.prec=768
    data=json.loads(args.input.read_text())
    cycle=json.loads(args.cycle.read_text())
    assert data['status']=='replication_transverse_instability_certified'
    assert data['input_sha256']==sha(args.cycle)
    assert cycle['status']=='attracting_five_cluster_cycle_certified'
    assert cycle['spectrum']['all_nonflow_moduli_below_point13']
    assert set(data['singleton_blocks'])=={'2','4'}
    blocks={name:certify(arb_mat([[read_ball(v) for v in row] for row in value['full_return']]))
            for name,value in data['singleton_blocks'].items()}
    passed=all(value['passed'] for value in blocks.values())
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                status='exact_replication_index_certified' if passed else 'not_certified',
                bits=ctx.prec, input_sha256=sha(args.input), cycle_sha256=sha(args.cycle),
                source_sha256=sha(__file__), blocks=blocks,
                conclusion='For every integer k>=1 at fixed normalized parameters and weights, N=8k has 2k-2 unstable, 1 neutral and 22k stable Floquet directions.' if passed else None,
                decomposition='13 stable quotient directions; 3(2k-1) copies of stable 3D paired blocks; 2(k-1) copies of singleton blocks with index (1,0,2).',
                mathematical_basis='For a real root u in the certified bracket, p(z)=(z-u)[z^2+(p2+u)z+p1+u(p2+u)]. Interval Jury inequalities hold for the entire enclosing quadratic coefficient box.',
                limitations=['Depends on the preserved validated periodic orbit and variational matrix enclosures.',
                             'Spectral replay, not a fresh integration.',
                             'No least period, nonlinear destination, physical N-variation or historical-priority assertion.'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    args.output.with_suffix('.source.py').write_bytes(Path(__file__).read_bytes())
    print(record['status'])
    for name,value in blocks.items():
        print(name,value['root_enclosure']['display'],value['unit_disk_jury']['passed'],value['sharper_jury']['passed'])
    if not passed:
        raise SystemExit(2)


if __name__=='__main__':
    main()
