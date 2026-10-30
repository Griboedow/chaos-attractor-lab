"""Enclose leading saddle rates from an existing validated monodromy."""
import argparse
import datetime
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from flint import arb, ctx
import validate_motor_cycle as base
import validate_motor_paired_cycle as paired


def evaluate(coefficients, x):
    value = arb(0)
    for a in reversed(coefficients):
        value = value*x+a
    return value


def rational(x):
    return arb(x.numerator)/arb(x.denominator)


def isolate(coefficients, left, right):
    left, right = Fraction(left), Fraction(right)
    if not (evaluate(coefficients, rational(left)) < 0 and evaluate(coefficients, rational(right)) > 0):
        raise ArithmeticError('Strict opposite endpoint signs required')
    history = []
    for _ in range(80):
        middle = (left+right)/2
        value = evaluate(coefficients, rational(middle))
        sign = -1 if value < 0 else 1 if value > 0 else 0
        history.append(dict(midpoint=str(middle), sign=sign))
        if sign == 0:
            break
        if sign < 0:
            left = middle
        else:
            right = middle
    root = rational(left).union(rational(right))
    # The intermediate value theorem gives an exact real root in this bracket.
    # Division by an interval containing that root encloses the exact quotient.
    quotient = [arb(0)]*(len(coefficients)-1)
    quotient[-1] = coefficients[-1]
    for k in range(len(quotient)-2, -1, -1):
        quotient[k] = coefficients[k+1]+root*quotient[k+1]
    remainder = coefficients[0]+root*quotient[0]
    schur = paired.schur_certificate(quotient, arb('0.11'))
    passed = bool(remainder.contains(0) and schur['passed'])
    return root, dict(passed=passed, exact_left=str(left), exact_right=str(right), root=base.ball_record(root),
        left_value=base.ball_record(evaluate(coefficients, rational(left))),
        right_value=base.ball_record(evaluate(coefficients, rational(right))),
        bisections=history, division_remainder=base.ball_record(remainder), other_roots_schur=schur)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output; preserve evidence')
    data = json.loads(args.input.read_text(encoding='utf-8'))
    if data.get('status') != 'six_cluster_saddle_certified':
        raise ValueError('Requires an existing complete saddle certificate')
    ctx.prec = int(data['arguments']['bits'])
    if not data['spectrum']['similarity_certificate']['passed']:
        raise ValueError('Quotient stability is required')
    core_radius = base.read_ball(data['spectrum']['similarity_certificate']['nonflow_schur']['radius'])
    records, roots = {}, {}
    for cluster, lo, hi in [('2', '0.8005', '0.8007'), ('5', '1.12', '1.13')]:
        coefficients = [base.read_ball(x) for x in data['spectrum']['transverse'][cluster]['characteristic_coefficients']]
        roots[cluster], records[cluster] = isolate(coefficients, lo, hi)
    period = base.read_ball(data['positive_period'])
    expanding = roots['5'].log()/period
    contracting = -roots['2'].log()/period
    ratio = -roots['2'].log()/roots['5'].log()
    gates = dict(two_simple_isolated_roots=all(x['passed'] for x in records.values()),
        expanding_root_greater_than_one=roots['5']>1,
        contracting_root_between_zero_and_one=roots['2']>0 and roots['2']<1,
        contracting_root_dominates_quotient=roots['2']>core_radius,
        contracting_root_dominates_other_contrasts=roots['2']>arb('0.11'),
        positive_period=period>0, ratio_greater_than_one=ratio>1)
    sources = [Path(__file__), Path(base.__file__), Path(paired.__file__)]
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    output = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status='leading_saddle_rates_certified' if all(gates.values()) else 'rate_certificate_failed',
        input_path=str(args.input), input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        precision_bits=ctx.prec, source_sha256=hashes, roots=records,
        expanding_rate=base.ball_record(expanding), leading_contraction_magnitude=base.ball_record(contracting),
        contraction_expansion_ratio=base.ball_record(ratio), gates=gates,
        interpretation='Positive real leading stable and unstable multipliers, with ratio of contraction to expansion strictly above one.',
        limitations=['No repeated integration.', 'Local rates do not prove heteroclinic connection, global return map, network attraction, or an exact switching-time recurrence.'])
    planned = [args.output.with_name(args.output.stem+'.'+p.name) for p in sources]
    if any(p.exists() for p in planned):
        raise FileExistsError('A planned source snapshot already exists')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for source, snapshot in zip(sources, planned):
        snapshot.write_bytes(source.read_bytes())
    args.output.write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
    print(output['status'], {k:output[k]['display'] for k in ('expanding_rate','leading_contraction_magnitude','contraction_expansion_ratio')})
    if not all(gates.values()):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
