"""Certify an extra FastICA preimage in T and exclude it from F(T).

All deciding arithmetic is exact rational arithmetic. The proposal is a
rationalized finite trajectory sample, not a point certified to belong to K.
The written argument is in docs/fastica-preimage-rejection.md.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
E = Q(1, 25)
K2 = Q(4, 25)
PREVIOUS = list(map(Q, ['-0.38715568266959405', '-0.9020845179494923',
                       '-0.19040449850479912', '0.01000633807187778']))
CENTER = list(map(Q, ['-0.43931892', '-0.29018779', '-0.03294978']))
RADIUS = Q(1, 10**7)


class Interval:
    """Closed rational interval; no binary floating-point decisions."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = Q(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def lift(value):
        return value if isinstance(value, Interval) else Interval(value)

    def __add__(self, other):
        other = self.lift(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.lift(other)

    def __rsub__(self, other):
        return self.lift(other) + -self

    def __mul__(self, other):
        other = self.lift(other)
        values = (self.lo * other.lo, self.lo * other.hi,
                  self.hi * other.lo, self.hi * other.hi)
        return Interval(min(values), max(values))

    __rmul__ = __mul__

    def __pow__(self, degree):
        assert isinstance(degree, int) and degree >= 0
        if degree == 0:
            return Interval(1)
        if degree % 2:
            return Interval(self.lo**degree, self.hi**degree)
        lower = 0 if self.contains_zero() else min(self.lo**degree, self.hi**degree)
        return Interval(lower, max(self.lo**degree, self.hi**degree))

    def contains_zero(self):
        return self.lo <= 0 <= self.hi

    def absolute_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def pair(self):
        return [str(self.lo), str(self.hi)]


def raw(w):
    """Homogeneous numerator, independently written from the published P."""
    A, B, C, D = w
    e = E
    return [
        A**3 - 3*A*B**2 + 18*e*A**2*C + e**2*(36*A*C**2 - 18*A*D**2) + 18*e**3*C**3,
        B**3 - 3*B*A**2 + 18*e*B**2*D + e**2*(36*B*D**2 - 18*B*C**2) + 18*e**3*D**3,
        6*e*A**3 + e**2*(36*A**2*C - 18*B**2*C) + 54*e**3*A*C**2 - 18*e**4*C**3 - 108*e**4*C*D**2,
        6*e*B**3 + e**2*(36*B**2*D - 18*A**2*D) + 54*e**3*B*D**2 - 18*e**4*D**3 - 108*e**4*D*C**2,
    ]


def evaluate(poly, values):
    result = Interval(0)
    for powers, coefficient in poly.terms():
        term = Interval(Q(coefficient))
        for value, power in zip(values, powers):
            term *= value**power
        result += term
    return result


def certify_root(target):
    variables = sp.symbols('A B C D')
    parameter = sp.Symbol('e')
    prior_path = ROOT / 'results/fastica-tripling/exact-check.json'
    prior = json.loads(prior_path.read_text(encoding='utf-8'))
    polynomial = sp.sympify(prior['cumulant_polynomial'],
                           locals=dict(zip(['A', 'B', 'C', 'D', 'e'], [*variables, parameter])))
    for i, value in enumerate(raw(variables)):
        assert sp.expand(value - sp.diff(polynomial, variables[i]).subs(parameter, sp.Rational(E))/4) == 0
    x = sp.symbols('a c d')
    h = raw([x[0], sp.Integer(-1), x[1], x[2]])
    h0 = raw(PREVIOUS)
    f = sp.Matrix([h[i]*h0[1] - h[1]*h0[i] for i in [0, 2, 3]])
    jacobian = f.jacobian(x)
    at_center = dict(zip(x, CENTER))
    inverse = jacobian.subs(at_center).inv()
    assert inverse * jacobian.subs(at_center) == sp.eye(3)
    base = sp.Matrix(CENTER) - inverse*f.subs(at_center)
    derivative = [[evaluate(sp.Poly(jacobian[i, j], *x), target) for j in range(3)] for i in range(3)]
    defect = [[Interval(int(i == j)) - sum(Q(inverse[i, k])*derivative[k][j] for k in range(3))
               for j in range(3)] for i in range(3)]
    contraction = max(sum(v.absolute_upper() for v in row) for row in defect)
    assert contraction < 1
    image = [Interval(Q(base[i])) + sum(defect[i][j]*Interval(-RADIUS, RADIUS) for j in range(3))
             for i in range(3)]
    margins = [min(image[i].lo - target[i].lo, target[i].hi - image[i].hi) for i in range(3)]
    assert min(margins) > 0
    # Same direction, not merely equality of projective classes.
    boxed_h = raw([target[0], Interval(-1), target[1], target[2]])
    assert h0[1] < 0 and boxed_h[1].hi < 0
    previous_margin = K2*(PREVIOUS[0]**2 + PREVIOUS[1]**2) - PREVIOUS[2]**2 - PREVIOUS[3]**2
    root_margin = K2*(1 + (target[0]**2).lo) - (target[1]**2 + target[2]**2).hi
    assert min(previous_margin, root_margin) > 0
    previous_chart_c = PREVIOUS[2]/(-PREVIOUS[1])
    separation = previous_chart_c - target[1].hi
    assert separation > 0
    return {
        'method': 'Exact inverse at the center; contraction and strict self-map of the closed rational box.',
        'polynomial_equations': [str(sp.expand(v)) for v in f],
        'inverse_center_jacobian': [[str(inverse[i, j]) for j in range(3)] for i in range(3)],
        'contraction_infinity_bound': str(contraction),
        'contraction_decimal_for_display': float(contraction),
        'fixed_point_map_image': [v.pair() for v in image],
        'strict_self_map_margins': list(map(str, margins)),
        'previous_h_B': str(h0[1]), 'root_h_B_interval': boxed_h[1].pair(),
        'previous_tube_margin': str(previous_margin), 'root_box_tube_margin': str(root_margin),
        'separation_in_C_over_minus_B': str(separation),
        'previous_algebra_sha256': hashlib.sha256(prior_path.read_bytes()).hexdigest(),
    }


def exclude_predecessor(target, maximum_depth=32):
    # After positive scaling and, if necessary, the antipodal change, every
    # input in T has A=1 or B=1, |t|<=1 and |u|,|v|<3/5.
    # We test parallelism with either target sign, so the antipodal change is safe.
    counts = dict(processed=0, splits=0, outside_tube=0, excluded_equation=0, maximum_depth=0)
    stack = [(chart, '', [Interval(-1, 1), Interval(Q(-3, 5), Q(3, 5)),
                          Interval(Q(-3, 5), Q(3, 5))]) for chart in [0, 1]]
    leaves = []
    while stack:
        chart, code, box = stack.pop()
        t, u, v = box
        depth = len(code)//2
        counts['processed'] += 1
        counts['maximum_depth'] = max(depth, counts['maximum_depth'])
        tube_gap = (u**2 + v**2).lo - K2*(1 + (t**2).hi)
        if tube_gap > 0:
            reason = 'outside_tube'
            witness = {'positive_gap': str(tube_gap)}
        else:
            h = raw([Interval(1), t, u, v] if chart == 0 else [t, Interval(1), u, v])
            equations = [h[i] + target[j]*h[1] for j, i in enumerate([0, 2, 3])]
            indices = [i for i, value in enumerate(equations) if not value.contains_zero()]
            if indices:
                reason = 'excluded_equation'
                witness = {'equation': indices[0], 'interval': equations[indices[0]].pair()}
            else:
                if depth >= maximum_depth:
                    raise RuntimeError(f'Unresolved box: chart {chart}, address {code}')
                counts['splits'] += 1
                coordinate = max(range(3), key=lambda j: (box[j].hi-box[j].lo)/(2 if j == 0 else Q(6, 5)))
                middle = (box[coordinate].lo + box[coordinate].hi)/2
                for bit in [1, 0]:
                    child = box.copy()
                    child[coordinate] = (Interval(middle, box[coordinate].hi) if bit
                                         else Interval(box[coordinate].lo, middle))
                    stack.append((chart, code+str(coordinate)+str(bit), child))
                continue
        counts[reason] += 1
        leaves.append({'chart': chart, 'address': code, 'reason': reason, **witness})
    assert len(leaves) == counts['splits'] + 2
    assert counts['processed'] == 2*counts['splits'] + 2
    return {'counts': counts, 'unresolved_boxes': 0, 'leaves': leaves}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        parser.error('Use a new output directory outside the repository, normally on D:.')
    output.mkdir(parents=True, exist_ok=False)
    target = [Interval(value-RADIUS, value+RADIUS) for value in CENTER]
    report = {
        'epsilon': str(E), 'tube_squared_ratio': str(K2),
        'previous_raw_vector': list(map(str, PREVIOUS)),
        'target_chart': '(a,-1,c,d)', 'target_box': [v.pair() for v in target],
        'root_certificate': certify_root(target),
        'predecessor_exclusion': exclude_predecessor(target),
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__, 'status': 'passed',
        'scope': 'A certified collision in T whose extra preimage is outside F(T), hence outside K. Neither injectivity nor noninjectivity of F restricted to K is established.',
    }
    (output/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status': report['status'], 'contraction': report['root_certificate']['contraction_decimal_for_display'],
                      'exclusion': report['predecessor_exclusion']['counts']}))


if __name__ == '__main__':
    main()
