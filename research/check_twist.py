"""Independent algebra checks and arbitrary-precision finite-separation tests."""
import argparse
import json
from pathlib import Path
import time
import mpmath as mp
import numpy as np
from twist_map import step, step_tangent, spectra


def independent(z, p, cos=np.cos, sin=np.sin):
    """Explicit scalar composition; does not call the compiled evaluator."""
    r, a, b, c, k, l, m, u, v, w = p
    x, y, z = r * z[0] + 1 - r, r * z[1], r * z[2]
    t = a + k * z + u * z ** 3
    x, y = cos(t) * x - sin(t) * y, sin(t) * x + cos(t) * y
    t = b + l * x + v * x ** 3
    y, z = cos(t) * y - sin(t) * z, sin(t) * y + cos(t) * z
    t = c + m * y + w * y ** 3
    z, x = cos(t) * z - sin(t) * x, sin(t) * z + cos(t) * x
    return [x, y, z]


def inverse(z, p):
    z = z.copy()
    for s in [2, 1, 0]:
        i, j, k = s, (s + 1) % 3, (s + 2) % 3
        t = -(p[1+s] + p[4+s] * z[k] + p[7+s] * z[k] ** 3)
        z[i], z[j] = np.cos(t)*z[i]-np.sin(t)*z[j], np.sin(t)*z[i]+np.cos(t)*z[j]
    z[0] -= 1-p[0]
    return z / p[0]


def algebra(parameters):
    rng = np.random.default_rng(20261106)
    errors = dict(state=0., jacobian=0., determinant=0., inverse=0., bound=0.)
    for p in parameters:
        for z in rng.uniform(-1, 1, (256, 3)):
            value, jac = step_tangent(z, np.eye(3), p)
            cs = np.column_stack([np.imag(independent(z.astype(complex) + 1e-30j*np.eye(3)[i], p))/1e-30 for i in range(3)])
            errors['state'] = max(errors['state'], float(np.max(np.abs(value-independent(z,p)))))
            errors['jacobian'] = max(errors['jacobian'], float(np.max(np.abs(jac-cs))))
            errors['determinant'] = max(errors['determinant'], abs(float(np.linalg.det(jac)-p[0]**3)))
            errors['inverse'] = max(errors['inverse'], float(np.max(np.abs(inverse(value,p)-z))))
            errors['bound'] = max(errors['bound'], float(np.linalg.norm(value)-p[0]*np.linalg.norm(z)-(1-p[0])))
    assert errors['state'] < 1e-13 and errors['jacobian'] < 1e-11
    assert errors['determinant'] < 1e-11 and errors['inverse'] < 1e-10
    assert errors['bound'] < 1e-13
    control = np.array(parameters)
    control[:,4:] = 0
    result = spectra(control, np.tile([.1,.2,.3], (len(control),1)), 2000, 4000, 1000, 16)[0]
    error = float(np.max(np.abs(result - np.log(control[:,0,None,None]))))
    assert error < 1e-11
    errors['zero_twist_spectrum'] = error
    return errors


def precision_run(p, digits, warmup=10000, observed=30000):
    mp.mp.dps = digits
    # Exact promotion of the archived binary64 parameters, not new rounded values.
    p = list(map(mp.mpf, map(float, p)))
    z = list(map(mp.mpf, ['.1','.2','.3']))
    eps = mp.mpf('1e-30')
    direction = [mp.mpf(1),mp.mpf(0),mp.mpf(0)]
    logs = mp.mpf(0)
    windows = []
    mean = [mp.mpf(0)]*3
    for i in range(warmup + observed):
        plus = independent([z[j]+eps*direction[j] for j in range(3)],p,mp.cos,mp.sin)
        value = independent(z,p,mp.cos,mp.sin)
        delta = [plus[j]-value[j] for j in range(3)]
        distance = mp.sqrt(mp.fsum(x*x for x in delta))
        direction = [x/distance for x in delta]
        z = value
        if i >= warmup:
            logs += mp.log(distance/eps)
            for j in range(3):
                mean[j] += z[j]/observed
            if (i-warmup+1) % 5000 == 0:
                windows.append(float(logs/5000));logs=mp.mpf(0)
    return {'digits':digits,'separation':'1e-30','warmup':warmup,'observed':observed,'window':5000,'windows':windows,'mean_exponent':float(np.mean(windows)),'state_mean':list(map(float,mean))}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--directory',type=Path,default=Path('results/twist-geometry'))
    ap.add_argument('--precision',action='store_true')
    args=ap.parse_args()
    candidates=json.loads((args.directory/'candidates.json').read_text())['candidates']
    result={'algebra':algebra([c['parameters'] for c in candidates])}
    print(json.dumps(result),flush=True)
    if args.precision:
        result['precision']=[]
        for c in candidates[:2]:
            for digits in [70,90]:
                row={'id':c['id'],**precision_run(c['parameters'],digits)}
                result['precision'].append(row)
                print(json.dumps(row),flush=True)
    (args.directory/'independent-checks.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
