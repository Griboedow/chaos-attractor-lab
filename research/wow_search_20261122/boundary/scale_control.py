"""Test the local phase permutation by directly iterating concentric pixel pairs in Arb."""
import json
import numpy as np
from flint import acb, arb, ctx
from certificate import configure, nmap
from search import RAW, OUTPUT
from selected import D, K


def capture(z, parameter, root_centers):
    for step in range(900):
        for root, center in enumerate(root_centers):
            if abs(z - center) < arb('1e-7'):
                return root, step
        if abs(z) < arb('1e-7'):
            return D + (-step) % K, step
        z = nmap(z, parameter)
    return -1, -1


def main():
    ctx.prec = 656
    _, parameter, _, roots, _ = configure()
    pole = (parameter / D).root(D - 1)
    with np.load(RAW / 'case-6902-microscopic.npz') as src:
        labels = src['labels']
    ny, nx = labels.shape
    rng = np.random.default_rng(22006902)
    sample = list(rng.choice(nx * ny, 144, replace=False))
    for label in range(D + K):
        sample.extend(rng.choice(np.flatnonzero(labels == label), 4, replace=False))
    sample = np.unique(sample)
    c = arb(D - 1) / D
    tests = []
    for scale in ['1e-1', '1e-3', '1e-8', '1e-30']:
        radius = arb(scale)
        counts = np.zeros((D + K, D + K), dtype=int)
        extra_step = 0
        unknown = 0
        for point in sample:
            y, x = divmod(int(point), nx)
            offset = acb(arb(2*x + 1)/nx - 1, arb(4)/5 - arb(4*(2*y+1))/(5*ny))
            source, n = capture(pole + radius * offset, parameter, roots)
            target, m = capture(pole + c * radius * offset, parameter, roots)
            if min(source, target) < 0:
                unknown += 1
            else:
                counts[source, target] += 1
                extra_step += int(m == n + 1)
        expected = [j if j < D else D + (j-D-1) % K for j in range(D+K)]
        matched = [int(counts[j, t]) for j, t in enumerate(expected)]
        mismatch = [[int(j), int(t), int(counts[j,t])] for j,t in zip(*np.nonzero(counts)) if t != expected[j]]
        tests.append(dict(half_width=scale, classified_pairs=int(counts.sum()), unknown=unknown,
                          permutation_agreements=sum(matched), matched_by_source_label=matched,
                          one_extra_capture_iteration=extra_step, off_permutation=mismatch))
        print(scale, sum(matched), int(counts.sum()), extra_step, flush=True)
    result = dict(case=6902, precision_bits=ctx.prec, sample_size=len(sample), seed=22006902,
                  sampling='144 uniform pixel centers plus four chosen per phase in the finest view, deduplicated; identical offsets in all four windows.',
                  dilation='5/6 around the same exact pole in the native complex coordinate',
                  expected_permutation=expected, tests=tests)
    (OUTPUT / 'scale-control.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
