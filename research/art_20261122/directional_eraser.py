"""Move each particle along one lattice direction, inside its original pixel."""
import json
from pathlib import Path

import numpy as np

from double_exposure import LAG, MODULUS, SIDE, evolve, histogram, in_mask, masks, render, sample_centers


OUT = Path(__file__).resolve().parents[2] / 'results/art_20261122/double-exposure'
SEED = 20261123


def step_bounds(coord, stride):
    """Inclusive bounds on m such that coord + m*stride stays in its pixel."""
    pixel = coord * SIDE // MODULUS
    lo = (pixel * MODULUS + SIDE - 1) // SIDE
    hi = ((pixel + 1) * MODULUS - 1) // SIDE
    if stride < 0:
        return step_bounds_from_offsets(coord - hi, coord - lo, -stride)
    return step_bounds_from_offsets(lo - coord, hi - coord, stride)


def step_bounds_from_offsets(low, high, stride):
    return -((-low) // stride), high // stride


def move_in_pixel(x, y, direction, rng):
    dx, dy = direction
    lx, hx = step_bounds(x, dx)
    ly, hy = step_bounds(y, dy)
    amount = rng.integers(np.maximum(lx, ly), np.minimum(hx, hy) + 1)
    return x + amount * dx, y + amount * dy


def main():
    past_mask, future_mask = masks()
    x, y = sample_centers(past_mask, future_mask)
    # A^10 e1 and A^(-10) e2. They are perpendicular and have equal length.
    directions = {
        'forward_stretch': ((10946, 6765), 'Along the forward-stretching direction'),
        'backward_stretch': ((-6765, 10946), 'Along the backward-stretching direction'),
    }
    rng = np.random.default_rng(SEED)
    results = {'seed': SEED, 'particles': len(x), 'side': SIDE, 'lag': LAG, 'directions': {}}
    rows = []
    for name, (direction, caption) in directions.items():
        nx, ny = move_in_pixel(x, y, direction, rng)
        assert np.array_equal(x * SIDE // MODULUS, nx * SIDE // MODULUS)
        assert np.array_equal(y * SIDE // MODULUS, ny * SIDE // MODULUS)
        past_points = evolve(nx, ny, -LAG)
        future_points = evolve(nx, ny, LAG)
        past = in_mask(past_mask, *past_points)
        future = in_mask(future_mask, *future_points)
        rows.append((caption, [histogram(*state, SIDE) for state in (past_points, (nx, ny), future_points)]))
        results['directions'][name] = {
            'vector': direction, 'changed_center_pixels': 0,
            'unmoved_particles': int(np.count_nonzero(x == nx)),
            'past_hits': int(past.sum()), 'future_hits': int(future.sum()),
        }
    OUT.mkdir(exist_ok=True, parents=True)
    render(rows, filename='directional-eraser.png')
    (OUT / 'directional-eraser.json').write_text(json.dumps(results, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
