"""Two ensembles with identical pixel counts and different histories under Arnold's cat map."""

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


MODULUS = 2**31 - 1
SIDE = 320
MASK_SIDE = 1280
LAG = 15
COUNT = 2**18
SEEDS = {"conditioned": 2026112201, "twin": 2026112202, "control": 2026112203}
OUTPUT = Path(__file__).resolve().parents[2] / "results/art_20261122/double-exposure"


def masks():
    # Boundaries of the fixed-radius circles in the archived complex-plane observation.
    u = np.linspace(0, 2 * np.pi, 12001)
    center = np.exp(3j * u) + (2 / 3) * np.exp(-4j * u)
    derivative = 3j * np.exp(3j * u) - (8j / 3) * np.exp(-4j * u)
    normal = 1j * derivative / np.abs(derivative)
    image = Image.new("1", (MASK_SIDE, MASK_SIDE))
    draw = ImageDraw.Draw(image)
    for sign in (-1, 1):
        curve = (center + sign * 0.27 * normal) / 2.12
        xy = np.column_stack([1 + curve.real, 1 - curve.imag])
        draw.line(list(map(tuple, np.rint(MASK_SIDE * xy / 2).astype(int))),
                  fill=1, width=13)

    # An almond, an iris with radial striae, and a pupil with a single catchlight.
    coord = 2 * (np.arange(MASK_SIDE) + 0.5) / MASK_SIDE - 1
    x, y = coord[None, :], -coord[:, None]
    radius = np.hypot(x, y)
    lid = 0.43 * (1 - (x / 0.87)**2)
    outline = (abs(x) < 0.88) & (abs(abs(y) - lid) < 0.021)
    upper_fold = (abs(x) < 0.70) & (abs(y - lid - 0.085) < 0.011)
    iris = (abs(radius - 0.315) < 0.018) | (abs(radius - 0.186) < 0.010)
    striae = (radius > 0.199) & (radius < 0.29)
    striae &= np.cos(32 * np.arctan2(y, x) + 2 * np.sin(9 * radius)) > 0.69
    pupil = (radius < 0.145) & (np.hypot(x + 0.049, y - 0.060) > 0.047)
    return np.asarray(image, dtype=bool), outline | upper_fold | iris | striae | pupil


@lru_cache(maxsize=None)
def power(steps):
    matrix = [[2, 1], [1, 1]] if steps >= 0 else [[1, -1], [-1, 2]]
    return np.linalg.matrix_power(np.array(matrix, dtype=np.int64), abs(steps))


def evolve(x, y, steps):
    (a, b), (c, d) = power(steps)
    # For |steps| <= 15, both dot products fit int64 on the entire grid.
    if abs(steps) > LAG:
        raise ValueError("Propagate longer intervals in separate legs.")
    return (a * x + b * y) % MODULUS, (c * x + d * y) % MODULUS


def in_mask(mask, x, y):
    return mask[y * MASK_SIDE // MODULUS, x * MASK_SIDE // MODULUS]


def sample_centers(past_mask, future_mask):
    rng = np.random.default_rng(SEEDS["conditioned"])
    samples_x, samples_y = [], []
    remaining = COUNT
    while remaining:
        x, y = rng.integers(MODULUS, size=(2, 2**18), dtype=np.int64)
        past = in_mask(past_mask, *evolve(x, y, -LAG))
        x, y = x[past], y[past]
        future = in_mask(future_mask, *evolve(x, y, LAG))
        x, y = x[future][:remaining], y[future][:remaining]
        samples_x.append(x)
        samples_y.append(y)
        remaining -= len(x)
    return np.concatenate(samples_x), np.concatenate(samples_y)


def hidden_twin(x, y):
    rng = np.random.default_rng(SEEDS["twin"])

    def resample(values):
        cells = values * SIDE // MODULUS
        low = (cells * MODULUS + SIDE - 1) // SIDE
        high = ((cells + 1) * MODULUS + SIDE - 1) // SIDE
        return rng.integers(low, high, dtype=np.int64)

    return resample(x), resample(y)


def histogram(x, y, side):
    keys = (y * side // MODULUS) * side + x * side // MODULUS
    return np.bincount(keys, minlength=side**2).reshape(side, side)


def uniformity(x, y, side):
    counts = histogram(x, y, side)
    expected = COUNT / side**2
    return float(np.sum((counts - expected)**2 / expected) / (side**2 - 1))


def typeface(style, size):
    # Use the redistributable DejaVu faces that ship with Matplotlib.
    from matplotlib import get_data_path

    return ImageFont.truetype(str(Path(get_data_path(), "fonts", "ttf", style + ".ttf")), size)


def render(rows, filename="two-histories.png"):
    margin, gap, top, label_space = 24, 15, 53, 53
    plate = Image.new("RGB", (2 * margin + 3 * SIDE + 2 * gap,
                              top + 2 * (SIDE + label_space) + 8), "#101921")
    draw = ImageDraw.Draw(plate)
    column = typeface("DejaVuSansMono", 17)
    label = typeface("DejaVuSans", 18)
    for i, title in enumerate(("PAST   -15", "PRESENT   0", "FUTURE   +15")):
        draw.text((margin + i * (SIDE + gap), 19), title, font=column, fill="#f3eee6")
    colors = np.array([[5, 16, 27], [18, 34, 48], [52, 73, 91],
                       [144, 155, 157], [233, 195, 156], [255, 238, 212]])
    knots = [0, 1, 3, 7, 16, 33]
    scale = np.arange(256)
    palette = np.column_stack([np.interp(scale, knots, channel) for channel in colors.T])
    for r, (title, panels) in enumerate(rows):
        y = top + (SIDE + label_space) * r
        for c, counts in enumerate(panels):
            rgb = palette[np.minimum(counts, 255)].astype(np.uint8)
            plate.paste(Image.fromarray(rgb), (margin + c * (SIDE + gap), y))
        draw.text((margin, y + SIDE + 9), title, font=label, fill="#e3ded6")
    plate.save(OUTPUT / filename, optimize=True)


def main():
    past_mask, future_mask = masks()
    centers = sample_centers(past_mask, future_mask)
    twin = hidden_twin(*centers)
    rows, observations = [], {}
    for name, points, title in [
        ("conditioned", centers, "Original points"),
        ("twin", twin, "Same present pixels, new subpixel positions"),
    ]:
        past, future = evolve(*points, -LAG), evolve(*points, LAG)
        for endpoint, steps in [(past, LAG), (future, -LAG)]:
            recovered = evolve(*endpoint, steps)
            assert all(np.array_equal(got, want) for got, want in zip(recovered, points))
        hits_past = in_mask(past_mask, *past)
        hits_future = in_mask(future_mask, *future)
        observations[name] = {"past_hits": int(hits_past.sum()), "future_hits": int(hits_future.sum()),
                              "both_hits": int((hits_past & hits_future).sum())}
        rows.append((title, [histogram(*state, SIDE) for state in (past, points, future)]))
    mismatches = int(np.count_nonzero(rows[0][1][1] != rows[1][1][1]))
    assert mismatches == 0 and observations["conditioned"]["both_hits"] == COUNT
    control = np.random.default_rng(SEEDS["control"]).integers(MODULUS, size=(2, COUNT))
    summary = {
        "points": COUNT, "prime_grid": MODULUS, "visible_grid": SIDE,
        "mask_grid": MASK_SIDE, "steps_to_each_endpoint": LAG, "pcg64_seeds": SEEDS,
        "past_mask_area_fraction": float(past_mask.mean()),
        "future_mask_area_fraction": float(future_mask.mean()),
        "different_present_pixel_counts": mismatches, "endpoint_counts": observations,
        "present_pearson_chi_square_per_degree_of_freedom": {
            str(side): {"conditioned": uniformity(*centers, side),
                        "uniform_control": uniformity(*control, side)} for side in (32, 80, 160, SIDE)},
    }
    OUTPUT.mkdir(exist_ok=True, parents=True)
    (OUTPUT / "measurements.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    render(rows)


if __name__ == "__main__":
    main()
