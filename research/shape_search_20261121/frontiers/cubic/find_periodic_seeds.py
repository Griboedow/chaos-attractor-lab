"""Find one-turn orbit proposals for Case 48000001. NumPy and SciPy; no proof arithmetic."""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.linalg import block_diag
from scipy.ndimage import gaussian_filter1d
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

PERIOD = 2 * np.pi
GAIN = 0.34
NODES = 64
ROOT = Path(__file__).resolve().parents[4]
GRADIENT = ROOT / "results/shape_search_20261121/frontiers/cubic/gradient.json"
DEFAULT_SCOUTS = [(17, 96, 600), (27, 144, 600)]  # seed, midpoint mesh, total starts


def drive(t):
    return 0.25 * np.sin(t) - 4 * np.sin(3 * t) + 3j * np.sin(6 * t)


def velocity(z, b):
    return -GAIN * np.conj(3 * z * z - 3) * (z**3 - 3 * z + b)


def jacobian(z, b):
    a = -GAIN * abs(3 * z * z - 3) ** 2
    c = -6 * GAIN * z.conj() * (z**3 - 3 * z + b)
    return np.array([[a + c.real, c.imag], [c.imag, a - c.real]])


def periodic_midpoints(z, b, step, indices):
    def residual(curve):
        following = np.roll(curve, -1)
        return following - curve - step * velocity((following + curve) / 2, b)

    error = residual(z)
    for _ in range(55):
        if np.max(abs(error)) < 1e-10:
            return z
        middle = (z + np.roll(z, -1)) / 2
        a, bxy, byx, d = (-step / 2 * jacobian(middle, b)).reshape(4, -1)
        entries = np.array([a - 1, bxy, byx, d - 1, a + 1, bxy, byx, d + 1])
        matrix = coo_matrix(
            (entries.ravel(), indices), shape=(2 * len(z), 2 * len(z))
        ).tocsc()
        delta = spsolve(matrix, -error.view(float)).view(complex)
        if not np.all(np.isfinite(delta)):
            return None
        previous = np.linalg.norm(error)
        for j in range(12):
            trial = z + 0.5**j * delta
            error = residual(trial)
            if np.linalg.norm(error) < previous:
                z = trial
                break
        else:
            return None
    return None


def different(curve, found):
    return all(np.max(abs(curve - old)) > 1e-4 for old in found)


def scout(job):
    (seed, mesh, attempts), starts = job
    step = PERIOD / mesh
    phase = np.arange(mesh) * step
    b = drive(phase + step / 2)
    # Regenerate the three sinks from the recorded phase-zero states.
    sol = solve_ivp(
        lambda t, z: velocity(z, drive(t)),
        (0, PERIOD),
        starts,
        t_eval=phase,
        method="DOP853",
        rtol=2e-12,
        atol=2e-14,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    basis = np.vstack([sol.y, np.ones(mesh), -np.ones(mesh), np.zeros(mesh)])
    i = np.arange(mesh) * 2
    j = np.roll(i, -1)
    indices = (
        np.array([i, i, i + 1, i + 1, i, i, i + 1, i + 1]).ravel(),
        np.array([i, i + 1, i, i + 1, j, j + 1, j, j + 1]).ravel(),
    )
    rng = np.random.default_rng(seed)
    found = []
    converged = 0
    for trial in range(attempts):
        if trial < len(basis):
            initial = basis[trial]
        else:
            pieces = rng.choice([2, 3, 4, 6, 6, 6, 9, 12])
            offset = rng.integers(mesh)
            assign = ((np.arange(mesh) + offset) % mesh) * pieces // mesh
            levels = rng.integers(len(basis), size=pieces)
            stitched = basis[levels[assign], np.arange(mesh)]
            width = rng.uniform(0.12, 0.7) * mesh / pieces / 3
            initial = gaussian_filter1d(stitched, width, mode="wrap")
            initial += rng.uniform(0, 0.4) * np.exp(1j * rng.uniform(-3.14, 3.14))
        curve = periodic_midpoints(initial, b, step, indices)
        if curve is not None:
            converged += 1
            if different(curve, found):
                found.append(curve)
    print(
        f"Seed {seed}, mesh {mesh}: {converged}/{attempts} converged, {len(found)} midpoint curves.",
        flush=True,
    )
    return found, converged


def refine(curve):
    def rhs(t, y):
        f = velocity(y[0] + 1j * y[1], drive(t))
        return np.array([f.real, f.imag])

    def derivative(t, y):
        return jacobian(y[0] + 1j * y[1], drive(t))

    closed = np.r_[curve, curve[0]]
    mesh = np.linspace(0, PERIOD, len(closed))
    sol = solve_bvp(
        rhs,
        lambda a, b: a - b,
        mesh,
        np.array([closed.real, closed.imag]),
        fun_jac=derivative,
        tol=2e-8,
        max_nodes=12000,
    )
    if not sol.success:
        return None
    values = sol.sol(np.arange(4 * NODES) * PERIOD / (4 * NODES))
    return values[0] + 1j * values[1]


def polish_nodes(orbit):
    """Newton correction of the 64 short flow segments before interval verification."""
    nodes = np.ascontiguousarray(orbit[::4])
    phases = np.arange(NODES) * PERIOD / NODES
    starts = np.tile(np.eye(2).ravel(), (NODES, 1))
    shift = np.roll(np.eye(2 * NODES), 2, axis=1)

    def rhs(t, state):
        q = state.reshape(NODES, 6)
        z = q[:, 0] + 1j * q[:, 1]
        b = drive(phases + t)
        f = velocity(z, b)
        matrices = jacobian(z, b).transpose(2, 0, 1) @ q[:, 2:].reshape(-1, 2, 2)
        return np.column_stack((f.real, f.imag, matrices.reshape(-1, 4))).ravel()

    for _ in range(4):
        initial = np.column_stack((nodes.real, nodes.imag, starts)).ravel()
        sol = solve_ivp(
            rhs, (0, PERIOD / NODES), initial, method="DOP853", rtol=3e-13, atol=1e-14
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        end = sol.y[:, -1].reshape(NODES, 6)
        error = np.roll(nodes, -1) - (end[:, 0] + 1j * end[:, 1])
        derivative = shift - block_diag(*end[:, 2:].reshape(-1, 2, 2))
        correction = np.linalg.solve(derivative, error.view(float)).view(complex)
        nodes -= correction
        if np.max(abs(correction)) < 1e-13:
            return nodes
    raise ArithmeticError("Multiple shooting did not resolve the periodic nodes")


def main(args):
    if args.out.exists():
        raise FileExistsError(f"Choose a new output file: {args.out}")
    settings = args.scout or DEFAULT_SCOUTS
    if (
        any(n < 4 or count <= 0 or seed < 0 for seed, n, count in settings)
        or args.workers < 1
    ):
        raise ValueError(
            "Supply a nonnegative seed, at least four mesh nodes, positive starts and workers."
        )
    case = next(
        row for row in json.loads(args.gradient.read_text()) if row["case"] == 48000001
    )
    starts = next(r["starts_last_period"] for r in case["runs"] if r["gain"] == GAIN)
    starts = np.array([complex(x, y) for x, y in starts])
    jobs = [(config, starts) for config in settings]
    if args.workers == 1:
        surveys = list(map(scout, jobs))
    else:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
            surveys = list(pool.map(scout, jobs))

    curves = []
    refined = 0
    for found, _ in surveys:
        for midpoint in found:
            orbit = refine(midpoint)
            if orbit is None:
                continue
            refined += 1
            # b(t+pi) = -conj(b(t)), so the half-turn reflection is another solution.
            for candidate in (orbit, -np.conj(np.roll(orbit, len(orbit) // 2))):
                if different(candidate, curves):
                    curves.append(candidate)
    if not curves:
        raise RuntimeError("No refined periodic proposals.")
    nodes = np.asarray([polish_nodes(curve) for curve in curves])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("xb") as f:
        np.savez(
            f, nodes=nodes, scouts=np.array(settings), converged=[k for _, k in surveys]
        )
    print(
        f"{refined} refined midpoint curves; {len(nodes)} distinct one-turn proposals at {NODES} phase nodes."
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out", required=True, type=Path, help="New NPZ file; no overwrites."
    )
    p.add_argument("--gradient", type=Path, default=GRADIENT)
    p.add_argument(
        "--scout",
        nargs=3,
        type=int,
        action="append",
        metavar=("SEED", "MESH", "STARTS"),
        help="Repeat for independent scouts. Defaults: 17 96 600 and 27 144 600. Each total includes six deterministic starts.",
    )
    p.add_argument(
        "--workers", type=int, default=1, help="Concurrent independent scouts."
    )
    main(p.parse_args())
