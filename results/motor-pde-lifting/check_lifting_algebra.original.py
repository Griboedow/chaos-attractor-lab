"""Exact sign and normalization controls for the motor density lifting note.

These polynomial identities supplement the explicit function-space argument;
they do not simulate a PDE or prove a Floquet classification by discretization.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import flint
from flint import fmpq_mpoly_ctx


def main():
    ctx = fmpq_mpoly_ctx.get(
        ("C", "S", "V", "n", "a", "b", "m", "A", "c", "s", "B", "gamma", "rho", "X", "sbar")
    )
    C, S, V, n, a, b, m, A, c, s, B, gamma, rho, X, sbar = ctx.gens()
    Cdot = -C - n * V * S + a
    Sdot = -S + n * V * C + b
    checks = {}

    def zero(name, expression):
        checks[name] = expression == 0
        if not checks[name]:
            raise AssertionError(f"Nonzero {name}: {expression}")

    zero("cosine_coefficient_PDE_balance", Cdot + n * V * S + C - a)
    zero("sine_coefficient_PDE_balance", Sdot - n * V * C + S - b)
    zero("complex_plus_i_convention_real", Cdot - (-C - n * V * S + a))
    zero("complex_plus_i_convention_imaginary", Sdot - (-S + n * V * C + b))
    zero("forced_amplitude_identity", 2 * C * Cdot + 2 * S * Sdot + 2 * (C*C + S*S) - 2 * (a*C + b*S))
    zero("unforced_amplitude_decay", (2 * C * Cdot + 2 * S * Sdot + 2 * (C*C + S*S)).subs({"a": 0, "b": 0}))
    zero("moving_phase_real", (Cdot + n * V * S + C).subs({"a": 0}))
    zero("moving_phase_imaginary", (Sdot - n * V * C + S).subs({"b": 0}))
    mdot = ctx.constant(1) / 2 - m
    zero("mean_relaxation", mdot + m - ctx.constant(1) / 2)
    zero("mean_deviation_squared_decay", 2 * (m - ctx.constant(1) / 2) * mdot + 2 * (m - ctx.constant(1) / 2)**2)
    normalized_velocity = -rho * X - gamma * A * (s - sbar)
    raw_cdot = -A*c - normalized_velocity*A*s + A
    raw_sdot = -A*s + normalized_velocity*A*c + B
    zero("normalized_first_cosine", raw_cdot - A * (-c - normalized_velocity*s + 1))
    # Multiply by A so B/A does not require a rational-function context.
    zero("normalized_first_sine", raw_sdot - (A * (-s + normalized_velocity*c) + B))
    zero("mechanical_mean_no_feedback", normalized_velocity.derivative("m"))

    here = Path(__file__).resolve()
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Exact rational multivariate polynomial identities in Python-FLINT",
        "python_flint_version": flint.__version__,
        "script_sha256": hashlib.sha256(here.read_bytes()).hexdigest(),
        "identities": checks,
        "all_zero": all(checks.values()),
        "scope": "Checks signs, mean normalization, amplitude decay, moving-phase cancellation and first-harmonic scaling",
        "not_computed": ["No PDE truncation", "No trajectories", "No new interval orbit or spectrum certificate"],
        "analytical_argument": "docs/motor-pde-lifting.md"
    }
    output = here.with_name("algebra.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"identities": len(checks), "all_zero": result["all_zero"]}))


if __name__ == "__main__":
    main()
