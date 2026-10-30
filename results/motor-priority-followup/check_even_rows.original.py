"""Exact polynomial checks for the equal-cluster replication corollary.

Run with the repository's proof environment. This does not integrate a trajectory
or revalidate a Floquet certificate; it checks the algebra used to transfer one.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import flint
from flint import fmpq_mpoly_ctx


def main():
    ctx = fmpq_mpoly_ctx.get(
        ("h", "Y", "C", "r", "rho", "dc", "ds", "dW", "eps", "d", "e")
    )
    h, Y, C, r, rho, dc, ds, dW, eps, d, e = ctx.gens()
    core = [(-rho + r * C) * h - r * Y, h * C - Y, 1 - C - h * Y]
    checks = {}

    def require_zero(name, expression):
        checks[name] = expression == 0
        if not checks[name]:
            raise AssertionError(f"Nonzero residual for {name}: {expression}")

    def differentiate_core(expression):
        return sum(expression.derivative(i) * core[i] for i in range(3))

    # W = rho*X avoids rational functions; rho is constant along each orbit.
    a = [C + Y, C - Y, h - r * Y]
    b = [C - Y, C + Y, -h + r * Y]
    mean_s = (a[1] + b[1]) / 2
    require_zero("equal_cluster_mean_s_is_C", mean_s - C)
    require_zero("mechanical_constraint", a[2] + b[2])

    for label, row, velocity in [("A", a, -h), ("B", b, h)]:
        c, s, W = row
        V = -W + r * (s - mean_s)
        require_zero(f"{label}_velocity", V - velocity)
        field = [1 - c - V * s, 1 - s + V * c, rho * V]
        for k, coordinate in enumerate(row):
            require_zero(f"{label}_embedding_{k}", differentiate_core(coordinate) - field[k])

        # A zero-sum within-cluster contrast has zero variation of mean_s.
        cp, sp, Wp = c + eps * dc, s + eps * ds, W + eps * dW
        Vp = -Wp + r * (sp - mean_s)
        perturbed = [1 - cp - Vp * sp, 1 - sp + Vp * cp, rho * Vp]
        expected = [
            -dc + (-V - r * s) * ds + s * dW,
            V * dc + (-1 + r * c) * ds - c * dW,
            rho * r * ds - rho * dW,
        ]
        for k in range(3):
            linearized = perturbed[k].derivative("eps").subs({"eps": 0})
            require_zero(f"{label}_contrast_{k}", linearized - expected[k])
        require_zero(f"{label}_physical_radius", c * c + s * s - 2 * (C * C + Y * Y))

    # u=Y+d, v=C+e in the full five-state cluster restriction.
    d_prime = -(Y + d) + h * (C + e) - core[1]
    e_prime = -(C + e) - h * (Y + d) + 1 - core[2]
    require_zero("passive_d_equation", d_prime - (-d + h * e))
    require_zero("passive_e_equation", e_prime - (-e - h * d))
    require_zero("passive_norm_decay", 2 * d * d_prime + 2 * e * e_prime + 2 * (d * d + e * e))
    require_zero("core_radius_identity", 2 * C * core[2] + 2 * Y * core[1] - 2 * C + 2 * (C * C + Y * Y))

    # Counts are consequences of the archived target certificate, not new
    # eigenvalue computations. The argument for all even N is in the note.
    counts = []
    for N in [2, 4, 6, 8, 10, 16, 100]:
        m = N // 2
        row = {
            "N": N,
            "cluster_size": m,
            "cluster_dimension": 5,
            "contrast_blocks": 2 * (m - 1),
            "full_dimension": 3 * N - 1,
            "target_unstable": N - 2,
            "target_neutral": 1,
            "target_stable": 2 * N,
        }
        assert 5 + 3 * row["contrast_blocks"] == row["full_dimension"]
        assert row["target_unstable"] + row["target_neutral"] + row["target_stable"] == row["full_dimension"]
        counts.append(row)

    here = Path(__file__).resolve()
    root = here.parents[2]
    certificate = root / "results/motor-proof/target-recheck.json"
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Exact multivariate rational polynomial identities using FLINT; no numerical trajectories",
        "python_flint_version": flint.__version__,
        "script_sha256": hashlib.sha256(here.read_bytes()).hexdigest(),
        "inherited_target_certificate": "results/motor-proof/target-recheck.json",
        "inherited_target_certificate_sha256": hashlib.sha256(certificate.read_bytes()).hexdigest(),
        "identities": checks,
        "all_identities_zero": all(checks.values()),
        "illustrative_dimension_counts": counts,
        "limitations": [
            "Does not rerun interval orbit existence or Floquet isolation.",
            "All-even-N result also uses the direct-sum argument in docs/motor-followup-priority.md.",
            "No proof of chaos, scientific priority, or arbitrary biological row count."
        ],
    }
    output = here.with_name("even-rows-algebra.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"identities": len(checks), "all_zero": result["all_identities_zero"], "output": str(output)}))


if __name__ == "__main__":
    main()
