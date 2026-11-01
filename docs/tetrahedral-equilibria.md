# Certified equilibria near the tetrahedral candidates

The [tetrahedral Gaussian-vortex flow](tetrahedral-flow.md) has more structure than a stable center surrounded by a chaotic-looking cloud. Each of cases **421 and 565 has at least 33 distinct equilibria**, certified in local interval boxes. One is stable, sixteen have one unstable direction, and sixteen have two. This is a local mathematical result for the stored parameters; it does not establish chaotic attraction or historical novelty.

The initial inventory used 1,024 scrambled Sobol starts in `[-3,3]^3`, additional roots on two symmetry axes, an analytic Jacobian and completion under the twelve proper tetrahedral rotations. It is **not an exhaustive enumeration**. The accepted numerical roots form five radial groups:

| Case | Radius, approximately | Count | Unstable directions per root |
| --- | ---: | ---: | ---: |
| 421 | 0 | 1 | 0 |
| 421 | 0.1862144 | 4 | 1 |
| 421 | 0.6212786 | 12 | 2 |
| 421 | 1.2255702 | 12 | 1 |
| 421 | 1.2443980 | 4 | 2 |
| 565 | 0 | 1 | 0 |
| 565 | 0.0503526 | 4 | 1 |
| 565 | 0.9527685 | 12 | 2 |
| 565 | 1.5194316 | 12 | 1 |
| 565 | 1.5768110 | 4 | 2 |

For each numerical center `x0`, a box of coordinate radius `1e-8` encloses a unique root. The verifier evaluates the vector field and Jacobian with outward `mpmath.iv` arithmetic, using the exact integer ratios of every stored binary64 parameter. With a fixed approximate inverse `C`, it checks strict inclusion of

`K(X) = x0 - C f(x0) + (I - C J(X))(X - x0)`

inside `X`, and an interval infinity-norm bound below one for `I-CJ(X)`. The actual interval subtraction `X-x0` is used, including its rounding enclosure. The boxes are pairwise disjoint. Both 60- and 90-decimal runs pass all 66 roots.

The characteristic polynomial is `lambda^3 + a1 lambda^2 + a2 lambda + a3`. Strict interval signs of its Routh first column, `1, a1, (a1*a2-a3)/a1, a3`, give the number of eigenvalues in the open right half-plane. No sign is ambiguous in any box. Thus the unstable-direction counts are certified rather than inferred from floating-point eigenvalues alone. As with other computational certificates, the conclusion relies on the stated interval implementation and runtime; it is not a proof-assistant formalization.

The orbit passes near several outer saddles. That suggests objects worth investigating, but proximity does not prove a heteroclinic connection or a chaos mechanism. In particular, the stable origin alone does **not** establish a formally hidden attractor: the basins near the other unstable equilibria remain to be examined.

## Reproduction

The [packet](../results/tetrahedral-equilibria/README.md) records the inventory, interval bounds, exact original scripts and input hashes. Portable scripts require an explicit output directory and refuse to replace existing certificates. From the repository root, with NumPy, SciPy, Numba and mpmath installed:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:CHAOS_SEARCH_OUTPUT='D:/Codex/reproduce-tetra-equilibria/cache'
python research/tetrahedral_equilibria/inventory.py --output D:/Codex/reproduce-tetra-equilibria/inventory
python research/tetrahedral_equilibria/certify.py --precision 60 --inventory D:/Codex/reproduce-tetra-equilibria/inventory/inventory.json --output D:/Codex/reproduce-tetra-equilibria/certificates
python research/tetrahedral_equilibria/certify.py --precision 90 --inventory D:/Codex/reproduce-tetra-equilibria/inventory/inventory.json --output D:/Codex/reproduce-tetra-equilibria/certificates
```

To recertify the published centers directly, omit `--inventory`. The numerical root finder may return slightly different centers across library versions; the interval verifier must still pass. The portable verifier strengthens the original floating-point box-separation check to an outward interval comparison. It otherwise evaluates the same certificate.
