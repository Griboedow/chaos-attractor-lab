# FastICA preimage rejection

[Statement and proof](../../docs/fastica-preimage-rejection.md) ·
[Original construction](../../docs/fastica-tripling.md)

At exact epsilon `1/25`, [report.json](report.json) certifies a collision in
the trapping torus and excludes the entire box containing its extra preimage
from the image of that torus. The extra preimage is therefore outside the
attracting set. This is a rejected intrinsic-collision candidate, not an
intrinsic classification of the attractor.

The proof uses exact rational polynomial arithmetic. It verifies the explicit
cubic against the earlier cumulant polynomial, constructs a strict contracting
self-map of the root box, checks input membership and common-image sign, and
exhaustively excludes predecessor solutions in two charts covering the torus.
The 45 terminal box addresses and exact witnesses are retained in the report.

## Reproduce

Only Python and SymPy are required for the certificate. Use the
[original study environment](../fastica-tripling/environment.json) or a Python
environment with SymPy 1.14.0. From the repository root:

```powershell
python research/fastica_preimage_rejection.py --output D:/Codex/fastica-preimage-recheck
python research/verify_fastica_preimage_rejection.py
```

The output directory must be new and outside the checkout. The verifier checks
file hashes, independently audits the leaf-address partition, and replays the
rational certificate. The replay uses the same interval implementation;
it is not an independent arithmetic implementation or a formal proof assistant.

The frozen [check source](check-source.py) records the executed bytes. Run the
research entry point above; the snapshot retains its original relative-path
assumptions. [Manifest](manifest.json) lists the packet and dependencies.

[exploration.json](exploration.json) and [exploration-source.py](exploration-source.py)
retain the original 16-target, 448-solve proposal search. That archived script
uses its original D: paths and dynamic Numba loading; it is historical evidence,
not the portable reproduction entry point. Its approximate roots do not prove
existence, completeness or invariant-set membership. The new certificate is
self-contained and requires neither SciPy, Numba, trajectories nor a GPU.

The [search scope](search-scope.json) records the limited supplementary query
audit. No search result is used as a historical-uniqueness certificate.
