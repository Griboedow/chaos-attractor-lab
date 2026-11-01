# FastICA local-injectivity checks

[Full argument](../../docs/fastica-tripling-structure.md) ·
[Original construction](../../docs/fastica-tripling.md)

The [report](report.json) checks algebra used to prove that the full map
fails to be locally injective at some point of its attracting set for every
`0 < epsilon <= 1/25`. The proof combines those checks with connectedness,
invariance of domain and local orientations. It is not a machine-formalized
topological proof.

The retained computation contains:

- Two symbolic determinant factorizations and four exact endpoint identities.
- Ten positive rational margins, including periodic-point containment in the tube.
- Nine exact Bernstein positivity certificates covering the full parameter
  interval and the relevant `t/epsilon` rectangles; polynomial reconstruction
  is verified before accepting a coefficient bound.
- Eight supplementary 90-digit periodic-point controls at `epsilon=0.01`,
  `0.02`, `0.03`, `0.04`, with oriented tangent frames and two centered-difference
  scales. The largest derivative-entry discrepancy was below `5.3e-10`.

No sampled trajectory is used to infer membership of the special points in
the attracting set: they have exact periodic-orbit existence proofs inside
the trapping region. The reported decimal coordinates are illustrations of
those roots, not substitutes for their sign brackets.

## Reproduce

Use the [original study environment](../fastica-tripling/environment.json).
From the repository root:

```powershell
python research/fastica_local_injectivity.py --output D:/Codex/fastica-injectivity-recheck
```

The output directory must be new and outside the checkout. This short CPU
check requires NumPy, SymPy and mpmath; no trajectories or GPU are needed.
The executed [source snapshot](check-source.py), [dependency manifest](manifest.json)
and [primary-source review](source-review.json) preserve the scope of the
calculation. Run the current research entry point above; the source snapshot
retains its original repository-relative path assumptions.

The result excludes a conjugacy on a neighborhood of the entire attracting
set to a locally injective map. It does not prove noninjectivity of the
restriction to the invariant set, classify a Whitney fold, establish an
intrinsic solenoid obstruction or certify historical uniqueness.
