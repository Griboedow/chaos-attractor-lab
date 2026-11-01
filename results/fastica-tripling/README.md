# Structured FastICA evidence

[Construction and full proof](../../docs/fastica-tripling.md) ·
[Main-page portraits](../../README.md#structured-data-fastica)

The proof concerns the exact 720-row distribution and every
`0 <= epsilon <= 1/25`. The portrait at `epsilon=0.15` lies outside that
proved interval and has numerical evidence only. A factor onto the expanding
circle does not establish transitivity, a solenoid conjugacy or historical
novelty.

## Retained records

| Record | Scope |
| --- | --- |
| [Exact algebra](exact-check.json) | Symbolic moment identities and rational trapping inequalities |
| [Direct controls](controls.json) | 1,024 empirical-map comparisons and 16,384 circle angles |
| [Pilot protocol](pilot-01/protocol.json), [all outcomes](pilot-01/report.json) | All 80 parameters, including those failing the positive-growth gate |
| [Long protocol](long-01/protocol.json), [spectra](long-01/report.json) | Four five-million-step observations, ten windows each |
| [Fresh starts](long-01/fresh-starts.json) | 128 starts for each portrait parameter, all starts and window rates |
| [Precision](precision.json) | Independent explicit cubic at 70 and 100 digits |
| [Library](library.json) | Actual first deflation-unit loop in scikit-learn, rounded literal data |
| [Portrait 0.04](long-01/portrait-040.json), [portrait 0.15](long-01/portrait-150.json) | Separate unthinned 500,000-state inputs and their hashes |
| [Omitted arrays](omitted-data.json) | Six generated payloads, byte counts, shapes and recovery instructions |
| [Source review](source-review.json) | Read sources, search scope and unresolved priority questions |
| [Environment](environment.json), [manifest](manifest.json) | Executed source and compact evidence hashes |

The mathematical argument is written in the report. The symbolic check verifies
its finite algebraic identities and endpoint inequalities; it is not a formal
verification of topology in a proof assistant. Numerical roundoff checks are
separate from the exact real-arithmetic proposition.

## Reproduction

The executed sources are preserved byte for byte in [code](code). Use Python
3.12 with the versions in [environment.json](environment.json). The short
symbolic check requires NumPy, Numba and SymPy; numerical and library stages
also use mpmath and scikit-learn. The public figures use the repository's
Matplotlib/Pillow renderers.

From the repository root, choose a **new external directory**:

```powershell
python research/fastica_tripling.py --output D:/Codex/fastica-tripling-recheck --stage exact
python research/fastica_tripling.py --output D:/Codex/fastica-tripling-recheck --stage controls
```

The driver copies the frozen sources to that directory before running them.
`--stage all` runs `pilot`, `exact`, `controls`, `long`, `precision` and
`library`; each stage can also be selected separately. Pilot and long stages
refuse to reuse existing run directories. The long stage includes fresh-start
cohorts and separate consecutive-state portrait inputs. Set
`NUMBA_NUM_THREADS` to control CPU use. No GPU is required.

The baseline initial vector is `(0.83,0.29,-0.13,0.11)`, normalized before
iteration. Pilot clouds have stride 10 and are used only for screening;
published portraits use every consecutive state. For the presentation input,
project `long-01/portrait-040.npz` and `portrait-150.npz` using
`points = states[:, :3]`. This discards the fourth coordinate; it does not
reconstruct a three-dimensional state. The saved [camera metadata](presentation)
uses equal units and the same center for both panels.

All six large generated arrays are omitted from Git and recorded in the global
[data registry](../../records/generated-data.json). The archived bytes are
available in the local D: cache. They were never committed; an external clone
must recompute them. Chaotic trajectory bytes can differ with the numerical
environment. A hash mismatch against archived images is not a failure of the
exact theorem.

To check retained evidence without downloading trajectories:

```powershell
python research/verify_fastica_tripling.py
python tools/check_portraits.py
```

The first command verifies hashes, complete numerical windows and the stated
arithmetic margins. It does not rerun the long study. The normal portrait
checker verifies image receipts and counts without requiring omitted arrays.
