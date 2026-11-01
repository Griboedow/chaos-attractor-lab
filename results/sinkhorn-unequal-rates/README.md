# Rectangular and unequal-rate matrix-scaling evidence

This packet retains 4,864 declared inputs, every flagged continuation,
escape and recurrence accounting, full tangent spectra and independent
probability-domain checks. Balanced-row Cases 150 and 987 reach step
106,040,000 and supply the two selected three-state portraits.

[Report and equations](../../docs/sinkhorn-unequal-rates.md) ·
[Main-page portraits](../../README.md#matrix-scaling-with-unequal-rates)

## Evidence inventory

| Location | Contents |
| --- | --- |
| [Dataset index](datasets.json) | All 26 pilot, continuation and control datasets, with original hashes and bounded-size JSON parts |
| [Common-rate pilot](pilot01), [first continuation](followup01), [long continuation](rectangular-long01) | 1,536 inputs; all 20 flags followed through step 6,060,000 |
| [Unequal-rate pilot](dual01), [first continuation](dual-follow01), [long continuation](dual-long01) | 2,304 inputs; all 12 flags and the subsequent escape |
| [Balanced-row pilot](balanced01), [first continuation](balanced-follow01), [long continuation](balanced-long01) | 1,024 inputs; all 57 flags, five escapes and one periodic recurrence |
| [Selected extension](selected-100m) | Twenty more five-million-step blocks for both selected cases |
| [Final spectra](spectra-selected-106m.json) | Two full QR spectra after step 106,040,000 |
| [Earlier spectra](spectra-balanced-1m.json), [rectangular spectra](spectra-rectangular-6m.json) | Five balanced-row and two common-rate comparisons |
| [Long precision windows](precision-balanced-6m) | 200,000 observed steps at each of 70 and 100 digits for both selected cases |
| [Short precision windows](precision-balanced01), [rectangular precision](precision-rectangular01) | Earlier 20,000-step independent controls |
| [State and parameter controls](controls-balanced01) | 58 runs per selected case, including ordinary Sinkhorn |
| [Accounting and recurrence](analysis.json) | Complete screen outcomes and final-window lag tests |
| [Spectral factorization](schur-check.json), [actual-state checks](actual-checks.json) | 240 moderate-state and 18 retained-state implementation checks |
| [Presentation](presentation) | Rigid coordinate rotations, sampling intervals and input hashes |
| [Sources](sources.json), [environment](environment.json), [manifest](manifest.json) | Literature scope, package versions and retained-byte inventory |
| [External arrays](arrays.json) | Trajectory filenames, dimensions, byte counts and hashes |

Each dataset index contains metadata and a list of parts. Concatenating the
`runs` lists reconstructs the original JSON values. Original source-document
hashes preserve parent-child lineage; splitting does not claim to preserve
the original whitespace bytes. Case identifiers are local to a screen and shape.

The initial unequal-rate driver tested volume contraction before swapping
rates on odd cases. This leaves 76 actual noncontracting inputs in that
screen. The actual determinant is retained for every input; all twelve
flags contract volume. The later balanced-row driver has no swap and tests
the rates actually used.

## Verify the retained evidence

```bash
python research/verify_sinkhorn_unequal_rates.py
python research/verify_sinkhorn_unequal_rates.py --numerical --workdir D:/Codex/sinkhorn-dual-recheck
```

The first command verifies hashes, all-input and all-flag accounting,
continuation lineage, spectra and control counts. The optional numerical
mode replays 240 factorization checks and 18 actual-state checks in a new
external directory. Neither command certifies infinite-time chaos.

## Reproduce calculations

Copy `code/` into a new external working directory. Use the recorded
[environment](environment.json). Set `NUMBA_CACHE_DIR`, `MPLCONFIGDIR`,
`TEMP` and `TMP` to external directories, preferably on D:, and use one
BLAS thread. Drivers use four Numba threads; the selected long extension
uses two, and arbitrary-precision jobs use two worker processes.

Run each sequence in order:

```bash
python survey_rectangular.py --output pilot01
python follow_rectangular.py --pilot pilot01 --output followup01
python extend_rectangular.py --input followup01 --output rectangular-long01
python survey_dual.py --output dual01
python follow_dual.py --input dual01 --output dual-follow01
python follow_dual.py --input dual-follow01 --output dual-long01 --continuation --block 250000
python survey_balanced.py --output balanced01
python follow_balanced.py --input balanced01 --output balanced-follow01
python follow_balanced.py --input balanced-follow01 --output balanced-long01 --continuation --block 250000
python extend_selected.py --input balanced-long01/3x2.json --output selected-100m
python controls_balanced.py --input balanced-follow01/3x2.json --pilot balanced01/3x2.json --output controls-balanced01
python spectrum_dual.py --input selected-100m/3x2.json --output spectra-selected-106m.json --cases 150 987
python precision_dual.py --input balanced-long01/3x2.json --output precision-balanced-6m --cases 150 987 --observed 200000
python check_schur.py
python check_actual_dual.py
python analyze_runs.py
python prepare_presentation.py --input selected-100m --output presentation01
```

The other retained spectra and short precision windows use the same drivers
with the parent file and case numbers named in their records. A rerun can
change a chaotic history or capture time; preserve that outcome rather than
substituting archived coordinates. Input equality and structural identities
are stronger reproducibility targets than pointwise chaotic agreement.

## Data and image storage

Large numerical arrays remain outside Git. Pilot files contain endpoints
and at most 2,048 final states. Continuations retain all 50,000 final states;
the selected extension retains all 250,000 final states. Unreached windows
contain NaNs. Earlier full histories were not stored. Precision exports
contain every observed state converted to binary64, with final states also
retained as high-precision decimal strings. The control file retains only
its explicitly named final 2,048 states.

The two portrait inputs contain the entire selected final windows, translated
and rotated orthogonally. Their hashes and recovery instructions are also in
the repository [data registry](../../records/generated-data.json). No numerical
download is required to read the report or view the images. Recreate the
inputs externally with `prepare_presentation.py`; the central portrait and
animation tools consume these registered arrays. No PDF or generated point
cloud is committed.
