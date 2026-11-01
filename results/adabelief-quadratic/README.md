# Scalar-quadratic AdaBelief evidence

See the [report](../../docs/adabelief-quadratic.md) for equations, interpretation and limits. Numerical candidates are not certified chaotic attractors or priority claims. All spectra are rates per discrete iteration.

| Path | Contents |
|---|---|
| `pilot.npz` | All 4,096 Sobol lanes plus five controls; four growth windows, terminal full states and x extrema |
| `screen-replay/` | Exact-array replay of the exploratory screen with the published driver |
| `gallery.npz` | Twelve exploratory full-spectrum checks, 100,000 observed iterations each; 25,000-state tails |
| `long/` | Five 2,000,000-iteration observations, full spectra in 10,000-step windows and 200,000-state unthinned tails |
| `ensemble/` | 64 histories across two parameter pairs; 24 autonomous starts and eight corrected-start histories per pair, including sign-related starts |
| `neighborhood/` | Two 25-point coefficient neighborhoods; 100,000 observed iterations per point |
| `epsilon/` | Three parameter pairs times five original epsilon values; 200,000 observed iterations each |
| `followups/` | Saved-endpoint continuation of rejected neighborhood lane 22, period-three control, and unstable explicit two-cycle algebra |
| `analysis/` | Full-state lag screen and reference-support distance diagnostics |
| `controls/` | Complex-step derivative, short-path, inverse, determinant and independent-QR checks |
| `author-check/` | Actual pinned upstream PyTorch histories, including continuous bias correction |
| `precision70/`, `precision90/` | Original-unit finite-secant calculations at two precisions and perturbation sizes |
| `options/` | Upstream rectification/AMSGrad negative controls and exact symbolic determinant check |
| `source/` | Executed upstream code, repository BSD notice and package MIT notice; no upstream binaries |
| `exploration/` | Frozen exploratory scripts; their absolute paths document the original local run. Use the portable `research/` drivers for reproduction |
| `environment-lock.txt` | Observed optional full-check environment; not a cross-platform install prescription |
| `provenance.json`, `manifest.json` | Source identity, exact-array replay record, periodic follow-ups and SHA-256 file manifest |

NPZ run fields: `initial` and `end` are normalized `(x,m,s)` states; `parameters` is `(beta1,beta2,epsilon/alpha^2,epsilon/alpha)`; `lengths` is `(warmup,observed,window,tail_length)`; `windows` contains the three QR rates; `moments` stores window means of `(x,x^2,m^2,s)`; `minimum` and `maximum` cover every observed iterate; `tail` is consecutive, unthinned data. The `corrected_warmup` flag distinguishes corrected cold-start preparation from the autonomous main observation.

The upstream code differs bytewise only by an added terminal newline. The executed source hash is frozen in `author-check/summary.json`, and the original Git blob identity is in `provenance.json`. Original license notices are retained. All generated results are ours.
