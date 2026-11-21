# Deformed torus posterior-reflection maps

The GTX 1080 screened **2,359,296** new parameter configurations. Of 84 CPU finalists, **78** sustained three positive finite-time QR rates for five million iterations and passed the multistart and derivative checks. Sixty-seven passed all 17 starts. The earlier 31 three-positive records seeded the search. The original surface had identity scaling; both new streams change the physical surface. These counts describe parameter configurations. They do not establish 78 distinct attractor classes or priority in the literature.

| Stream | Distinct inputs | Short flags | Deep GPU evaluations | Deep flags | CPU finalists | Retained |
|---|---:|---:|---:|---:|---:|---:|
| Affine | 1,572,864 | 1,937 | 17,280 | 3,361 | 42 | 39 |
| Quadratic shear | 786,432 | 77 | 8,640 | 134 | 42 | 39 |

## Physical map

For sites $a_i$, use normalized Gaussian weights proportional to $\exp(-q\|x-a_i\|^2)$, where $q$ is the stored inverse temperature in `affine_model.py`. Their weighted mean is $m$. Let $P$ be the nearest projection onto the unit-major-radius circular ring torus with tube radius $r$. For $A=\operatorname{diag}(s_x,s_y,s_z)$ and $\Phi(u)=A(u_x,u_y,u_z+w(u_x^2-u_y^2))$, set $Q(v)=\Phi(P(\Phi^{-1}(v)))$. The iteration is $F(x)=\beta x+\beta Q(2m-x)+(1-2\beta)m$. The affine stream uses $w=0$. The sheared stream samples nonzero $w$ of either sign. Gaussian distances remain Euclidean in physical coordinates for both streams. The shear yields a retraction onto a deformed surface.

For $0<\beta<1$, putting $B=(\max s_i)((1+r)+|w|(1+r)^2)$ and $S=\max_i\|a_i\|$ gives $\|F(x)\|\le\beta\|x\|+\beta B+|1-2\beta|S$. Any radius greater than $(\beta B+|1-2\beta|S)/(1-\beta)$ therefore absorbs every initial state. This bound gives global boundedness. Existence of a specific hyperchaotic invariant probability measure remains open.

## Validation and comparisons

GPU screening used 1,024 burn-in steps and four 384-step windows. Each batch sent up to nine cases per anisotropy stratum to 24,576 burn-in steps and four 4,096-step windows. The table gives counts passing minimum-of-all-rates thresholds of 0.01 and 0.005, respectively. The deep queue retained quota winners even below the short threshold. A matched 42-feature distance to the previous clouds and within the new candidates informed selection. Physical parameter rows were hashed with BLAKE2b-128 and checked for duplicates per stream; zero/nonzero shear separates the two streams.

The independent binary64 CPU model propagated a full analytic Jacobian with orthogonalization. Each retained case passed ten consecutive 500,000-step windows after 100,000 burn-in steps (every rate >0.001), a centered derivative check with relative error below 1e-5, a lag-1-to-256 period screen and at least 9/17 multistart tests. Each start test used 32,768 burn-in steps and four 16,384-step windows. The case files include the deterministic start, 16 independent initial coordinates, all windows, shape variation, derivatives and an identity-chart counterfactual. The 250,000 consecutive XYZ samples are the tail of the five-million-step trace. Compact JSON covers all six rejected finalists as well.

| Case | Mean five-million-step spectrum per iteration | Positive starts | Additional evidence |
|---|---|---:|---|
| 409700, twin density plumes | 0.12052, 0.04923, 0.01663 | 17/17 | Minimum third rate 0.01642 in twenty five-million-step windows and 0.01432 in four 32,768-step MPFR windows |
| 2298665, lobes with intermittent phase slips | 0.12023, 0.03655, 0.00792 | 17/17 | Minimum third rate 0.00769 in twenty five-million-step windows and 0.00627 in four 32,768-step MPFR windows |
| 1403259, folded petals | 0.10496, 0.04601, 0.00877 | 17/17 | One short 8,192-step MPFR window had a negative third rate. All four 32,768-step MPFR windows were positive, minimum 0.00504; the 100-million-step window minimum was 0.00864 |

For 2298665 a plane estimated from the saved odd/even means separates the entire 250,000-point tail with a gap of 0.553 spatial RMS. A direct 50,000,000-step run found **791 same-halfspace consecutive pairs** among 49,999,999 comparisons. Sixteen other starts each alternated throughout 262,144 observed iterations. This gives a finite-time intermittent two-lobe pattern. [Halfspace control](parity-2298665.json)

Four matched 1,536-point windows from every retained case were compared to 108 archived regimes. Both headline cases ranked closest to JerkCircuit. Their distances were 1.580 (409700) and 1.499 (2298665), compared with within-pair sums of maximum deviations 1.184 and 1.076. The measurements depend on sampling, descriptors and feature scales. [Archive distances](archive-control.json)

[Cloud index and SHA-256](manifest.json) Â· [Search counts](summary.json) Â· [External checkpoint hashes](checkpoints.json) Â· [Long QR windows](extended) Â· [256-bit longer windows](precision-long)

## Reproduction

Run from the repository root using Python with NumPy, Numba, SciPy, PyOpenCL and gmpy2. GPU screening used NVIDIA OpenCL FP32. CPU decisions used binary64. The MPFR check propagates neighboring trajectories at 256 bits and reorthogonalizes every four steps.

```powershell
$src = 'research/shape_search_20261121/extension/projection'
$study = 'results/shape_search_20261121/projection'
$out = 'results/shape_search_20261121/extension/projection'
$cache = 'D:/Codex/chaos-attractor-cache/shape-search-20261121/extension/projection'
python "$src/survey.py" --study $study --cache "$cache/screen" --batches 384
python "$src/rescreen.py" --screen "$cache/screen" --study $study --out "$cache/deep" --batches 384
python "$src/survey.py" --study $study --cache "$cache/screen-warp" --batches 192 --warp --parent-shortlist "$cache/deep/shortlist.npz"
python "$src/rescreen.py" --screen "$cache/screen-warp" --study $study --out "$cache/deep-warp" --batches 192
python "$src/confirm.py" --shortlist "$cache/deep/shortlist.npz" --output "$out/cases" --clouds "$cache/clouds"
python "$src/confirm.py" --shortlist "$cache/deep-warp/shortlist.npz" --output "$out/cases" --clouds "$cache/clouds"
python "$src/replay.py" --record "$out/cases/case-2298665.json" --cloud "$cache/replay-2298665.npz" --reference "$cache/clouds/affine-2298665.npz"
python "$src/replay.py" --record "$out/cases/case-409700.json" --extended "$cache/100m-409700.json"
python "$src/precision.py" --records "$out/cases" --out "$cache/mpfr" --cases 409700 2298665 --burn 32768 --window 32768
python "$src/band_control.py" --record "$out/cases/case-2298665.json" --cache-root "$cache/../.." --out "$cache/parity-replay.json" --length 50000000 --other 262144
python "$src/archive_control.py" --records "$out/cases" --cache-root "$cache/../.." --shortlist "$cache/deep/shortlist.npz" --warped-shortlist "$cache/deep-warp/shortlist.npz" --out "$cache/archive-replay.json"
python "$src/inventory.py" --study $study --cache $cache --output $out
```

On the original workstation, independent replays of both 409700 and 2298665 from their case JSON produced identical SHA-256 cloud bytes and zero difference across all 30 stored binary64 window rates. For fresh confirmations, point the output argument at an empty directory: the confirmation stage resumes by skipping existing case records.
