# Scalloped field amplitudes and coexisting cavity regimes

Case 45941 produces a fan of folded arcs in the physical complex-field plane. One trajectory draws the entire fan. Color records the time since the most recent detuning pulse. Case 86780, at a different setting of the same smooth map, supports two numerically distinct chaotic regimes with different detuning ranges and fast return times. Both have rigorous finite-time enclosures. The search used numerical trajectories throughout.

![Scalloped field plane and detuning recurrence](scalloped-cavity-fan.png)

<!-- drawing-date -->
*Drawn: **2026-11-22**.*
<!-- /drawing-date -->

### The cavity map

Let the field immediately before the next round trip be $`E=x+iy`$, the dynamic detuning be $`T`$, and the normalized intensity response be $`f`$. All quantities are nondimensional. A constant real pump adds to a damped, phase-shifted field each round trip. The detuning relaxes toward an intensity-dependent value. In Cartesian coordinates,

```math
\begin{aligned}
I&=x^2+y^2,& f&=I/(1+I),&\phi&=\delta+kf+T,\\
x'&=P+a(x\cos\phi-y\sin\phi),&
y'&=a(x\sin\phi+y\cos\phi),&
T'&=(1-h)T+hgf .
\end{aligned}
```

Here $`P`$ is the pump, $`a`$ the field retention, $`\delta`$ the static detuning, $`k`$ the instantaneous shift, $`g`$ the memory gain and $`h`$ the relaxation fraction. The optical ring-cavity map and delayed-response instabilities have established antecedents. This calculation explores the particular simultaneous fast and relaxing saturable terms above. No experimental calibration or priority claim is available for these parameter sets. [Ikeda (1979)](https://doi.org/10.1016/0030-4018(79)90090-7), [Ikeda, Daido and Akimoto (1980)](https://doi.org/10.1103/PhysRevLett.45.709), [Ikeda and Akimoto (1982)](https://doi.org/10.1103/PhysRevLett.48.617)

For $`0<a,h<1`$, the triangle inequality gives $`|E_n|\le a^n|E_0|+P(1-a^n)/(1-a)`$. The detuning's distance to the interval with endpoints $`0,g`$ decays at least as $`(1-h)^n`$. The map is a diffeomorphism. From $`E',T'`$ recover $`I=|E'-P|^2/a^2`$, then $`T=(T'-hgf)/(1-h)`$ and $`E=(E'-P)e^{-i(\delta+kf+T)}/a`$. Before the final translation, the transformation in squared radius, angle and detuning has constant Jacobian determinant $`a^2(1-h)`$. Continuity extends the formula to the origin.

### Cases and evidence

| Case and orbit | $`P`$ | $`a`$ | $`\delta`$ | $`k`$ | $`g`$ | $`h`$ | Mean Lyapunov spectrum per iteration |
|---|---:|---:|---:|---:|---:|---:|---|
| 45941, scalloped fan | 2.20306635 | .96300721 | −1.33054030 | 9.05175877 | −1.67616093 | .69594193 | (.09617, −.15696, −1.20514) |
| 86780, lower detuning | .68380272 | .47694415 | .64927274 | −17.94019890 | 19.38252068 | .00387439 | (.03432, −.37055, −1.14837) |
| 86780, upper detuning | .68380272 | .47694415 | .64927274 | −17.94019890 | 19.38252068 | .00387439 | (.21472, −.42000, −1.27932) |
| 56729, alternating components | 2.07087350 | .92762107 | 2.81838250 | 3.30097151 | −17.03116226 | .01012690 | (.05810, .03249, −.25103) |

The identifiers index the original Sobol sequence from zero. The standard initial condition is $`(.173,-.341,.271)`$. The upper regime of Case 86780 starts at $`(.7124512033204118,-.4674104975241916,9.823930680698279)`$. Each table row uses 131,072 transient iterations and twelve 262,144-step full-tangent windows. The archived cloud stores the final 786,432 consecutive native Cartesian states without rescaling the axes. A second implementation uses complex exponentials and the alternative expression $`f=1-1/(1+I)`$. Its statistics may be compared with the two halves of the primary orbit. [Exact parameters, initial states, hashes and cloud filenames](manifest.json)

**Case 45941.** In 786,432 samples the detuning produced 38,493 peaks at the 0.10 prominence and ten-step minimum spacing used in the figure. Of 38,492 return intervals, 21,734 lasted exactly twenty steps (56.46%). The median was 20 and the 1st and 99th percentiles were 11 and 29. Moving the height cutoff from −1.55 to −1.45 left the counts unchanged. All sixteen additional starts from an independently sampled bounding box reached similar measures. The leading exponent spanned .09471–.09789 across the twelve long windows and their sums agreed with the analytic log determinant to $`4.5\times10^{-14}`$. Base-2 Jensen–Shannon distances between native pairwise histograms of independent arithmetic implementations were .0384, .0247 and .0202. Splitting the original cloud in half gave .0492, .0323 and .0263. Holding the phase's detuning at the observed mean still produced positive growth, so the numerical evidence does not attribute all chaos here to memory. [Window and start controls](case-45941.json)

**Case 86780.** The same parameters give a thin, mostly three-step spiral at low detuning and a thicker spiral at high detuning. In the field plane, the median one-step and three-step chord lengths are .546 and .0017 for the lower range, compared with .0030 and .0085 for the upper range. The lower trajectory produced 1,828 slow peaks with a median spacing of 424 iterations. Freezing the phase's detuning at its observed lower-regime mean gave a negative leading exponent in the control. [Lower regime](case-86780.json), [upper regime and starts](case-86780-basins.json)

![Two field-plane trajectories in common coordinates and separate detuning traces](coexisting-cavity-regimes.png)

<!-- drawing-date -->
*Drawn: **2026-11-22**.*
<!-- /drawing-date -->

Across 2,313 initial conditions (nine common field initializations crossed with 257 equally spaced detunings), all post-transient 8,192-step windows stayed on one side of $`T=7.5`$. The counts were 822 below and 1,491 above. These equal-spaced slices measure sample counts and do not estimate basin volumes. With eight fresh runs of the complex recurrence, every one of the 134,217,728 measured steps remained in its respective range. The largest lower-regime detuning was 6.10748 and the smallest upper-regime value was 7.95709. The two chosen exact trajectories also have **rigorous finite-time separation**. Arb ball arithmetic encloses every iterate from 2,048 through 8,192 inclusive inside $`5<T<6.3`$ and $`7.6<T<11`$, respectively. Every coordinate-ball radius stayed below $`2^{-64}`$. Inputs are the exact binary rational values of the archived floats. These tests leave permanent separation of the chaotic invariant sets open. [Numerical scan](case-86780-basins.json), [24,576-bit lower enclosure](arb-lower.json), [40,960-bit upper enclosure](arb-upper.json)

There is a unique fixed point, with $`6.88810<T<6.88811`$ in the observed gap. Its intensity is the positive root of

```math
I\left[1+a^2-2a\cos\left(\delta+(k+g)\frac{I}{1+I}\right)\right]=P^2 .
```

For these parameters, $`k+g>0`$ and $`0<\delta<\delta+k+g<\pi`$. The left side is strictly increasing from zero to infinity. Arb bisection and characteristic-polynomial sign changes bracket its three multipliers in $`(.05780,.05783)`$, $`(1.00159,1.00162)`$ and $`(3.91303,3.91306)`$. Thus the equilibrium has one stable direction and two unstable directions. Its stable manifold has codimension two. The calculation does not identify a basin-dividing surface. [Equilibrium certificate](case-86780-equilibrium.json)

**Case 56729.** Two disjoint field lobes alternate on every stored iterate. The observed real-field ranges of the two parity classes were $`[-.456,.321]`$ and $`[1.890,2.805]`$ (outward-rounded). The second exponent remained positive in every long window, from .03141 to .03330. Four of the sixteen additional starts followed the two-positive-exponent regime and twelve had negative leading exponents. The latter runs show sensitivity to initial conditions when selecting the observed regime. [Spectrum and starts](case-56729.json), [parity extrema](manifest.json)

### Selection and replication

A scrambled six-dimensional Sobol design generated 131,072 distinct parameter vectors. A GTX 1080 OpenCL screen used 3,584 burn-in iterations and four 512-step windows. Of 54,340 with leading growth above .012 in all windows, 51,563 also passed a minimum covariance-eigenvalue fraction of .0003. Deterministic farthest-point selection in growth, covariance, detuning and fluctuation statistics reduced these to 2,400. A second GPU run used 24,576 transient iterations and four 4,096-step windows. The minimum leading rate exceeded .01 for 2,308 cases. A separate double-precision full-spectrum screen found 165 of the 2,400 above .01 in the second exponent in each window. The orbit geometry, return lags and projection occupancy produced a 144-case visual comparison. Eighteen cases received new 131,072-point visual checks and the three parameter settings above received long controls. The search targeted no silhouette. Finite exponents are numerical evidence; a global trapping proof for either chaotic regime remains open. [Search and validation code](../../../research/wow_search_20261122/orbit), [selection and measurements](manifest.json)

From the repository root, set the cache path to an external writable directory. The cloud-only commands consume the committed JSON records directly. They require NumPy, Numba and SciPy. The scan also requires PyOpenCL and Pillow, the plots use Matplotlib, and the enclosure checks use python-flint 0.9.0.

```powershell
$code = 'research/wow_search_20261122/orbit'
$data = 'results/wow_search_20261122/orbit'
$cache = 'D:/Codex/chaos-attractor-cache/wow-search-20261122/orbit'
foreach ($case in @(45941, 86780, 56729)) {
    python "$code/validate.py" --record "$data/case-$case.json" --cache $cache --cloud-only
}
python "$code/basins.py" --record "$data/case-86780.json" --cache $cache --out $data --cloud-only
python "$code/plot.py" --cloud "$cache/case-45941.npz" --out "$data/scalloped-cavity-fan.png"
python "$code/plot_basins.py" --cache $cache --out "$data/coexisting-cavity-regimes.png"
python "$code/finite_enclosure.py" --record "$data/case-86780.json" --branch lower --bits 24576 --out "$data/arb-lower.json"
python "$code/finite_enclosure.py" --record "$data/case-86780.json" --branch upper --bits 40960 --out "$data/arb-upper.json"
python "$code/fixed_point.py" --record "$data/case-86780.json" --out "$data/case-86780-equilibrium.json"
```

For the complete screen, run `survey.py`, `morphology.py` and `spectra.py` with the same `--cache` argument. Pass the eighteen IDs in `manifest.json` as positional arguments to `inspect_survivors.py`. Omitting `--cloud-only` from the validation commands repeats the numerical controls. Supply `--out` to direct the numerical records. Running `summarize.py` with `--cache` and `--out` regenerates the compact manifest. The external NPZs store a float32 array named `points` with chronological native coordinates (real field, imaginary field, detuning), plus the exact parameter vector and initial state as float64. Bitwise cloud replay was tested on the original Windows software stack. New platforms can compare statistics after chaotic trajectories diverge.
