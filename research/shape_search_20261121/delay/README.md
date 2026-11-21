# Delayed population maps

The numerical search sampled 24,000 delayed crowding maps, 48,000 delayed recruitment maps and 32,000 recruitment maps around selected seeds. The study also includes three points of a parameter rectangle with a symbolic entropy certificate. The 22 long-run records include two stable controls. Each record contains 2,000,000 measured iterates in eight windows, a separate calculation in logarithmic coordinates, a QR determinant check and twelve additional starting states. The saved portraits use consecutive native abundances. [Parameters, samples and numerical records](../../../results/shape_search_20261121/delay/manifest.json)

## Dynamics

The state records three successive abundances, in order of decreasing time. The delayed crowding map sends $`(x,y,z)`$ to $`(xE,x,y)`$, with $`E=\exp(r(1-ax-by-cz))`$. The delayed recruitment map takes the form

```math
F(x,y,z)=(f(x,y,z),x,y),\qquad
f(x,y,z)=sx+(1-s)z\exp\!\left[r(1-ax-by-cz)\right],\quad c=1-a-b .
```

Case 1008712 gives a folded ring with finite-time exponents approximately 0.068, 0.041 and âˆ’0.263 per update. Case 1025516 gives three positive rates, approximately 0.105, 0.077 and 0.029. The third direction remained positive in each of five independent 224-bit finite-difference windows. Case 1031204 approaches three interleaved copies of a one-dimensional Ricker map. At the exact uncoupled limit $`a=b=s=0`$, the third iterate acts independently on the three coordinates. Several sampled positive-spectrum clouds are consequently close relatives. Delayed Ricker models and rigorous chaos results for related recurrences predate this experiment. [Independent precision checks](../../../results/shape_search_20261121/delay/mpfr-local-long.json) Â· [Earlier work on a three-dimensional delayed Ricker map](https://www.sciencedirect.com/science/article/pii/S0377042719300214)

## A positive-entropy local attracting set

For every point in the four-dimensional closed rectangle below, the map has a connected local attracting set with topological entropy at least $`\log((1+\sqrt5)/2)`$ per update.

| Parameter | Closed interval |
| --- | --- |
| Growth parameter $`r`$ | [4.48, 4.52] |
| Weight $`a`$ | [0.034, 0.035] |
| Weight $`b`$ | [0.024, 0.025] |
| Survival fraction $`s`$ | [0.00045, 0.00055] |

Let $`B=[10^{-22},10]^3`$. Throughout $`B`$, the new first coordinate lies strictly between its two faces. In bounding its minimum, the function $`z\exp(-rcz)`$ takes its minimum at an endpoint. Its maximum over positive $`z`$ equals $`1/(rce)`$. The shift coordinates then give $`F^3(B)\subset\mathrm{int}(B)`$ and $`F(B)\subset B`$. The nested sets $`F^n(B)`$ are compact and connected, and their intersection $`A`$ is invariant and attracts a neighborhood. [Rational endpoint and trapping inequalities](../../../results/shape_search_20261121/delay/entropy-certificate.json)

Put $`L=[0.0037,0.209]`$, $`H=[0.38,2.85]`$ and $`C=[0.0037,2.85]`$. For every $`x,y\in C`$, the function $`z\mapsto f(x,y,z)`$ increases strictly on $`L`$ and decreases strictly on $`H`$. Its endpoint bounds cover each of the transitions $`L\to H`$, $`H\to L`$, and $`H\to H`$. A transition refers to coordinates three times apart. To see this, write $`u_{j+3}=f(u_{j+2},u_{j+1},u_j)`$.

Consider a cyclic word of length $`3k`$ satisfying these transitions from position $`j`$ to $`j+3`$. The Cartesian product of its prescribed closed intervals is a compact convex box. Define a continuous self-map by replacing coordinate $`j`$ by the unique point on its prescribed branch which maps to $`u_{j+3}`$ with inputs $`u_{j+2},u_{j+1}`$. Strict monotonicity and uniform coverage give existence and continuity of this inverse. Brouwer's fixed-point theorem supplies a periodic orbit with the requested word. Its full orbit lies in $`A`$. Different words give separated starting states in the $`(3k,\epsilon)`$ dynamical metric for every $`0<\epsilon<0.171`$.

Each residue class modulo three can choose its cyclic word independently. The number of words equals $`\operatorname{tr}(M^k)^3`$, where $`M`$ is the two-state golden-mean transition matrix. Taking the exponential growth rate over $`3k`$ steps gives the stated bound. The closure of these periodic points gives a compact positive-entropy subsystem inside $`(L\cup H)^3`$. All nine strict inequalities use rational arithmetic with rigorously bounded exponential series. An independent 200-bit Arb implementation verified the same coverage over 256 separate boxes in the two conditioning coordinates, with 1,029 strict comparisons. [Arb recheck](../../../results/shape_search_20261121/delay/entropy-arb.json) [Verifier](certificate.py)

## Coexistence with a period-39 sink

At $`r=4.52`$, $`a=0.035`$, $`b=0.025`$, $`c=0.94`$ and $`s=0.00055`$, a separate contraction certificate isolates an attracting orbit of least period 39. The return map sends a closed parallelepiped strictly inside itself. The verifier bounds its derivative in the stored basis with 380-bit Arb balls and excludes periods 1, 3 and 13 by disjoint enclosures. One phase of the cycle has a coordinate above 6.73, separating the cycle from the symbolic subsystem. In each of three cohorts, all 4,096 randomly drawn starts approached the sink within the specified numerical tolerances, including starts in the symbolic intervals. The result illustrates the difference between sampling typical trajectories and detecting all invariant behavior. [Period certificate](../../../results/shape_search_20261121/delay/period-39.json) Â· [Basin experiment](../../../results/shape_search_20261121/delay/period-39-basins.json)

Along the diagonal joining opposite corners of the rectangle, the scan measures six QR windows after a 40,000-step transient at 1,201 points. The finite-time return test found 20 periodic samples with candidate periods 7, 21, 39, 78 and 156. The uniform entropy bound applies throughout that scan. [Parameter scan](../../../docs/assets/shape-search-20261121-delay-line.png)

## Long-period windows and transient growth

Ten parameter values in the entropy rectangle give independently isolated sinks with least periods 7, 21, 39, 78, 156, 312, 624, 1,248, 2,496 and 4,992. Each cycle has at least one phase above the symbolic box. The exact root for every proper divisor is excluded by disjoint enclosures. For periods through 1,248, the direct verifier encloses the full return in an affine basis. The longest two cycles use eight-step interval derivative blocks, approximate block centers checked by rigorous residuals, and balanced matrix products. Every intermediate block is checked against its proposed tube. The resulting bound on the Euclidean derivative norm is below 0.022 for the representative of period 2,496 and below 0.014 for period 4,992. [Periods 7–156](../../../results/shape_search_20261121/delay/period-windows.json) · [Periods 312–1,248](../../../results/shape_search_20261121/delay/high-period-windows.json) · [Long-cycle proofs](../../../results/shape_search_20261121/delay/long-cycle-blocks.json)

For the periodic branches of periods 39, 78 and 156, interval Newton calculations isolate transversal crossings of the multiplier −1 at diagonal coordinates approximately 0.9940966470, 0.9857790046 and 0.9842227814. Four unknowns describe the three cycle coordinates and diagonal parameter. Second derivatives of the map bound the derivative of the characteristic determinant, and strict Jury inequalities put the other two multipliers inside the unit circle. These certificates establish the multiplier crossings. The stability and minimal periods of the listed daughter cycles are independently certified at their stated parameter values. [Multiplier crossing certificates](../../../results/shape_search_20261121/delay/multiplier-crossings.json)

At the period-7 parameter, five of 256 uniformly sampled starts had still not matched the cycle after 1.5 million iterates. Three of those continued to give positive leading exponents in eight consecutive 100,000-step windows with independent tangent propagation. Extended direct iteration nevertheless captured all five, with the latest numerical match to the certified cycle at about 2.61 million updates. This experiment shows how long chaotic-looking transients can enter a finite-time search. [Recorded trajectories and capture times](../../../results/shape_search_20261121/delay/window-transients.json)

## Reproduction

The numerical computations use NumPy and Numba. The sink check uses python-flint 0.8.0 and mpmath. The independent finite-difference check uses gmpy2. RUN is an external scratch directory. Use Python with assertions enabled. The saved manifest is sufficient to regenerate any finalist without repeating the scout.

```text
python research/shape_search_20261121/delay/confirm.py results/shape_search_20261121/delay/manifest.json --cases 1008712 1025516 --long --out RUN/replay
python research/shape_search_20261121/delay/period_windows.py --verify --out results/shape_search_20261121/delay/high-period-windows.json
python research/shape_search_20261121/delay/block_cycles.py --verify --out results/shape_search_20261121/delay/long-cycle-blocks.json
python research/shape_search_20261121/delay/flips.py --verify --from-dir results/shape_search_20261121/delay --out results/shape_search_20261121/delay/multiplier-crossings.json
python research/shape_search_20261121/delay/certificate.py --out RUN/entropy.json
python research/shape_search_20261121/delay/sink.py --certificate results/shape_search_20261121/delay/period-39.json
python research/shape_search_20261121/delay/basins.py results/shape_search_20261121/delay/period-39.json --out RUN/basins
python research/shape_search_20261121/delay/continuation.py --out RUN/line --count 1201 --window 9000 --plot RUN/line.png
python research/shape_search_20261121/delay/precision.py results/shape_search_20261121/delay/manifest.json --cases 1025516 1031204 --out RUN/precision.json
python research/shape_search_20261121/delay/search.py --count 24000 --out RUN/scout
python research/shape_search_20261121/delay/search.py --count 48000 --recruitment --out RUN/recruits
python research/shape_search_20261121/delay/neighborhood.py RUN/recruits/screen.json --count 32000 --out RUN/local
```
