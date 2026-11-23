# A driven cubic changes from unlinked cycles to a Brunnian link

A periodically forced gradient equation in one complex variable produces three attracting cycles whose collective link is selected by the forcing coefficients. For the primitive-period Case 48000001, the three cycles at gain 0.34 form the Borromean rings. Every pair can be separated by a sphere when the third cycle is omitted. The full triple cannot be split by a sphere. The [braid forces chaos in the gradient equation](../../../cubic-forced-orbits/). An associated expanding return map on the same cubic sheets has three chaotic solenoids with the same collective nonsplitting property.

The original search used 1,261,568 randomized Fourier parameter draws and a separate grid of 320,544 one- and two-term real drives. The two screens respectively returned 458 and 2,552 nonzero triple-invariant candidates with a margin from the cut threshold. The grid can represent phase shifts or repetitions of the same physical pattern. Eighteen representative drives were checked by interval arithmetic and by extracting geometric crossings from their actual cubic roots, including five deliberate frequency continuations. No strand trajectories or tube centerlines enter either search. [Coefficients and interval receipts](certificates.json) · [Independent geometric checks and point-cloud hashes](clouds.json)

## A single gradient field

Write θ for phase on the circle, z for a complex coordinate, κ for the positive response gain and b for a complex Fourier drive. The following autonomous flow on the phase circle times the complex plane uses one formula everywhere.

```math
p_\theta(z)=z^3-3z+b(\theta),\qquad
\dot\theta=1,\qquad
\dot z=-\kappa\,\overline{(3z^2-3)}\,p_\theta(z).
```

Case 48000001 uses the primitive-period drive below. The independent sparse-grid Case 47087142 has real part −4 sin(3θ) + 4 sin(4θ) and imaginary part 3 sin(7θ). It also has primitive period 2π. The simpler parent Case 47053641 omits the first harmonic of Case 48000001 and repeats its forcing three times per displayed revolution. The first harmonic in the primitive case removes that time symmetry.

```math
b(\theta)=\tfrac14\sin\theta-4\sin(3\theta)+3i\sin(6\theta).
```

All root and orbit diagrams use the cylindrical embedding ((R + Re z) cos θ, (R + Re z) sin θ, Im z), with R = 3.25 for these two primitive examples. The trajectory has three real state coordinates. The critical points of the cubic are the fixed sites z = ±1, where the critical values are b(θ) ∓ 2.

## Recovering the link from the coefficients

In the b-plane use the two real slits extending outward from ±2. A loop crossing the right slit winds around the double-root value +2 and exchanges the two rightmost roots of the base cubic at b = 0. The left slit similarly exchanges the two leftmost roots. The remaining slit domain is simply connected. With the two ray orientations reversed relative to each other, the crossing word therefore determines the braid up to choice of base path and orientation. The imaginary drive is B sin(kθ), so all possible slit crossings occur at the known phases jπ/k and are transverse. Arb decides which side of −2 and +2 contains the real drive at each such phase. The certificate then encloses b over the entire circle by adaptive subdivision and verifies a positive distance δ from both branch values.

For each retained case the strand permutation is the identity. All pairwise linking numbers vanish. Deleting one strand of a pure three-strand braid therefore leaves a pure two-strand braid of total exponent zero, which is the identity. In all eighteen cases the third longitude has nonzero Milnor invariant, computed in the integral Heisenberg group. In a separate check the actual roots of the cubic were tracked once around the circle. Two generic projection directions reproduced the pairwise values and the absolute triple invariant in every case. The waveforms include core absolute invariants 1, 2, 3, 5 and 7. The arguments use established braid and link invariants. [Source search and geometric check](../../../../research/shape_search_20261121/frontiers/cubic) · [Polynomial root braids](https://arxiv.org/abs/2006.00396) · [DeTurck et al. on the triple invariant](https://arxiv.org/abs/1207.1793)

## Rigorous gain regimes for the gradient flow

For a sufficiently small residual radius d below δ, the region |p| ≤ d has three disk components above each phase. Since it contains no critical point and has trivial sheet permutation, these form three solid tori. The disk fibers are inverse images under individual univalent branches of the cubic. On their boundary the upper derivative of |p| is at most max |b′| − κd min |p′|². A disk bound |z| < R₀ follows from R₀³ − 3R₀ > max |b| + d. The factorizations of p(z) − p(±1), together with max(|z−1|,|z+1|) ≥ 1, imply

```math
|p'(z)|^2\ \ge\frac{9(\delta-d)}{R_0+2}.
```

The real Jacobian in each phase has Euclidean logarithmic norm at most κ(−|p′|² + 6|z||p|). For the recorded d and gain this is strictly negative throughout the tubes and the boundary points strictly inward. Each one-turn map contracts its closed disk fiber in the intrinsic path metric: the flow maps connecting paths inside the moving tube and contracts their lengths. Completeness gives a unique attracting periodic orbit in each tube. Contracting the residual radially to zero in each inverse branch gives an isotopy from the three cycles to the three root graphs. This proves their collective link at each recorded large gain. For Case 48000001 the analytic threshold used in the archive is 16,384. It is a conservative sufficient value.

The small-gain limit gives a different link. On compact sets the one-turn stroboscopic map Pκ has (Pκ − id)/κ converging in the first derivative and in value to the averaged field multiplied by 2π. The mean of b vanishes in Case 48000001. The averaged field has three simple sinks at the roots of z³ − 3z and two simple saddles at ±1. The implicit function theorem gives one fixed point near each for all sufficiently small positive gains. No periodic orbit escapes the common disk |z| < 3 in this case. On its exterior, writing r = |z| and B = max |b| ≤ 7.25, the inward radial sign follows from r⁶ − 4r⁴ + 3r² − B(r³+r) > 0. Compactness and uniform convergence then exclude additional fixed points of the one-turn map in the small-gain limit. The three sinks are isotopic to the unlinked constant-root graphs.

It follows that some finite positive gain has a fixed point of the two-real-dimensional stroboscopic map with multiplier +1. Otherwise over a compact gain interval connecting the two regimes its fixed points would form a proper finite covering by the implicit function theorem and the uniform radial bound. The two fixed-point branches at ±1 are permanent. The other three would continue without collisions, yielding an isotopy between the unlink and the nonsplit triple, a contradiction. This implication does not determine the gain, bifurcation type or number of bifurcations.

## Finite-gain orbit checks

At gain 0.34 the three attracting cycles force [positive entropy and at least 21 one-period solutions in this flow](../../../cubic-forced-orbits/). An [advected ellipse](../../../cubic-material-lines/) shows the resulting folds.

Case 48000001 was scanned at 57 gains with 15 initial conditions. The observed responses at gains 0.03, 0.15, 0.23 and 0.34 are respectively three one-turn cycles, one observed one-turn cycle, one observed three-turn cycle and three one-turn cycles. Independent adaptive integration gives return residuals and Floquet multipliers for these examples. The fifteen starts do not exhaust the possible basins. [Gain survey and independent controls](transition.json) · [Other primitive and repeating drives](gradient.json)

For finite-gain proving, a piecewise-linear periodic center curve supplies a moving disk at each phase. On each phase cell Arb bounds the defect of the interpolant and the logarithmic norm of the true field throughout a disk of fixed maximum radius. If the defect is at most E and the logarithmic norm at most M, the upper radius over an interval Δ is bounded by exp(MΔ) times the input radius plus max(1, exp(MΔ)) EΔ. The exact scalar comparison solution is monotone, so checking both endpoint bounds below the disk cap validates the enclosure throughout the cell by first exit. A strict return of the starting disk and a negative integral upper bound for M give a unique stroboscopic sink within that disk. True orbits are homotopic to their center polygons inside disjoint moving disks. The center coordinates are interpreted as exact binary numbers for crossing order and the signs of projected crossings. The numerical centers are proposals; all inequalities are verified in Arb. [Finite-gain certificates and exact crossing words](finite-gains.json) · [Verifier](../../../../research/shape_search_20261121/frontiers/cubic/finite_gain.py)

At gains 0.03 and 0.34, all three initial disks return strictly inside themselves. Using 65,536 cells per turn, the largest final radius is below 0.00860 from 0.0088 at the lower gain and below 0.000796 from 0.0055 at the higher gain. The integral log-norm bounds are below −0.248 and −86 respectively. The exact polygon word is empty at gain 0.03 and is (σ₁⁻¹ σ₂)³ at 0.34. These certify an unlink and the Borromean rings for the actual attracting cycles at the two specified gains.

## Chaotic return on the same emergent sheets

Use a larger residual radius D = δ/2. The base phase expands and the residual contracts on each sheet, with one common rule for every initial state:

```math
\theta'=2\theta+0.12\sin(2\theta)+0.137\pmod{2\pi},\qquad
p_{\theta'}(z')=0.2p_\theta(z)+0.6D e^{i\theta}.
```

Choose the lifted solution on the same connected component of |p| ≤ D. The three components are found from the coefficients. This lift is single-valued since each component covers the phase circle times the residual disk with degree one. Input residual disks map strictly inside the target since 0.6 + 0.2 < 1. The two target images from different preimages of the same output angle are disjoint, with center distance 1.2D and combined radii 0.4D. The angular derivative exceeds one. The three attracting sets are degree-two solenoids, each with entropy log 2. This is a discrete map; its expanding base is specified explicitly. Chaotic attraction in the unswitched gradient flow is an open question.

The three compact solenoids inherit the Brunnian nonsplitting property from their traps: each depth multiplies the triple core invariant in absolute value by 2³. An enclosing separating sphere for the three limiting sets would separate the cores of some finite depth. The odd part of the original absolute triple invariant is an ambient-homeomorphism obstruction because the nesting degrees here all have only the prime factor 2. The certified sample therefore contains at least four distinct ambient types, distinguished by odd parts 1, 3, 5 and 7. For the proof of the winding law and the nested-torus comparison, see the [preceding three-component construction](../../extension/brunnian/README.md#passing-from-tubes-to-the-compact-sets). The historical priority of this particular oscillator and parameter transition is open.

## Reproduction

The computational packet uses Python with NumPy, SciPy, Numba and python-flint. Set OPENBLAS_NUM_THREADS to 1 for the screening run. Use an external RUN directory for numerical payloads.

```text
python research/shape_search_20261121/frontiers/cubic/search.py --per-pair 16384 --out RUN/search.json
python research/shape_search_20261121/frontiers/cubic/sparse.py --out RUN/sparse.json
python research/shape_search_20261121/frontiers/cubic/certify.py RUN/search.json RUN/sparse.json --out RUN/certificates.json
python research/shape_search_20261121/frontiers/cubic/sample.py RUN/certificates.json --out RUN/frontiers/cubic --manifest RUN/clouds.json
python research/shape_search_20261121/frontiers/cubic/gradient.py RUN/certificates.json --out RUN/gradient.json --adaptive
python research/shape_search_20261121/frontiers/cubic/transition.py RUN/certificates.json --out RUN/frontiers/cubic --manifest RUN/transition.json
python research/shape_search_20261121/frontiers/cubic/finite_gain.py RUN/certificates.json RUN/gradient.json --out RUN/finite-gains.json --cache RUN/frontiers/cubic
python research/shape_search_20261121/frontiers/cubic/figure.py RUN/frontiers/cubic/transition-curves.npz RUN/transition.json --out RUN/cubic-transition.png
```
