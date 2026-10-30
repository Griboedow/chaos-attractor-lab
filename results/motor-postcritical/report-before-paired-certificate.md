# What follows the motor-cluster transverse instability

A stable periodic state with **four synchronized pairs** intervenes between the loss of stability of the four-plus-four cycle and the irregular eight-row motion. Its period approaches twice that of the parent cycle. The four pair trajectories are quarter-period shifts of one another. On further increasing source activity, this paired state itself loses transverse stability, with one negative multiplier crossing minus one for each pair.

These are numerical observations in the unchanged Alouges et al. motor equations. Periodic shooting, independent solver checks and complete Floquet calculations support the statements below; they are not interval proofs of the nonlinear bifurcations or of a chaotic attractor. The [earlier interval certificate](motor-cycle-certificate.md) concerns the parent four-plus-four cycle. Scientific priority remains unresolved. The broad mechanism of equivariant cluster splitting and period doubling already occurs in other oscillator models [2].

![Physical paired-cycle traces, transverse stability and near-onset splitting](assets/motor-postcritical-branches.png)

## Exact model and cluster restrictions

The [source equations](motor-rows-source-audit.md) have eight identical rows and the constraint $\sum_j X_j=0$. We vary only $\rho$, keep $b=1$, and use the source gain

```math
r=1+\frac{19/20}{\sqrt{32001/500000}},\qquad
\delta=\frac{(r-1)^2}{\rho^2}-1.
```

The simulations evaluate the gain in binary floating point. The displayed exact expression identifies the intended source path; the computations in this report do not bound its rounding error rigorously. Time is the original normalized time $\tau$. Increasing activity means decreasing $\rho$.

If a cluster contains a fraction $w_a$ of the rows, its exact restriction is

```math
\begin{aligned}
V_a&=-\rho X_a+r\left(s_a-\sum_b w_b s_b\right),\\
c_a'&=1-c_a-V_as_a,\qquad
s_a'=1-s_a+V_ac_a,\qquad X_a'=V_a,\\
\sum_a w_aX_a&=0.
\end{aligned}
```

Four clusters of size two have $w_a=1/4$ and eleven independent coordinates. One three-dimensional zero-sum row contrast inside each pair supplies the other twelve tangent directions. Thus the complete 23-dimensional spectrum consists of an eleven-dimensional quotient block and four three-dimensional transverse blocks. The implementation checks this restriction and its Jacobian against the original eight-row field. The reconstructed full monodromy agrees with an independent full-state integration to $8.2\times10^{-14}$ at the reference cycle.

## First split: a stable doubled periodic branch

The parent four-plus-four cycle loses transverse stability near $\rho_1=3.370218399579$. At

```math
\rho=3.3652183995790423,
```

shooting finds a nonconstant four-pair periodic solution with

```math
T=6.69044151030482.
```

Its parent four-plus-four cycle at the same parameter has period $3.34339822910285$. The ratio tends to two as the branch approaches the first crossing; away from the crossing the periods need not be in an exact factor-of-two relation.

At the reference four-pair cycle, the complete Floquet calculation gives one neutral multiplier and twenty-two multipliers strictly inside the unit circle within numerical resolution. The largest nonneutral multiplier is approximately $0.96958005$. The leading within-pair multiplier is approximately $0.87559992$, repeated four times. The computed cycle is therefore locally attracting in the full eight-row model, not just inside the paired restriction.

The upward $X_1=0$ section alternates two distinct points: the lag-one state difference is about $0.158$, while the lag-two difference tends to zero. Shooting closes the full period to $3.7\times10^{-15}$. Radau reproduces closure to $1.8\times10^{-14}$ and the eleven-dimensional fundamental matrix to $1.5\times10^{-13}$. A separate full 23-state variational integration closes to $1.1\times10^{-14}$ and checks the flow direction to $3.3\times10^{-14}$.

For cluster labels $0,1,2,3$, let $q_a=(c_a,s_a,X_a)$. The cycle obeys numerically

```math
(q_0,q_1,q_2,q_3)(t+T/4)
=(q_2,q_3,q_1,q_0)(t).
```

This is an order-four permutation. After half a period it swaps clusters 0 with 1 and 2 with 3. Direct integration over a quarter period gives a permutation-closure error below $2.5\times10^{-13}$. If $P$ is the corresponding eleven-dimensional permutation matrix and $A$ is the quarter-period fundamental matrix, the full quotient monodromy agrees with $(P^{-1}A)^4$ to $2.8\times10^{-12}$. These are checks of the spatiotemporal symmetry, not interval enclosures.

Continuation towards the first crossing gives the expected small-amplitude behavior. Define $A$ as the sampled root-mean-square distance in $(c,s,X)$ between the two pairs created from one former cluster, and $\epsilon=\rho_1-\rho$. Then:

| Epsilon | Pair-splitting RMS amplitude | A squared divided by epsilon |
| ---: | ---: | ---: |
| 0.0050 | 0.13423011 | 3.60354 |
| 0.0020 | 0.08493580 | 3.60705 |
| 0.0010 | 0.06006844 | 3.60822 |
| 0.0005 | 0.04247827 | 3.60881 |
| 0.0002 | 0.02686695 | 3.60917 |
| 0.0001 | 0.01899812 | 3.60929 |

The square-root scaling and stable multipliers support a supercritical period-doubling interpretation. They do not replace an equivariant normal-form calculation or a rigorous branch-existence theorem.

## Second split of the four-pair state

The paired branch can be continued past a second instability. A numerical root solve locates it at

```math
\rho_2\simeq3.27487471656754,\qquad
\delta_2\simeq0.314813893279832,\qquad
T_2\simeq6.81054204975887.
```

The eleven-dimensional paired quotient remains stable there. A negative multiplier in each within-pair block passes through $-1$, supplying four destabilizing directions in the full system. The quarter-period permutation makes the four transverse spectra identical. Centered differences give $dm/d\rho\simeq18.8271790$; halving the difference step from $10^{-4}$ to $5\times10^{-5}$ changes it by $6.6\times10^{-7}$. This is a numerical slope, not a validated derivative bound.

| Rho | Largest quotient nonflow modulus | Leading within-pair multiplier |
| ---: | ---: | ---: |
| 3.285218400 | 0.589636 | -0.802351 |
| 3.280218400 | 0.616256 | -0.898632 |
| 3.275218400 | 0.644626 | -0.993526 |
| 3.270218400 | 0.673307 | -1.087106 |
| 3.265218400 | 0.702343 | -1.179434 |

The finite survey also finds irregular positive-growth trajectories farther along the activity path. Establishing the intervening branches and their basins requires additional evidence; the second multiplier crossing alone does not prove a route to hyperchaos.

## Beyond the second split: periodic saddles and a connection candidate

A second targeted screen contains fifteen parameters near the second crossing and five perturbations of the paired cycle at each parameter. All 75 lanes complete 3,000 discarded and 2,000 observed time units, followed by another 10,000 observed units for the entire fixed cohort. No lane is selected out. The two GPU stages take 156.81 and 321.19 seconds, respectively, and report no invariant-domain or numerical-failure flags. These flags do not detect the synchronization problem described below.

An apparently regular period-four state at $\rho=3.268$ has two synchronized pairs and four individual rows. In the zero-based data convention its groups are $(0),(1),(2,3),(4),(5),(6,7)$. Its period is $13.6331457678307$, and the 17-dimensional cluster restriction is attracting. However, the **full 23-dimensional spectrum has one expanding multiplier**:

```math
m_u\simeq1.125871594581,\qquad
\lambda_u\simeq0.008696267791.
```

The expanding direction splits rows 7 and 8 in one-based notation. The other synchronized pair, rows 3 and 4, has a leading contraction exponent approximately $-0.01631448594$. Thus the apparent period-four state is a saddle in the full model. The same calculation near $\rho=3.274$ finds a weaker saddle with expansion $0.00112263355$. A clean section plot is insufficient to infer full-system attraction.

We then keep rows 3 and 4 synchronized in an **exact seven-cluster invariant restriction**, and perturb the expanding contrast of rows 7 and 8. Both signs are tested, with perturbation sizes $10^{-7}$ and $10^{-9}$ in the original reduced-state norm. All four trajectories are followed for 400 saddle periods, or 5,453.26 time units. They leave the initial saddle and approach a phase-shifted, row-permuted copy of it. At the destination, rows 3 and 4 are now the expanding pair; rows 5 and 6 are the contracting pair. Because the expanding destination pair is held equal by the chosen invariant restriction, this experiment approaches the destination instead of following a subsequent departure.

![Growth, saddle matching, pair distances and the destination periodic motion](assets/motor-postcritical-saddle-connection.png)

The local growth and endpoint checks are quantitative:

- Stroboscopic splitting grows at $0.008696250$ per time unit, agreeing with the Floquet exponent to about $1.9\times10^{-8}$.
- Reducing the perturbation by a factor of 100 delays a fitted departure threshold by 529.5564–529.5565 time units. The linear Floquet prediction is $\log(100)/\lambda_u=529.5571$.
- Direct ODE phase matching puts the final states within $4.2\times10^{-13}$ per physical coordinate of a permuted reference cycle. The larger floor in the plotted matching distance comes from cubic interpolation of the reference samples.
- Halving the maximum integration step from 0.05 to 0.025 gives the same pair transition and destination. The maximum coordinate difference over the complete traces is $1.38\times10^{-6}$; the final-state difference is $1.06\times10^{-10}$. The trajectory comparison is not an error enclosure.

For one outgoing sign, a representative permutation from reference rows to current rows is

```math
g=(5,4,6,7,2,3,0,1),\qquad q_j^{\rm destination}=q_{g_j}^{\rm reference},
```

using zero-based indices. Its fourth power swaps the two members of all four original pairs. This operation agrees with a half-period shift of the saddle to $6.5\times10^{-14}$ in a direct integration. Thus the selected permutation generates a candidate four-saddle circuit **modulo phase**, with the expanding pair following this sequence in one-based row numbers:

```text
7-8 -> 3-4 -> 5-6 -> 1-2 -> 7-8
```

This is a numerical heteroclinic-network candidate: symmetry provides related candidate edges, and finite trajectories approach the expected periodic saddles inside the relevant invariant subspace. We have not certified the infinite-time connections or attraction of the network. The local contraction-to-expansion ratio is about 1.876; that ratio alone is not a network-stability theorem. Other outgoing signs and permutation choices may give a larger network than the displayed circuit.

## Floating-point synchronization can hide the next departure

The long GPU continuation gives an instructive failure of a naive classification. At $\rho=3.268$, one lane eventually has **bitwise equal** states for rows 3–4 and 5–6. Its carried tangent has bitwise equal components in these pairs too. The reported late leading exponent is then near zero, and its section looks periodic. Yet full periodic shooting and the complete monodromy recover the same expanding multiplier $1.125871594581$ after accounting for the new row permutation.

Once subtraction of nearly equal row states rounds a contrast to exactly zero, permutation symmetry makes that numerical subspace invariant. A single carried tangent can also lose its expanding component through the same cancellation. This may produce apparent permanent capture by a saddle. The observed low exponent therefore does not establish a stable second daughter, quasiperiodicity, or absence of switching. The exact-zero escape experiments above deliberately use an invariant restriction and are labeled as such; they are not generic full-system trajectories.

A controlled diagnostic restart uses the **same physical initial state twice**, changing only the tangent used to estimate growth. After 500 discarded and 1,000 observed time units, the old carried tangent gives approximately $0.000248$, while a fresh sine tangent gives $0.007552$. The latter recovers the expanding direction; neither finite-window value is expected to equal the asymptotic Floquet exponent exactly. This restart does not perturb the physical trajectory or turn its periodic motion into chaos. [Data and code](../results/motor-postcritical/tangent-trap-control/summary.json).

This issue is localized to the nearly synchronized states studied here. It does not reclassify the previously checked strongly separated hyperchaotic candidate, whose full independent CPU spectrum is a separate result. Reliable unrestricted switching calculations require coordinates that retain tiny row contrasts directly, or appropriately controlled higher precision.

## Survey design and evidence

The first screen includes twelve declared parameter values and five starts per value: three seeded generic contrasts of the parent cycle, one contrast preserving four pairs, and one independently specified source-style initial condition. Every lane uses 2,000 discarded and 2,000 observed time units. All sixty lanes finish without an invariant-domain or numerical-failure flag. The GPU integration takes 126.95 seconds on the GTX 1080, with relative tolerance $10^{-10}$ and half the existing adaptive step cap.

All five starts give four synchronized pairs immediately below the first crossing and at the sampled points down to $\rho=3.2902184$. CPU section samples distinguish periodicity from a merely small leading exponent. The closest point below the first crossing still has a small finite-time drift after 4,000 time units; the shooting and continuation, rather than that transient alone, establish the numerical periodic solution there. At $\rho=3.2502184$ and lower, some survey states have positive finite-time growth and incomplete clustering. A trajectory confined by exact or nearly exact initial symmetry can remain on an unstable lower-dimensional set, so its positive ambient exponent is not by itself evidence that its own motion is chaotic.

The GPU implementation and its independent field/Jacobian controls are inherited unchanged from the [row survey](motor-rows-gpu-survey.md). CPU shooting uses DOP853 at relative tolerance $2\times10^{-11}$ and maximum step 0.025, with tighter/implicit controls at selected cycles. The preserved data include all declared starts, final states, tangent vectors, growth windows, pair-distance statistics, shooting attempts and solver comparisons.

## Reproduce and inspect

New output directories are required because the scripts refuse to overwrite existing evidence. Run from the repository root with the simulation dependencies installed:

```bash
python research/analyze_motor_postcritical.py prepare --output outputs/postcritical-reproduction
python research/analyze_motor_postcritical.py survey --output outputs/postcritical-reproduction
python research/analyze_motor_postcritical.py sample --output outputs/postcritical-reproduction --lanes 18 --duration 150
python research/motor_postcritical_cycles.py --sample outputs/postcritical-reproduction/samples/lane-018.npz --groups '[[0,1],[2,3],[4,5],[6,7]]' --lag 2 --output outputs/paired-cycle-reproduction
python research/motor_postcritical_branch.py continue --source outputs/paired-cycle-reproduction --output outputs/paired-branch-reproduction
python research/motor_postcritical_escape.py run --output outputs/saddle-escape-reproduction
python research/motor_postcritical_escape.py analyze --output outputs/saddle-escape-reproduction
```

- [Complete first survey and design](../results/motor-postcritical/initial-survey/).
- [Reference doubled cycle, complete spectrum and independent controls](../results/motor-postcritical/paired-cycle-rho-3p3652184-retry/analysis.json).
- [Quarter-period symmetry and derivative reconstruction](../results/motor-postcritical/paired-cycle-rho-3p3652184-retry/symmetry.json).
- [Full monodromy block reconstruction](../results/motor-postcritical/paired-cycle-rho-3p3652184-retry/block-control.json).
- [Near-onset and subsequent paired-branch continuation](../results/motor-postcritical/paired-branch/analysis.json).
- [Numerical second-crossing refinement](../results/motor-postcritical/paired-second-crossing/analysis.json).
- [Complete second survey and fixed-cohort continuation](../results/motor-postcritical/second-survey/).
- [Six-cluster saddle with one expanding direction](../results/motor-postcritical/six-cluster-cycle-rho-3p268/analysis.json).
- [Same saddle after the observed pair transition](../results/motor-postcritical/final-six-cluster-rho-3p268-correct-groups/analysis.json).
- [Four outgoing-perturbation experiments](../results/motor-postcritical/saddle-escape/summary.json), [phase/permutation matches](../results/motor-postcritical/saddle-escape/diagnostics.json), [growth and step controls](../results/motor-postcritical/saddle-escape/growth-controls.json), and [permutation circuit](../results/motor-postcritical/saddle-escape/network-algebra.json).
- [Half-step trajectory repeat](../results/motor-postcritical/saddle-escape-half-step/summary.json).
- [Plot renderer](../research/report_motor_postcritical.py) and [final figure manifest](../results/motor-postcritical/figure-manifest-v2.json).
- [Identical-state carried/fresh tangent comparison](../results/motor-postcritical/tangent-trap-control/summary.json).
- [Exact-byte manifest and saved-data audit](../results/motor-postcritical/manifest-verified.json).

One initial command failed before integration because its provenance routine mixed relative and absolute paths. The failure record remains in the results; path resolution was corrected before the successful run. A subsequent shooting attempt used an outdated cluster grouping for a later trajectory and converged to a different, more symmetric saddle with two unstable directions. That projected result is retained with an interpretation note and excluded from the connection evidence. Executed script snapshots accompany the numerical outputs. No failed numerical trajectory was removed from the surveys.

## Relation to earlier work

The equations and source activity prescription come from Alouges et al. [1]. Haugland and Krischer [2] already demonstrated equivariant splitting of periodic clusters, a quarter-period rotating four-oscillator branch, and connections to chaotic dynamics in a globally coupled Stuart–Landau system. That system has different equations and state variables. Our numerical result identifies a concrete sequence and parameter values in the motor-row equations; it does not establish a new general route to chaos or a new attractor type.

1. F. Alouges et al., “Some mathematical models for flagellar activation mechanisms,” *Mathematical Models and Methods in Applied Sciences* **35**, 2395–2424 (2025). [DOI](https://doi.org/10.1142/S0218202525500423); [primary manuscript](https://arxiv.org/html/2409.03506v2).
2. S. W. Haugland and K. Krischer, “Connecting minimal chimeras and fully asymmetric chaotic attractors through equivariant pitchfork bifurcations” (2021). [Primary manuscript](https://arxiv.org/html/2102.10138v1), especially Figure 1 and the discussion of equivariant period doubling and the rotating four-oscillator state.
