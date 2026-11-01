# Irregular weight dynamics in finite-sample FastICA

[Research index](../README.md#research-index) · [Reproduction](../research/fastica/README.md) · [Evidence](../results/fastica/README.md)

## Result

The standard one-unit FastICA update has persistent irregular trajectories on
the unit sphere for four fixed, small data sets in four dimensions. Case 204
is the strongest numerical candidate: its leading tangent growth rate is
0.14919 per iteration over five million steps; all 64 fresh starts and all 64
small input perturbations tested below retain positive finite-time growth.
Independent library execution and 70- and 100-digit arithmetic recover the
irregular regime. These observations do not certify an asymptotic attractor
or establish historical priority.

![FastICA weight trajectories for Cases 204 and 351](assets/portraits/fastica.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

`Case` identifies a frozen data set and initialization, not a measured quantity.
Both panels project the literal normalized weight onto $(w_1,w_2,w_4)$, followed
by the same rigid orthographic camera. Each contains one million consecutive
iterates after 10,000 discarded steps. The spatial scale is isotropic within
each panel and differs between panels; density saturation is independent.
Case 204 alternates between two antipodal components. Neither component is
created by reflecting data. Case 351 occupies one component for this start;
other starts often converge to regular motion.

## Map and input conventions

For fixed centered observations $x_i\in\mathbb R^4$, empirical whitening gives
$N^{-1}\sum_i x_ix_i^T\simeq I$. The stored binary64 observations, rather than
their random seeds alone, define the selected numerical problems. The update is

```math
u_i=x_i^Tw,\qquad
h(w)=\frac1N\sum_i x_i g(u_i)-\left(\frac1N\sum_i g'(u_i)\right)w,
\qquad F(w)=\frac{h(w)}{\|h(w)\|}.
```

The state lies on $S^3\subset\mathbb R^4$ wherever $h\ne0$. There is no changing
data, external forcing, momentum, feedback term or sign canonicalization.
The nonlinearities are the usual $\tanh(\alpha u)$, $u\exp(-u^2/2)$ and $u^3$.
This is the first deflation unit, before projections against previously
extracted components. The underlying algorithm is established
[FastICA](https://www.princeton.edu/~sswang/awake_cerebellar_ms_refs/hyvarinen99_ieee_trans_neural_netw.pdf).
We study its continued iteration, with the stopping rule separately audited.

The 576-case screen uses six finite-data families, six sample counts and three
nonlinearities. All designs and outcomes are retained. Only five cases passed
the initial four-window gate. Cases 204 and 246 use finite samples of independent
normal entries; Cases 15 and 351 use Student-$t$ entries. These are tests of a
finite empirical map, not assertions that identifiable non-Gaussian independent
sources exist in the underlying population.

## Long runs and independent controls

The spectrum is measured in the three-dimensional tangent space, using the
analytic normalized-map Jacobian and QR accumulation. Rates are per iteration.
The long run discards 10,000 steps and measures ten windows of 500,000 steps.

| Case | Samples | Nonlinearity | Mean tangent spectrum | Fresh starts passing / 64 |
| --- | ---: | --- | --- | ---: |
| 15 | 12 | Cube | $0.35135,-0.09946,-0.16526$ | 10 |
| 40 | 16 | Tanh, $\alpha=1.7961360840150344$ | Settles to a stable fixed point | 0 |
| 204 | 32 | Gauss | $0.14919,-0.47070,-2.39643$ | 64 |
| 246 | 64 | Cube | $0.32468,-0.10922,-0.96509$ | 64 |
| 351 | 16 | Gauss | $0.24069,-0.49633,-0.65848$ | 10 |

Each fresh start is independently sampled on $S^3$, followed for 5,000 discarded
and 200,000 measured steps. Passing requires four positive growth windows,
no detected recurrence up to sign at lags 1–128, and a nonvanishing update.
The recurrence threshold is $10^{-8}$. Case 15 has 50 detected regular outcomes
and four additional starts that fail the gate; Case 351 has 54 regular outcomes.
These fractions describe this finite ensemble, not proven basin measures.
Case 15 has a positive sum of tangent rates; no contracting-volume or fractal
dimension claim is made for it.

For Case 204, independently perturbing every data entry and then recentering
and whitening gives 32/32 passing runs at perturbation amplitudes $10^{-6}$ and
$10^{-3}$, and 31/32 at $10^{-2}$. Every outcome, including the failure, is
recorded. The tested budgets match the fresh-start experiment. This supports
local numerical persistence but does not prove an open parameter region.

The 320-state independent NumPy control differs from the compiled step by at
most $3.92\times10^{-15}$. Central differences on the sphere agree with tangent
Jacobian products to relative error $5.40\times10^{-10}$ at step $2\times10^{-6}$.
Independent mpmath calculations retain positive leading growth in all four
2,000-step windows for Cases 15, 204, 246 and 351 at 70 digits. Cases 204 and
351 also pass at 100 digits. These runs lift the stored binary64 input values
exactly and perform all subsequent normalization, updates and tangent growth
in high precision; they do not rewhiten a different data set.

## Stopping is a separate dynamical question

The [scikit-learn implementation](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.FastICA.html)
was exercised through `FastICA.fit`, with `algorithm='deflation'`,
`whiten=False`, explicit initialization and a callable nonlinearity delegating
to the installed implementation. An observation hook reads the first unit's
weight and terminates after its declared budget or completion. The library
iteration and nonlinearities are unchanged. Version and source SHA-256 are
recorded; the hook does not validate subsequent deflation units.

At tolerance $10^{-15}$, Cases 15, 204 and 351 retain positive growth throughout
50,000 measured library steps. Case 246 exits the first unit after 37,963
updates; continuing its recorded stopping state for another million measured
steps retains positive growth in every window. Thus the software's convergence
flag is not evidence that this state is asymptotically attracting. The final
partial window in its library receipt is explicitly excluded from full-window
comparisons.

At the usual tolerance $10^{-4}$, Case 204 exits after **six updates**. Continuing
that actual stopping state yields four rates between 0.14897 and 0.15073 over
200,000 steps. The experiment raises `max_iter` for the other cases; it is not a
claim that default software runs indefinitely. Continued map dynamics and the
output of a stopped numerical procedure are different experimental objects.

Case 40 illustrates the opposite error: 200,000 measured steps initially look
chaotic, but the five-million-step run reaches a stable fixed point. Its
one-step residual is $1.84\times10^{-16}$ and its three numerical tangent
multiplier moduli are approximately 0.02834, 0.21516 and 0.65719. All 64 fresh
starts eventually pass the recurrence test within the declared budget. This
case is excluded from the portrait and retained as a long-transient control.

## Prior art and equivalence limits

Audit date: 2026-11-01. The evidence supports particular finite empirical
regimes of an existing algorithm. It does not support “the first chaotic
FastICA,” a new dynamical class, or a numerical probability of historical
novelty.

| Source inspected | Relevance and inspection scope |
| --- | --- |
| [Hyvärinen, 1999](https://www.princeton.edu/~sswang/awake_cerebellar_ms_refs/hyvarinen99_ieee_trans_neural_netw.pdf) | Original update and nonlinearities; algorithm and local convergence sections. |
| [Wei, arXiv:1408.6693](https://arxiv.org/abs/1408.6693) | One-unit fixed points, spurious solutions and false convergence are established subjects. Model assumptions, fixed-point analysis and sampling discussion inspected. No blanket finite-sample convergence conclusion is inferred. |
| [Regalia and Kofidis, 2003](https://pubmed.ncbi.nlm.nih.gov/18238073/) | Author abstract describes step-size conditions for monotonic convergence. Full text was not recovered; an unconditional claim for the present update would exceed the inspected source. |
| [Kolda and Mayo, 2011](https://arxiv.org/abs/1007.1267) | Section 3.2 and Example 3.5 reproduce a nonconvergent higher-order power iteration from Kofidis and Regalia. This is a direct mechanism-level predecessor for the cube cases. |

The tensor connection follows algebraically. Let $M=N^{-1}\sum_i x_i^{\otimes4}$
and $C=N^{-1}\sum_i x_ix_i^T$. For the cube nonlinearity,

```math
h(w)=Mw^3-3(w^TCw)w.
```

With ideal whitening $C=I$ and $\|w\|=1$, this is $(M-3E)w^3$, where the
symmetric fourth-order identity tensor satisfies $Ew^3=\|w\|^2w$. The update is
therefore normalized fourth-cumulant tensor power iteration. For the stored
binary64 data, retain the displayed formula with $C$; whitening residuals are
small but not exactly zero. Cube-case geometry alone cannot establish a new
mechanism. The Gauss cases are nonpolynomial, which also does not establish
inequivalence to previously studied dynamical systems.

Searches included `FastICA chaotic`, `FastICA chaos convergence`, `FastICA
bifurcation`, `FastICA period doubling`, `FastICA strange attractor` and `tensor
power method chaotic`. Retrieved material did not identify these data-specific
portraits. Search coverage is incomplete; Douglas's early convergence papers
and a broader nonconvergent tensor-iteration comparison remain open. No global
conjugacy or shape-uniqueness exclusion has been established.
