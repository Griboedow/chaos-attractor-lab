# Prior art for the Gaussian-posterior circle experiment

**Review date: 2026-11-01.** The averaging map, its Bayesian interpretation,
its proximity-operator representation, and the use of nonprojectors in
reflection algorithms are established ideas. The [current experiment](soft-constellations.md)
does not introduce these ingredients. Its contribution is a reproducible
temperature continuation and selected spatial numerical regimes. Historical
priority for those regimes remains unresolved.

## Equivalent descriptions of the averaging map

For sites `b_i` and positive temperature `tau`, define

```math
m_\tau(x)=\sum_i p_i(x)b_i,\qquad
p_i(x)=\frac{e^{-\|x-b_i\|^2/\tau}}{\sum_j e^{-\|x-b_j\|^2/\tau}}.
```

This is exactly the equal-weight, isotropic Gaussian **mean-shift update**,
with bandwidth `sigma^2=tau/2`. Equation 7 of
[Carreira-Perpinan, *Gaussian mean shift is an EM algorithm* (2007)](https://doi.org/10.1109/TPAMI.2007.1057)
gives this formula. The [author manuscript](https://faculty.ucmerced.edu/mcarreira-perpinan/papers/pami07.pdf)
was reviewed on PDF pages 1–4; page 2 was also checked visually. Its EM
interpretation and monotonicity concern iteration of `m_tau` itself, not
our composition with a circle reflector. Fractal convergence domains in
that paper must not be confused with a strange attracting orbit.

The same formula is the **posterior mean / MMSE denoiser** for a uniform
discrete prior on the sites observed through independent Gaussian noise
of covariance `(tau/2) I`. This is an exact identification, not a visual
analogy. Gaussian posterior means can be proximity operators of generally
nonconvex penalties; see Proposition 3 of
[Gribonval and Nikolova, *On Bayesian estimation and proximity operators*](https://arxiv.org/abs/1807.04021)
(journal publication 2021). The abstract, Section 1.3 and Proposition 3
were inspected; PDF page 8 was checked visually. The full proof and the
other noise models were not audited.

There is also a direct **softmax attention** expression. With columns
`B=[b_1 ... b_N]` and biases `c_i=-||b_i||^2/tau`,

```math
m_\tau(x)=B\mathrm{softmax}(2B^Tx/\tau+c).
```

[Ramsauer et al., *Hopfield Networks is All You Need* (ICLR 2021)](https://arxiv.org/abs/2008.02217)
uses the unbiased update `X softmax(k X^T x)` in equation 3 and relates it
to attention. Reviewed PDF pages 2–3 and 15–16, including a visual check
of page 3. If all sites have exactly equal norm, our biases cancel and
the updates coincide after identifying `k=2/tau`. For general sites the
biases remain; our stored floating-point polyhedra are not asserted to
have exact equal norms. The full reflected iteration is not the standard
Hopfield update. The paper's energy convergence statements do not prove
convergence of our composed map, nor do they assert convergence of every
state sequence without additional hypotheses.

## Explicit proximal representation

The following calculation makes the known connection concrete for this
experiment. Put

```math
\psi(x)=\frac{\tau}{2}\log\sum_i
\exp\!\left(\frac{2b_i\cdot x-\|b_i\|^2}{\tau}\right).
```

Then `grad psi=m_tau` and

```math
D m_\tau(x)=\frac{2}{\tau}\sum_i p_i(x)
(b_i-m_\tau(x))(b_i-m_\tau(x))^T\succeq0.
```

Let `psi*` be the convex conjugate and set
`phi(z)=psi*(z)-||z||^2/2`. Minimizing

```math
\tfrac12\|x-z\|^2+\phi(z)
=\tfrac12\|x\|^2-x\cdot z+\psi^*(z)
```

is equivalent, by Fenchel equality, to `z in partial psi(x)`. Since
`psi` is differentiable, the unique minimizer is `z=m_tau(x)`.
Thus the Gaussian mean is an exact proximal map for this penalty.
The penalty need not be convex. This is an application of established
convex-potential characterizations, not a new representation theorem;
see Corollary 1 of
[Gribonval and Nikolova, *A characterization of proximity operators* (2020)](https://arxiv.org/abs/1807.04014).
The statement on PDF page 3 was checked in text and visually; the
general Bregman proof was not audited.

For the standalone iteration `y=m_tau(x)`, convexity gives

```math
E(y)-E(x)\le-\tfrac12\|y-x\|^2,
\qquad E(x)=\tfrac12\|x\|^2-\psi(x).
```

This excludes nonconstant finite cycles of the standalone update.
All its iterates after the first lie in the compact convex hull of the
sites. The reflected map in this experiment does not inherit this descent
inequality merely by containing the same averaging operator.

## Reflection and denoising precedents

- [Geppert and Plonka, *Frame Soft Shrinkage Operators are Proximity Operators*](https://arxiv.org/abs/1910.01820)
  (2019 preprint; 2022 journal publication) explicitly discusses replacing
  a constraint projector by a non-idempotent shrinkage map in RAAR.
  Reviewed the [author PDF](https://na.math.uni-goettingen.de/pdf/prox.pdf),
  pages 1–3 and the references on page 16. Its cited applications are
  Loock and Plonka (2014) and Loock's thesis (2016); those application
  bodies were not reviewed here. The shrinkage map differs from our
  finite-prior Gaussian mean, but the substitution principle is prior art.
- [Ryu et al., *Plug-and-Play Methods Provably Converge with Properly Trained Denoisers* (2019)](https://proceedings.mlr.press/v97/ryu19a.html)
  defines plug-and-play DRS/ADMM and provides contraction conditions.
  Reviewed the abstract and PDF pages 3–4. Its displayed DRS places
  the denoiser on the other side of the reflection composition from ours.
  Its strong-convexity and denoiser assumptions do not follow for a
  circle indicator and our posterior map. No theorem contradiction or
  direct operator identity with that particular formulation is asserted.
- [Lee et al., *Denoising low-intensity diffraction signals using k-space deep learning* (2021)](https://doi.org/10.1103/PhysRevResearch.3.043066)
  describes RAAR with a constraint modified using denoised diffraction
  data. Reviewed the primary PDF's page 3 text and the authors' repository
  README; local PDF retrieval returned HTTP 403 and no figure audit is
  claimed. Denoising fixed measurement data differs from evaluating our
  state-dependent finite-prior mean at every iteration.

The earlier [finite-alphabet audit](raar-constellations.md#prior-art-and-equivalent-terms)
and [projection-dynamics audit](raar-priority.md) remain applicable.
They already rule out a first-discovery claim for chaos in nonconvex
reflection algorithms.

## Alternating lobes are established dynamics

Strict alternation between two separated chaotic bands is not a new
mechanism. Section 2 and Figure 2(a) of
[Mori, Hata, Horita and Kobayashi, *Statistical Mechanics of Dynamical Systems* (1989)](https://doi.org/10.1143/PTPS.99.1)
describe this behavior for `H(x,y)=(1-a*x^2+b*y,x)` at `a=1.10357`,
`b=0.3`: successive iterates alternate between bands and the second iterate
stays in one band. The primary publisher PDF text on printed pages 5–7,
including equation (2.1), the surrounding description and figure caption,
was reviewed. The screenshot request failed; no visual inspection of the
original figure is claimed. The rest of the 63-page article was not audited.

Case 29's sampled two-lobe alternation therefore supplies a comparison
target, not an originality claim. Even a rigorous trapping and alternation
proof would establish that property for this model without making the
property itself new. Equivalence to a two-band map through nonlinear
coordinates, or to a higher-dimensional extension, remains unresolved.

## Remaining novelty question

Search terms included Gaussian mean shift, mean update, MMSE denoisers,
nonconvex proximal maps, plug-and-play DRS/ADMM, soft-projection RAAR,
soft quantization, softmax/Hopfield memory, and chaotic reflected dynamics.
The reviewed sources do not establish that our particular parameters and
portraits were previously published. They also do not establish that they
were not. Unread full bodies, other parameter regimes, nonlinear changes
of coordinates, return maps and symbolic equivalence remain open.
The finite 108-regime geometry comparison is a separate diagnostic, not
a substitute for these equivalence questions.
