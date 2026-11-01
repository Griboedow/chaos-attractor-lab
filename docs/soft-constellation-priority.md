# Prior art for the Gaussian-posterior circle experiment

**Review date: 2026-11-01.** The averaging map, its Bayesian interpretation,
its proximity-operator representation, and the use of nonprojectors in
reflection algorithms are established ideas. The [current experiment](soft-constellations.md)
does not introduce these ingredients. Its contribution includes a reproducible
temperature continuation and, for Case 156, a [certified attracting
positive-entropy set with mixed-index periodic points](soft-symbolic-dynamics/manuscript.md).
Historical priority for this realization remains unresolved.

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

Two further convergence results delimit the comparison; our experiment
does not refute either theorem.

- [Teodoro, Bioucas-Dias and Figueiredo, *Scene-Adapted Plug-and-Play
  Algorithm with Guaranteed Convergence* (2018)](https://arxiv.org/pdf/1801.00605)
  fixes the posterior mixture weights during ADMM. Its resulting linear
  denoiser is a convex proximity operator. The text explicitly distinguishes
  this modification from a pure MMSE denoiser, which may be expansive.
  Corollary 2 additionally assumes convexity of the other penalties. Our
  weights depend on the current state, and the circle or sphere indicator
  is nonconvex. Reviewed the abstract, introduction, Theorem 2, and Section
  IV's Corollaries 1–2 and surrounding remarks; no full proof or experimental
  figure audit is claimed.
- [Shoushtari et al., *Prior Mismatch and Adaptation in PnP-ADMM with a
  Nonconvex Convergence Analysis* (2023 preprint)](https://arxiv.org/html/2310.00133v1)
  allows expansive MMSE denoisers and nonconvex data fidelity, subject to
  additional assumptions. Assumption 2 requires continuously differentiable
  data fidelity; a circle or sphere indicator does not satisfy it. The
  theorem is therefore not a convergence certificate for our substitution.
  Reviewed the abstract, Algorithm 1, Section 3.2's assumptions and Theorems
  1–2; the full proofs and experiments were not audited. No exact identity
  between our relaxed reflected map and that ADMM formulation is asserted.

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

## Periodic structure and limited equivalence obstructions

Case 156 has certified periodic points with one and two unstable directions
inside the same attracting invariant set. Mixed unstable dimensions are
known; see [Viana and Grebogi (2000)](https://fisica.ufpr.br/viana/artigos/2000/mky.pdf)
and [Saiki, Takahasi and Yorke (2021)](https://arxiv.org/html/1811.01009v4).
The latter's hetero-chaos requires density properties not established here.
The [article](soft-symbolic-dynamics/manuscript.md) gives the reviewed scopes
and distinguishes our located cycles from those stronger known results.

Their certified multipliers exclude C1 neighborhood conjugacy to a
dissipative planar Hénon base with a uniformly contracting extra coordinate,
and to a fixed-time map of an autonomous three-dimensional flow. These are
explicit, restricted obstructions. They leave general three-dimensional
Hénon-like maps outside the scope of the multiplier argument alone.

The subsequent connectedness and orientation certificate strengthens the
comparison: F fails local injectivity at some point of K. Consequently
there is no ambient neighborhood topological conjugacy to a locally
injective map. This excludes invertible three-dimensional Hénon models,
including [Gonchenko and coauthors (2005), model (1)](https://www.ma.imperial.ac.uk/~dturaev/mypapers/bc2005.pdf)
and [Hampton and Meiss (2022), equations (1)–(2)](https://arxiv.org/html/2206.07855v1),
whose stated constant Jacobians are nonzero. The [article, Section 6.1](soft-symbolic-dynamics/manuscript.md#61-an-obstruction-under-continuous-changes-of-ambient-coordinates),
gives the proof using invariance of domain and locally constant orientation.
This is an established topological argument, not a new general mechanism.
It does not prove a collision inside K or exclude a conjugacy defined only
on K. Known noninvertible systems and higher-dimensional realizations
remain open comparisons.

## Canonical spin and neural formulations

For ideal cube sites $`b_\sigma=c+aQ\sigma`$, with $`Q^TQ=I`$ and
$`\sigma\in\{-1,1\}^3`$, set $`z=Q^Tx`$, $`d=Q^Tc`$ and
$`h=2a(z-d)/\tau`$. Expanding the squared Gaussian distances removes a
common factor; the remaining weight is proportional to $`e^{h\cdot\sigma}`$.
The eight-state partition function factorizes, giving the exact identity

```math
Q^Tm(x)=d+a\tanh h.
```

The hyperbolic tangent acts componentwise. Rotating the entire state,
including the site's center, preserves the unit sphere and conjugates the
ideal reflected map to

```math
z^+=\beta z+\beta\frac{2\mu(z)-z}{\|2\mu(z)-z\|}
 +(1-2\beta)\mu(z),\qquad \mu(z)=d+a\tanh h.
```

For tetrahedral sign vertices, restrict to $`\sigma_1\sigma_2\sigma_3=1`$.
Writing $`t_i=\tanh h_i`$ gives

```math
Z=4\left(\prod_i\cosh h_i+\prod_i\sinh h_i\right),\qquad
\mu_i(z)=d_i+a\frac{t_i+t_jt_k}{1+t_1t_2t_3},
\quad \{i,j,k\}=\{1,2,3\}.
```

This is an exact three-spin posterior with one parity constraint. It is
also the even-parity three-bit code under the correspondence
$`\sigma_i=(-1)^{q_i}`$. The quotient is an algebraic identity; stable
computation should retain the four shifted exponential weights, since
saturated floating-point tangents can produce a spurious zero denominator.
These reductions do not make the full reflected iteration identical to
standalone neural retrieval or belief propagation.

The identities concern ideal regular sites. The original sphere search
stores separately rounded coordinates; fitting a rotation or rounding its
parameters gives nearby inputs, not an exact conjugacy certificate for
those binary64 arrays. Translation alone also moves the sphere and is not
an allowed simplification with the sphere fixed at the origin.

Additional primary comparisons delimit the novelty claim:

- [Chen and Aihara, *Chaos and Asymptotical Stability in Discrete-time
  Neural Networks* (1997)](https://arxiv.org/html/chao-dyn/9701020v1)
  studies leaky sigmoid recurrences and snap-back repellers. Reviewed the
  abstract and Section 3.1, equations (2)–(6). Our explicit normalized
  feedback is not present in the displayed recurrence. Sigmoid feedback
  and chaos are established ingredients; no nonlinear conjugacy exclusion
  or independent audit of their chaos proof is claimed.
- [Verzelli, Alippi and Livi, *Echo State Networks with Self-Normalizing
  Activations on the Hyper-Sphere* (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6761167/)
  uses the autonomous map $`z^+=Wz/\|Wz\|`$. Reviewed equations (1)–(6)
  and their surrounding assumptions. Its normalized-linear power-iteration
  identity does not hold for our nonlinear posterior plus residual terms.
  The paper establishes prior use of radial neural normalization; its
  later Lyapunov claims were not independently verified here.
- [Sibel, Reynal and Declercq, *Evidence of chaos in the Belief Propagation
  for LDPC codes*](https://perso.etis-lab.fr/declercq/PDF/ConferencePapers/Sibel_2012_CMSIM.pdf)
  reports periodic and chaotic numerical behavior in iterative decoding.
  Reviewed the abstract and Sections 1–3.3 of the author-hosted conference
  manuscript. Its state consists of Tanner-graph messages, unlike our
  three-coordinate exact posterior and sphere feedback. A parity-code
  interpretation therefore does not establish a first chaotic decoder.
  Neither an equation-by-equation equivalence nor a figure audit was made.

## Polyhedral appearance has direct precedents

The comparison must include symmetric discrete maps, beyond familiar
three-dimensional flows. [Reiter's 1996 paper, *Attractors with the Symmetry
of the n-Cube*](https://emis.de/ft/51972), constructs equivariant polynomial
maps and searches their coefficients. Its abstract and Sections 1–2
establish that cube symmetry is prior art. The author's abstracts for
[the three-dimensional cube study (1996)](https://webbox.lafayette.edu/~reiterc/abstracts/mv_1996c.html)
and [tetrahedral study (1997)](https://webbox.lafayette.edu/~reiterc/abstracts/mv_1997c.html)
likewise describe earlier spatial attractors. The latter two full papers
have not been reviewed here.

[Reiter's dodecahedral study (1999)](https://webbox.lafayette.edu/~reiterc/dodec/dodec02_pp.pdf)
uses group-averaged polynomial maps and explicitly describes diffuse clouds,
concentrated threads, and ball-like regions. Reviewed the abstract,
Sections 1–3 and the reproducibility appendices, and visually inspected
Figures 1–2 on the final page. These visual motifs and
polyhedral symmetry cannot support a general first-discovery claim for our
portraits. Formula differences and a finite cloud-distance comparison
still leave nonlinear equivalence and other published parameters open.

## Intrinsic collisions and homoclinic geometry

The [canonical cube certificate](canonical-sphere/manuscript.md#5-certified-attraction-merged-histories-and-a-homoclinic-return)
proves noninjectivity on the maximal attracting set of a specified trapping
box, and a homoclinic connection to a symmetric fixed point. The former
excludes topological conjugacy of that restricted map with any injective
map, including a fixed-time flow map. It does not identify that whole set
with the numerical portrait, or exclude known noninjective models.

[Jelleyman and Osinga, *Matching geometric and expansion characteristics
of wild chaotic attractors* (2022)](https://link.springer.com/article/10.1140/epjs/s11734-022-00440-w)
provides a direct conceptual precedent: Sections 1–2 describe a complex
angle-doubling map, reflection symmetry, multiple preimages and possible
self-intersections of unstable sets. Reviewed those sections and equation (1);
the map is planar and its formula differs from the three-coordinate posterior
reflection. Its general merging and homoclinic mechanisms are prior art.
Neither its wild-chaos criteria nor its attractor identification is imported
into our proof. The present audit did not independently reproduce its
numerics or visually review its figures; the publisher PDF download returned
HTML and the page-image retrieval did not yield an inspected figure.

The concrete addition is a reproducible certificate for the exact Case 653
inputs, with an explicit invariant-set scope. It is not a claim to have
invented symmetry-induced collisions, homoclinic orbits, or a universal
test of historical novelty.

## Product alphabets and feedback spin models

The [product-alphabet extension](product-alphabet/manuscript.md#2-model-factorization-and-established-interpretations)
uses independent two-, three- and four-level coordinates. Factorization
does not introduce a new posterior. Its ideal binary response is the
independent-spin tanh formula; its ideal ternary response is the standard
spin-1/Blume–Capel formula. Earlier chaotic spin-1 recursions and feedback
Potts flows therefore belong in the comparison, even when their papers do
not use denoising or projection terminology.

The extension's [prior-art section](product-alphabet/manuscript.md#6-prior-art-and-geometric-comparison)
and [review receipt](../results/product-alphabet/literature-review.json)
record the exact source scope. The supported contribution is a reproducible
four-wing parameter regime and controlled numerical geometry, with regular
nearby outcomes retained. A different displayed recurrence and a poor finite
cloud fit do not establish historical uniqueness or nonlinear inequivalence.

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
