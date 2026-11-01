# RK45 geometry: priority and equivalence audit

Audit date: 2026-11-01. **Historical uniqueness is unresolved.** These are
implementation-specific numerical candidates in the existing Dormand–Prince
solver. A visually different cloud is not sufficient evidence of a new
dynamical mechanism or a new equivalence class.

The [RK23 audit](rk23-priority.md) supplies the shared Hall–Higham–Humphries
reference chain and access gaps. Its conclusions apply here too: adaptive
solution/step maps, unstable numerical equilibria, norm dependence and
tolerance-scale oscillations near a stable physical equilibrium are prior work.

## Closer sources checked for this extension

| Primary source and inspected scope | Relevant prior result | What remains open |
| --- | --- | --- |
| [Higham and Hall, 1990, *Embedded Runge–Kutta formulae with stable equilibrium states*](https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P6.pdf). All nine pages visually inspected. | Dormand–Prince family polynomials, coupled solution/step equilibrium stability, computer-selected pair coefficients and complex-dominant linear examples. Page 32 discusses norms and nonnormality. | The inspected paper does not classify our selected long-lived real-coordinate clouds. The framework and the setting are already known. |
| [Humphries and Christodoulou, *Phase Space Error Control for Dynamical Systems II*, Sussex report 2000/13](https://math.mcgill.ca/humphries/research/papers/psthetarep.pdf). Introduction and Sections 2, 5.2 and 7 inspected. | Section 5.2 connects near-half-turn numerical rotations with irregular-looking oscillations and analyzes DOPRI5(4). | This is close mechanistic prior art for D335. Alternating rosettes are not claimed as a newly discovered mechanism. |
| [Iserles, Peplow and Stuart, 1991, *A Unified Approach to Spurious Solutions Introduced by Time Discretisation. Part I: Basic Theory*](https://epubs.siam.org/doi/10.1137/0728086). Publisher abstract and reference chain; full text not obtained. | General theory of spurious numerical fixed points, period-two solutions and invariant curves, with bifurcations as the discretization step varies. | Not an identified exact match to this adaptive map; broader numerical bifurcation theory must be credited. |
| [Govorukhin, Tsybulin and Karasözen, 2001, *Dynamics of numerical methods for cosymmetric ordinary differential equations*](https://avesis.metu.edu.tr/yayin/5f4593c4-7ec2-4f63-af78-88ecb6c8b341/dynamics-of-numerical-methods-for-cosymmetric-ordinary-differential-equations), DOI 10.1142/S0218127401003504. Institutional abstract and author publication list; not full text. | Numerical methods can create spurious dynamics and transitions to chaos in cosymmetric systems. | A different model and discretization context; it rules out presenting numerical creation of chaos as a new broad phenomenon. |
| [Ranocha and Giesselmann, 2024, *Stability of step size control based on a posteriori error estimates*](https://arxiv.org/html/2307.12677v3). Introduction and Sections 2–3.2 inspected in the shared audit. | Explicit solution/step dynamics and I/PI/PID stability analysis for complex linear test equations. | Their limiter and studied pairs differ. The methodology is close prior art, not a new framework here. |
| [SciPy RK45 documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html) and the archived installed SciPy 1.18.1 source. | Existing Dormand–Prince 5(4) method and acceptance/rejection implementation. | Our candidate geometry is an observation about this existing algorithm, not a proposed solver. |

The 1990 inspection extends the earlier RK23 review, which had covered only
the first two scanned pages. The frozen RK23 packet is left intact; the new
reading scope and file hash are retained in this packet.

## Narrow contribution under investigation

The test problem is a linearly damped spring. With the standard original-axis
component scales and RMS norm, its numerical state can exhibit persistent
folded clouds near the requested absolute tolerance. D335 is especially
interesting because the two visible rosettes belong to a single orbit, and
two independent expanding directions survive long runs, high precision and
unmodified-library replay. D448 supplies a different, strongly folded visual
case with one positive finite rate.

A potentially useful contribution is a reproducible characterization of these
specific regimes, accompanied by an exact accepted-solution bound, declared
search outcomes and careful controls. **The evidence does not yet establish
that these regimes, or equivalent reduced maps, were previously unknown.**
Changing the pair from RK23 to RK45, choosing new parameters, or giving a
shape a nickname does not by itself establish novelty.

## Equivalences deliberately excluded from novelty claims

- Scaling absolute tolerance together with solution coordinates is an exact
  equivalence of the normalized map.
- Reflection, exchanging two orbit bands and taking a return map describe
  the same underlying sampled dynamics; they are not independent discoveries.
- A rigid camera, color map, projection or nonlinear visual deformation cannot
  create a new attractor. The published portraits use actual `(q/atol,p/atol,h)`
  coordinates, equal units and rigid cameras.
- Diagonalizing the ODE while transforming the whole error rule is a change of
  coordinates. Replacing the norm instead changes the algorithm; this control
  is explicitly labeled and the known norm dependence is credited.
- An accepted-step map is different from a physical flow. Numerical
  hyperchaos here is not chaos of the damped spring.

No global conjugacy, invariant-set inequivalence, topological classification
or certified full-state attractor has been established. The two-positive
finite spectrum is a useful discriminator, not a proof of topological novelty.

## Search coverage and gaps

Additional queries included `Dormand Prince chaotic linear`, `Dormand-Prince
stepsize equilibrium stability`, `Runge-Kutta chaotic stepsize`, `adaptive
time stepping bifurcation`, `step-size border collision`, `Runge-Kutta
spurious chaos linear`, `stepsize Lyapunov`, `adaptive integration hyperchaos`,
`Dormand-Prince hyperchaos`, and `stepsize chaotic attractor`. Many hits
concern integrating an already chaotic ODE, which is a different claim.
The bounded search did not identify an exact match; that is not a certificate
of absence.

Material gaps remain: Hall (1986) full text, Wain's 1998 Dundee thesis,
Söderlind (2002) full text, the full Iserles–Peplow–Stuart and cosymmetric
papers, Humphries's 1993 Bath thesis (the author-linked PDF returned 403),
N. S. Christodoulou's Sussex thesis *Fixed and Variable Time-Stepping
Numerical Methods for Dynamical Systems* (title identified in the
[author's publication entry](https://www.christodoulou.edu.gr/Publications.html),
full text not obtained), and exhaustive older numerical-dynamics book coverage. Stuart and Humphries's
1994 stability review was inspected in its introduction and Sections 5 and 9;
it supplies additional prior context for the solution-bound calculation. Previously
denied endpoints were not bypassed. Equivalent reduced maps remain
unclassified. Use **“numerically supported RK45 controller candidates”**,
not **“definitely unprecedented attractors.”**
