# RK23 controller geometry: priority and equivalence audit

Audit date: 2026-11-01. **Historical uniqueness is unresolved.** The selected
figures are numerical candidates in a published solver, not a new differential
equation or a certified new attractor class. A search that fails to locate the
same figure cannot establish that nobody has described an equivalent system.

## What is already known

| Prior work | Relevant result | Consequence for this study |
| --- | --- | --- |
| [Hall, 1986, *Equilibrium states of Runge–Kutta schemes: part II*](https://doi.org/10.1145/7921.7922) | Extends the equilibrium analysis of adaptive algorithms to complex dominant eigenvalues. Publisher abstract inspected; full PDF access returned 403. | Complex linear test equations and algorithmic equilibrium analysis are established territory. |
| [Hall and Higham, 1988, *Analysis of stepsize selection schemes for Runge–Kutta codes*](https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P1.pdf) | Studies coupled solution/step recurrences, stability-boundary equilibria and a transformed norm. Section 3 explicitly distinguishes normal from nonnormal problems and practical error norms. | Neither the polynomial-map formulation nor the importance of coordinates/error norms is new here. |
| [Higham and Hall, 1990, *Embedded Runge–Kutta formulae with stable equilibrium states*](https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P6.pdf) | Designs embedded 4/5 pairs using the earlier equilibrium theory. The inspected opening pages address smooth step control when stability restricts the step. | Deliberately choosing better controller/solver stability is prior work. This source is adjacent, not an identified match to our RK23 clouds. |
| [Aves, Griffiths and Higham, 1997, *Does error control suppress spuriosity?*](https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P28.pdf) | Shows that error control need not remove stable spurious fixed points, even at arbitrarily small tolerances. | Small tolerance does not make nonphysical numerical dynamics a new phenomenon or physical evidence. |
| [Higham, Humphries and Wain, 2000, *Phase space error control for dynamical systems*](https://webhomes.maths.ed.ac.uk/~dhigham/Publications/P35.pdf) | Describes persistent tolerance-scale numerical oscillations near stable equilibria and develops a phase-space criterion. | Failure to settle at the physical equilibrium is known. That broad observation must not be presented as our discovery. |
| [Humphries and Vigneswaran, 2004, *Phase space stability error control with variable time-stepping Runge–Kutta methods for dynamical systems*](https://www.math.mcgill.ca/humphries/research/papers/PhaseVicky04.pdf) | Extends phase-space/stability control and discusses standard error-per-step controllers, rejections and stable linear systems. Introduction and the displayed standard algorithm inspected. | Adaptive rejection and stable linear examples are not new ingredients. |
| [Ranocha and Giesselmann, 2024, *Stability of step size control based on a posteriori error estimates*](https://arxiv.org/html/2307.12677v3) | Explicitly develops solution/step dynamics for I and PI/PID control and analyzes embedded third-order Runge–Kutta stability on complex linear equations. Their limiter differs from SciPy's policy. | This is close methodological prior art. Our complex-scalar quotient must be compared with it, not marketed as a new framework. |
| [Izgin-Melnikov and Ranocha, 2026, *Using machine learning to design time step size controllers for stable time integrators*](https://doi.org/10.1007/s10915-026-03496-1) | Uses Bayesian optimization to tune DSP controllers for MPRK and Rosenbrock methods. Abstract and introduction inspected; published 2026-10-08. | Automated parameter searches for controllers are also established. This is not an exact method match. |
| [SciPy RK23 documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK23.html) and the archived installed implementation | Implements the Bogacki–Shampine 3(2) pair and its acceptance/rejection policy. | The method, coefficients, factors and error scales are prior work. No new-integrator claim. |

The review distinguishes instability, irregular step histories, spurious cycles
and a positive finite Lyapunov estimate. These are different assertions. The
sources above establish extensive prior work on the underlying phenomenon;
their existence alone is not proof that the particular chaotic invariant sets
in our figures have already been classified.

## The exact question still worth investigating

For the standard real oscillator `q'=p`, `p'=-q-2*zeta*p`, the mathematical
accepted-step map has state `(q/atol, p/atol, proposed h)`. We retain SciPy's
componentwise relative scales, weighted RMS error, rejection loop and factor
limits. Selected parameter pairs produce folded, branched full-state clouds
and positive finite expansion estimates across independent calculations.

The narrow contribution under investigation is a reproducible description of
those specific regimes and their geometry, with exact implementation matching,
high-precision replay, all screened outcomes and a solution-coordinate bound.
Whether that description adds a new dynamical mechanism remains open. A new
parameter value, a new visual nickname, or a more polished render is insufficient.

## Equivalent descriptions that must not be counted twice

1. **Tolerance scaling.** Writing `u=y/atol` removes the absolute tolerance from
   the ideal recurrence. Scaled copies obtained by changing `atol` and scaling
   the initial solution are the same normalized dynamics.
2. **Reflection.** The real map commutes with `(q,p,h)->(-q,-p,h)`. Reflected
   clouds are not separate discoveries.
3. **Complex phase.** A complex scalar test equation using an absolute-value
   error norm has continuous phase symmetry. Its amplitude/step quotient and
   the rotated three-dimensional cloud are two descriptions of one algorithm.
4. **Changing oscillator coordinates.** A linear change that diagonalizes the
   ODE also transforms the error norm and component scales. Retaining the old
   Euclidean/componentwise rule changes the numerical algorithm. Conversely,
   transforming the complete rule is only a conjugate representation.
5. **Clock and warm-up.** The autonomous mathematical map omits the finite
   floating-point clock. Actual SciPy runs preserve it. Agreement of finite
   statistics supports the implementation interpretation; it does not turn the
   floating-point computation into a proof about infinite-time dynamics.
6. **A flow versus an accepted-step map.** The three plotted coordinates include
   a numerical step proposal. This is not a three-dimensional physical flow,
   nor evidence that a linearly damped spring has become chaotic.

The [structure note](rk23-controller-structure.md) gives the precise equation,
policy, normalization and limits of the boundedness certificate. No global
conjugacy claim, invariant-set inequivalence proof, or topological classification
has been established. Appearance and Lyapunov spectra cannot provide those.

## Search coverage and remaining gaps

Searches included `adaptive time stepping chaos`, `stepsize strange attractor`,
`Runge-Kutta controller chaotic linear`, `Bogacki chaotic linear`, `adaptive
Runge-Kutta spurious attractors`, `stepsize chaotic behaviour`, and exact titles
from the Hall–Higham–Humphries reference chain. Searches for algorithmic chaos
in Rayleigh quotient iteration already located older work, so that direction
was not promoted as an unexplored setting.

The bounded review above did not identify an exact published match to the
selected real-coordinate, implementation-specific clouds. **This is a search
result, not a priority certificate.** Several gaps remain material:

- Hall's 1986 full text was inaccessible through the attempted publisher route.
- R. J. Wain's 1998 Dundee thesis, *Long Term Dynamics of Adaptive Algorithms for
  the Numerical Solution of Ordinary Differential Equations*, appears in the
  primary reference chain but its full text was not obtained.
- Söderlind's 2002 *Automatic Control and Adaptive Time-Stepping* was inspected
  through its institutional abstract; the attempted full-text endpoint denied access.
- The broad books and older numerical-dynamics literature have not been checked
  exhaustively, and equivalent reduced maps have not been classified.

Until those gaps and the equivalence problem are resolved, use **"numerically
supported RK23 controller candidates"**, not **"previously unseen attractors"**.
