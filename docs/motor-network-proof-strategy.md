# What would prove the motor switching network?

The [postcritical trajectories](motor-postcritical.md) suggest repeated visits to symmetry-related six-cluster periodic saddles at `rho=3.268`. The finite connections and long switching trajectories are numerical evidence. This note identifies the missing global theorem and the role of the exact invariant subspaces. It does not certify a heteroclinic network or its attraction.

## Exact invariant subspace and local dimensions

Use zero-based row labels and write each row state as `q_j=(c_j,s_j,X_j)`, with the total displacement constrained to zero. The full state dimension is 23. The set

```math
I_{23}=\{q_2=q_3\}
```

is invariant under the polynomial motor field: the two rows have identical equations and identical mean input. It has dimension 20. At the reference six-cluster orbit, pairs `(2,3)` and `(6,7)` are synchronized. The numerical spectrum has one expanding direction, which splits `(6,7)`; all three contrasts in `(2,3)` contract.

The observed destination is a row-permuted copy of the same orbit. Its expanding pair is `(2,3)`, while its other synchronized pair is `(4,5)`. Restricting to `I_23` therefore removes different spectral blocks at the two ends:

| Periodic orbit | Removed block | Spectrum remaining in the 20D invariant subspace |
| --- | --- | --- |
| Source | Three stable contrasts of pair `(2,3)` | One expanding, one flow, eighteen stable |
| Destination | One expanding and two stable contrasts of pair `(2,3)` | One flow, nineteen stable |

These dimension counts follow algebraically once the full block spectra are certified. The existing numerical calculations support the required counts. They do not themselves establish a connection.

Conditional on those strict spectra, the source is a hyperbolic periodic saddle in `I_23`, and the destination is an orbitally attracting periodic orbit in that subspace. The local unstable manifold of the source lies inside `I_23`: restriction to this invariant space has the same unstable direction and flow direction, and uniqueness of the local unstable manifold identifies it with the ambient one. Here the unstable manifold has dimension two, including time along the periodic orbit.

## A concrete remaining existence problem

It would suffice to enclose a point on the **actual unstable manifold** of the source and prove that its forward orbit enters a certified attraction neighborhood of the destination inside `I_23`. A point displaced along a computed unstable eigenvector is only an approximation to that manifold. Integrating such a point close to the destination, even very accurately, does not settle its infinite past or supply a global connection proof.

A possible validated computation has three parts:

1. Bound a local unstable manifold in a Poincaré section, including its nonlinear deviation from the Floquet eigendirection. Its phase dependence must be controlled.
2. Integrate an enclosed outgoing set through the finite transition with controlled wrapping. Starting after the long, nearly periodic residence would avoid validating thousands of unnecessary time units.
3. Prove entry into a quantitative contraction neighborhood of the destination's return map in the 20D invariant space. Stable Floquet eigenvalues alone do not specify the neighborhood's size; derivative bounds in a suitable norm are needed.

The first validated connection would have all of its row-permuted images. To identify the whole observed network, one must still determine which outgoing branches and phases connect to which destinations. The displayed four-node circuit is taken modulo a half-period symmetry; its numerical permutation match should also be upgraded to an exact statement if used to close a rigorous circuit.

## Robustness is conditional; attraction is another question

If an unstable-manifold point enters the destination's open basin in `I_23`, that entry persists under sufficiently small perturbations preserving the invariant subspace and the local hyperbolic periodic orbits. This is the relevant source-to-sink mechanism **within the restricted space**. Perturbing individual row parameters can destroy equality subspaces, so arbitrary physical heterogeneity is not covered.

For cycles containing the required unstable manifolds, Ashwin, Karabacak and Nowotny give a sufficient robustness criterion using invariant subspaces in which each destination is attracting. Their hypothesis already includes the connecting manifolds: the criterion cannot be used to infer their existence from the spectral table alone. Their paper also separates robustness from attraction. [Primary article, Section 2.3 and Theorem 3](https://link.springer.com/article/10.1186/2190-8567-1-13).

The measured contraction-to-expansion ratio near one saddle is about 1.876. A ratio above one is suggestive, but it is not a complete stability test for this network. Return maps, other transverse directions, competing outgoing paths and the coverage of the unstable manifolds still matter. Increasing residence times or very small pair differences alone cannot establish an asymptotic attractor.

## Established precedents and the remaining novelty question

Slow switching among conjugate cluster symmetries is established prior art. Ashwin and Borresen studied it in globally coupled phase oscillators and a phase-energy extension in 2004. [Publisher abstract](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.70.026203). Wordsworth and Ashwin studied a five-oscillator network with two synchronized pairs and one singleton, including a network of thirty related periodic states and sensitivity to input and noise. Their phase coupling differs from the motor field. [Primary author-hosted paper, Sections I–II](https://empslocal.ex.ac.uk/people/staff/pashwin/PAPERS/WorAsh_2008.pdf).

Consequently, even a complete motor-network theorem would establish a result in these specific motor equations, rather than a new general switching mechanism. The [exact-model priority audit](motor-followup-priority.md) still has unread author theses. A bounded follow-up search found no matching exact motor-network report in its retrieved material; that is not evidence that none exists. [Search scope and source record](../results/motor-priority-followup/network-literature-20261030.json).
