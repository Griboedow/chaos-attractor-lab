# What would prove the motor switching network?

The [postcritical trajectories](motor-postcritical.md) suggest repeated visits to symmetry-related six-cluster periodic saddles at `rho=3.268`. The finite connections and long switching trajectories are numerical evidence. This note identifies the missing global theorem and the role of the exact invariant subspaces. It does not certify a heteroclinic network or its attraction.

The subsequent [280,000-unit logarithmic continuation](motor-log-seventh-departure.md) resolves seven departures in two numerical settings, including recovery after a contrast falls below `10^-486`. Their seventh departure times differ by about `4.27e-6`. This longer history strengthens the finite switching evidence, while leaving the global existence and attraction problems below unchanged.

## Exact invariant subspace and local dimensions

Use zero-based row labels and write each row state as `q_j=(c_j,s_j,X_j)`, with the total displacement constrained to zero. The full state dimension is 23. The set

```math
I_{23}=\{q_2=q_3\}
```

is invariant under the polynomial motor field: the two rows have identical equations and identical mean input. It has dimension 20. At the reference six-cluster orbit, pairs `(2,3)` and `(6,7)` are synchronized. The [subsequent interval certificate](motor-saddle-certificate.md) proves one expanding direction, which splits `(6,7)`; all three contrasts in `(2,3)` contract.

The observed destination is a row-permuted copy of the same orbit. Its expanding pair is `(2,3)`, while its other synchronized pair is `(4,5)`. Restricting to `I_23` therefore removes different spectral blocks at the two ends:

| Periodic orbit | Removed block | Spectrum remaining in the 20D invariant subspace |
| --- | --- | --- |
| Source | Three stable contrasts of pair `(2,3)` | One expanding, one flow, eighteen stable |
| Destination | One expanding and two stable contrasts of pair `(2,3)` | One flow, nineteen stable |

These dimension counts now follow from the certified full block spectra. They do not themselves establish a connection.

By those strict spectra, the source is a hyperbolic periodic saddle in `I_23`, and the destination is an orbitally attracting periodic orbit in that subspace. The local unstable manifold of the source lies inside `I_23`: restriction to this invariant space has the same unstable direction and flow direction, and uniqueness of the local unstable manifold identifies it with the ambient one. Here the unstable manifold has dimension two, including time along the periodic orbit.

## A concrete remaining existence problem

It would suffice to enclose a point on the **actual unstable manifold** of the source and prove that its forward orbit enters a certified attraction neighborhood of the destination inside `I_23`. A point displaced along a computed unstable eigenvector is only an approximation to that manifold. Integrating such a point close to the destination, even very accurately, does not settle its infinite past or supply a global connection proof.

A possible validated computation has three parts:

1. Bound a local unstable manifold in a Poincaré section, including its nonlinear deviation from the Floquet eigendirection. Its phase dependence must be controlled.
2. Integrate an enclosed outgoing set through the finite transition with controlled wrapping. Starting after the long, nearly periodic residence would avoid validating thousands of unnecessary time units.
3. Prove entry into a quantitative contraction neighborhood of the destination's return map in the 20D invariant space. Stable Floquet eigenvalues alone do not specify the neighborhood's size; derivative bounds in a suitable norm are needed.

The first validated connection would have all of its row-permuted images. To identify the whole observed network, one must still determine which outgoing branches and phases connect to which destinations. The subsequent [exact group calculation](motor-network-combinatorics.md) supplies eight distinct saddle images and sixteen formal directed edges in the generated group, including the observed mixed four-step circuit. The [full permutation extension](motor-saddle-isotropy.md) counts 5,040 labeled saddle copies and 630 formal eight-node components. These exact combinatorial statements establish candidate node and edge labels; an edge still requires an actual connecting-orbit proof.

## Robustness is conditional; attraction is another question

If an unstable-manifold point enters the destination's open basin in `I_23`, that entry persists under sufficiently small perturbations preserving the invariant subspace and the local hyperbolic periodic orbits. This is the relevant source-to-sink mechanism **within the restricted space**. Perturbing individual row parameters can destroy equality subspaces, so arbitrary physical heterogeneity is not covered.

For cycles containing the required unstable manifolds, Ashwin, Karabacak and Nowotny give a sufficient robustness criterion using invariant subspaces in which each destination is attracting. Their hypothesis already includes the connecting manifolds: the criterion cannot be used to infer their existence from the spectral table alone. Their paper also separates robustness from attraction. [Primary article, Section 2.3 and Theorem 3](https://link.springer.com/article/10.1186/2190-8567-1-13).

The [certified leading contraction-to-expansion ratio](motor-saddle-certificate.md) near the saddle is about 1.876033067695236. A ratio above one is suggestive, but it is not a complete stability test for this network. Return maps, other transverse directions, competing outgoing paths and the coverage of the unstable manifolds still matter. Increasing residence times or very small pair differences alone cannot establish an asymptotic attractor.

## Keep the distinct research questions separate

The lower-rho activity-path results expose a different problem. A [five-cluster periodic attractor](motor-five-cycle-certificate.md) is now certified at exact `rho=129/55`. Its [numerical branch](motor-cluster-cycle-branch.md) loses quotient stability through a complex pair near `rho=2.29823502`, while its pair contrasts remain stable at the sampled points. Validating that crossing and computing the nonlinear normal-form coefficient would be concrete next steps toward a local torus theorem. Neither would by itself connect the resulting motion to the separate irregular cluster candidate.

The [one-positive cluster candidate](motor-lane46-spectrum.md) has a full conditional spectrum and small full-space perturbation controls, but still needs evidence distinguishing sustained irregular dynamics from very long transients near invariant cluster states. The [original-start logarithmic diagnostic](motor-original-seed-log.md) demonstrates recovery of tiny differences and disagreement between distant numerical histories. It therefore rules out using a small terminal pair distance as a stand-alone capture test, while leaving the exact initial state's asymptotic outcome open.

For the [six-positive full-space candidate](motor-activity-continuation.md), the priority is an invariant-set or return-map argument that proves asymptotic chaotic behavior and attraction. Further finite-time exponent accumulation alone would not supply that theorem. Scientific priority is a separate source problem: the unresolved exact-model literature gaps remain even if an attracting chaotic set is proved.

## Established precedents and the remaining novelty question

Slow switching among conjugate cluster symmetries is established prior art. Ashwin and Borresen studied it in globally coupled phase oscillators and a phase-energy extension in 2004. [Publisher abstract](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.70.026203). Wordsworth and Ashwin studied a five-oscillator network with two synchronized pairs and one singleton, including a network of thirty related periodic states and sensitivity to input and noise. Their phase coupling differs from the motor field. [Primary author-hosted paper, Sections I–II](https://empslocal.ex.ac.uk/people/staff/pashwin/PAPERS/WorAsh_2008.pdf).

Consequently, even a complete motor-network theorem would establish a result in these specific motor equations, rather than a new general switching mechanism. The [exact-model priority audit](motor-followup-priority.md) still has unread author theses. A bounded follow-up search found no matching exact motor-network report in its retrieved material; that is not evidence that none exists. [Search scope and source record](../results/motor-priority-followup/network-literature-20261030.json).
