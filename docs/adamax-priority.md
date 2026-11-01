# Adamax priority audit: what the pictures do and do not establish

Audit date: 1 November 2026. **The displayed crown and serrated wings are well-supported numerical candidates in a known algorithm. This audit does not establish that nobody has seen them, or that they define a new topological attractor class.** The narrower contribution to investigate is the full-state geometry created by maximum-memory resets on a scalar quadratic, with an exact sheet description and reproducible numerical evidence.

## Direct optimizer precedents

| Primary source | Material inspected | Consequence for this study |
| --- | --- | --- |
| Kingma and Ba, *Adam: A Method for Stochastic Optimization*, ICLR 2015, [Section 7.1 and Algorithm 2](https://arxiv.org/pdf/1412.6980v9) | Adamax definition and recommended parameters | Adamax and its maximum-memory recurrence are established. The algorithm is not our invention. Original pseudocode omits epsilon; the exact placement used here comes from the implementation below. |
| [PyTorch 2.15 Adamax documentation](https://docs.pytorch.org/docs/2.15/generated/torch.optim.Adamax.html) and installed implementation | Actual `_single_tensor_adamax` code, version and source hash; short and 300,000-step optimizer runs | The experiment uses a standard implementation with epsilon added to the absolute gradient before taking the maximum. It is not an ad hoc nonlinear replacement of Adamax. |
| Bock and Weiss (2022), [*Non-Convergence and Limit Cycles in the Adam optimizer*](https://arxiv.org/html/2210.02070), Section 4 | Full numerical section, parameter table, fountain description and chaotic example | Scalar-quadratic optimizer chaos and expressive three-state geometry already occur in Adam. Neither broad claim is new. Their square-moment recurrence differs from Adamax. |
| Loh, Mireles James, Xu and Tao (2026), [*The Dynamics of RMSProp Beyond the Stability Threshold*](https://arxiv.org/abs/2610.30078) | Full-text introduction/related work, trapping results, conclusion and Adam extension | Bounded quadratic-optimizer dynamics, invariant circles, tori, coexistence and chaos have recent direct precedents. The inspected full text contains no Adamax or AdaBelief occurrence. This lexical observation is not an exhaustive priority exclusion. |
| Huber and Waclawek (2023), [*Ck-continuous Spline Approximation with TensorFlow Gradient Descent Optimizers*](https://arxiv.org/abs/2303.12454), Section 4 | Full experimental section; model and loss definitions | Recurring instability in Adamax training was already observed, approximately every 500 epochs in their polynomial-fit example. They show loss behavior and move to other optimizers; they do not supply the scalar three-state attractor geometry analyzed here. |
| Sveleba et al. (2023), *The Influence of Sampling Parameters on the Learning Error of a Multilayer Neural Network*, [ELIT proceedings, pp. 51–54](https://onu.edu.ua/pub/bank/userfiles/files/fmvps/kaf_svit_gos_mij_eko_vidnosin/mijnarodni_proekti/2023_IEEE_13th_International_Conference_on_Electronics_and_Information_Technologies.pdf) | Complete four-page article; parameter-plane figures visually inspected | Adamax is already discussed using dynamical-mode, periodicity and transition-to-chaos terminology. Its examples concern a multilayer digit-recognition network and short learning runs. The figures are parameter-plane regime diagrams, not three-dimensional optimizer-state portraits. |

The 2026 RMSProp paper was submitted on 21 October 2026 and is especially relevant to any manuscript written now. This audit reads more than its abstract; the earlier AdaBelief report's abstract-only caveat remains accurate for that earlier packet.

## Equivalent terms and established maximum maps

The search also covered `AdaMax`/`AdamMax`, non-convergence, recurrent instability, limit cycles, bifurcations, Lyapunov exponents, maximum-memory recurrences, max-type difference equations, generalized Lyness equations and Lozi maps. It did not identify a direct published match to the displayed scalar Adamax candidates in the retrieved material. A missing search hit is not evidence of historical absence. The search is bounded, and a complete citation graph, thesis search and mathematical classification remain open.

An important comparison is Linero Bas and Nieves Roldan, [*On the relationship between Lozi maps and max-type difference equations*](https://arxiv.org/abs/2209.03582), DOI [10.1080/10236198.2022.2120393](https://doi.org/10.1080/10236198.2022.2120393). Sections 2 and 5 were inspected. They give logarithmic changes of variables between generalized Lozi recurrences

```text
y[n+1] = A*abs(y[n]) + B*y[n] + C*y[n-1] + D
```

and multiplicative maximum recurrences on positive variables. Their examples explicitly show how different-looking plots and different parameter values can represent conjugate dynamics. Therefore the presence of `max`, sharp corners or unfamiliar silhouettes is not a novelty certificate.

Adamax's signed position, additive first-moment update, positive epsilon and separate maximum-memory state do not directly fit the displayed two-variable logarithmic transformation. Its full map collapses open reset fibers and is not locally one-to-one there. This distinguishes that full-state operation from an invertible Lozi map on an open domain. It **does not** exclude conjugacy of restricted invariant sets or induced return maps; no such exclusion has been proved.

## Exact reductions that prevent artificial novelty

The plotted coordinates are a common scaling of all three physical optimizer variables: `(x,m,u)=(theta,M,U)/alpha` for `f(theta)=theta^2/2`. No anisotropic rescaling or nonlinear display transformation is used. Rigid rotations expose the same state cloud.

For the more general scalar quadratic `h*(theta-theta0)^2/2`, set `x=(theta-theta0)/alpha`, `m=M/(alpha*h)`, `u=U/(alpha*h)`. The normalized recurrence is the same with `e=epsilon/(alpha*h)`. Changing curvature, learning rate and epsilon while preserving this combination does not produce a new normalized system. Translation of the optimum also creates no new dynamics.

The exact sign symmetry identifies reflected states. Reset surfaces and their memory-step images explain much of the visible layering. These facts should be stated before comparing silhouettes or naming apparent substructures.

## Defensible manuscript claim

The present evidence supports a reproducible **case study of full-state chaotic candidates in scalar-quadratic Adamax**, with two visually distinct geometries, arbitrary-precision and actual-library checks, explicit competing regular regimes, and elementary algebraic results describing reset surfaces and bounded cold starts. The [structural note](adamax-structure.md) also proves that a bounded autonomous limit set lies on finitely many transported sheets when epsilon is positive.

The following stronger claims are unsupported: first chaos in an optimizer; first non-convergence in Adamax; a newly invented dynamical system; a unique attractor for every initial state; a proven strange attractor; or a new topological class. The phrase “previously unseen attractor” should remain a research objective, not a conclusion.
