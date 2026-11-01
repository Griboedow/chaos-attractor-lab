# Attractor portraits

[Selected results](../README.md) · [Research index](../README.md#research-index) · [Figure provenance](../tools/README-portraits.md)

The gallery shows finite numerical trajectories and one separately certified periodic attractor. Each entry links to its evidence and limitations. Different pictures do not by themselves establish different dynamical classes.

All portraits use the same background, palette, dimensions and caption convention. Higher-dimensional systems are identified as projections; explicit physical unit changes are stated. Dates retain the earliest recorded portrait of a study. A later rendering is dated separately when it falls on another day.

Labels such as `Case 177` identify saved parameter sets in the linked reports.
They are identifiers, not measured quantities; physical parameter labels name
the quantity being varied.

## Constructed systems

### Icosahedral Gaussian flow

![Icosahedral Gaussian flow](assets/portraits/icosahedral.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Cases 177 and 183 are three-state autonomous flows. Long integrations and independent solvers support chaotic motion; absorbing balls and local equilibria are interval certified. The directed geometric skeleton has a known integrable counterpart.

Full three-state coordinates; archived orthographic cameras.

[Report](icosahedral-flow.md)

### Tetrahedral Gaussian flow

![Tetrahedral Gaussian flow](assets/portraits/tetrahedral.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Cases 421 and 565 have a stable central equilibrium and sustained irregular outer motion. Three-dimensional tetrahedral chaos is established prior work.

Full three-state coordinates; archived orthographic cameras.

[Report](tetrahedral-flow.md)

### Localized Gaussian vortices

![Localized Gaussian vortices](assets/portraits/gaussian-vortex.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Case 77 survives a 30,000-unit spectrum and independent integration. A competing periodic orbit and resemblance to known threefold geometries limit interpretation.

Full three-state coordinates; archived orthographic camera.

[Report](gaussian-vortex.md)

### Dissipative three-axis twist map

![Dissipative three-axis twist map](assets/portraits/twist.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Cases 1709 and 3254 retain positive finite-time growth under long and high-precision controls. The first case also has a regular regime.

Full three-state coordinates; archived orthographic cameras.

[Report](twist-geometry.md)

## Numerical algorithms

### RAAR between spatial circles

![RAAR between spatial circles](assets/portraits/raar.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Case 1051 has a certified attracting invariant set with entropy at least log(2)/49 per iteration, containing the symbolic set. Transitivity and equality with this finite portrait remain unresolved. Case 1483 is a numerical companion.

Full spatial state (x, y, z); archived orthographic cameras.

[Trapping report](raar-trapping.md) · [Fold in the attracting set](raar-critical-membership.md) · [Symbolic certificate](raar-symbolic-dynamics.md)

### Scalar iteration with acceleration memory

![Scalar iteration with acceleration memory](assets/portraits/anderson.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Plain scalar iteration is globally convergent for the two selected maps. Regularized three-history acceleration instead produces persistent irregular memory states, checked at 80 and 110 digits.

Three-history state; archived orthographic cameras.

[Report](anderson-memory.md)

### RK45 step-size controller

![RK45 step-size controller](assets/portraits/rk45.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

The physical linear oscillator decays. The adaptive solver has a persistent irregular internal state near its error tolerance, reproduced with the original library and high-precision arithmetic.

Position / atol, velocity / atol, proposed step; archived orthographic cameras.

[Report](rk45-oscillator.md)

### RK23 step-size controller

![RK23 step-size controller](assets/portraits/rk23.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

The physical linear oscillator decays. The adaptive solver has a persistent irregular internal state near its error tolerance, reproduced with the original library and high-precision arithmetic.

Position / atol, velocity / atol, proposed step; archived orthographic cameras.

[Report](rk23-oscillator.md)

### Smooth Yogi memory

![Smooth Yogi memory](assets/portraits/yogi-tanh.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

The published tanh variant produces two positive finite growth rates in selected scalar-quadratic regimes. Independent high-precision and TensorFlow controls support the calculations.

Normalized optimizer state (x, m, v); archived orthographic cameras.

[Report](yogi-tanh.md)

### Yogi parameter-plane candidate

![Yogi parameter-plane candidate](assets/portraits/yogi-arc.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Case 670 is retained after five million steps and independent controls. Case 125 has stronger evidence for a second expanding direction; their rates remain finite-time quantities.

Normalized optimizer state (x, m, v); archived orthographic camera.

[Report](yogi-parameter-plane.md)

### Yogi memory

![Yogi memory](assets/portraits/yogi.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Cases 44 and 7590 are visually distinct numerical regimes of the same published optimizer. Initialization and precision controls are retained in the report.

Normalized optimizer state (x, m, v); archived orthographic cameras.

[Report](yogi-geometry.md)

### Yogi: additional expanding regimes

![Yogi: additional expanding regimes](assets/portraits/yogi-expansion.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Cases 4214 and 125 have numerical evidence for two expanding directions. Window lengths, library differences and branch discontinuities are reported separately.

Normalized optimizer states; archived preview cameras. Panels have separate scales.

[Report](yogi-parameter-plane.md) · [Related result](yogi-geometry.md)

### Adamax memory

![Adamax memory](assets/portraits/adamax.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

Four-million-step, high-precision and PyTorch controls support the selected regimes. Case 4099 coexists with two reflected periodic orbits.

Normalized optimizer state; archived orthographic cameras.

[Report](adamax-geometry.md)

### AdaBelief memory

![AdaBelief memory](assets/portraits/adabelief.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

A scalar quadratic objective supports irregular full optimizer states at the selected memory parameters. Finite-time numerical evidence is separated from convergence and novelty claims.

Normalized state (x, m, s); equal displayed coordinate units.

[Report](adabelief-quadratic.md)

## Physical and biological models

### Four motor rows

![Four motor rows](assets/portraits/motor-four.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

At fixed mechanical parameters, b=1.5 gives a narrow braid and b=1.545 a broader regime. Long full-spectrum integrations support chaos in both; the displayed geometry is a projection of an eleven-dimensional flow.

Three orthonormal displacement contrasts; final 1,200 saved time units; common scale.

[Report](motor-four-b-study.md)

### Five motor rows

![Five motor rows](assets/portraits/motor-five.png)

<!-- drawing-date -->
*First recorded portrait: **2026-11-01**.*
<!-- /drawing-date -->

The same saved fourteen-dimensional trajectory is shown in two projections of its four displacement contrasts. Neither view defines a three-dimensional embedding.

Fourier displacement contrasts; final 800 saved time units.

[Report](motor-five-geometry.md)

### Three motor rows with unequal rates

![Three motor rows with unequal rates](assets/portraits/motor-three-asymmetric.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-31** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

Case 2566 is an eight-dimensional numerical candidate. The projection uses a fixed orthogonal rotation of the three cosine occupancies.

Rotated cosine occupancies; final 600 saved time units; stride 2.

[Report](motor-asymmetry-geometry.md)

### Eight motor rows: periodic and irregular motion

![Eight motor rows: periodic and irregular motion](assets/portraits/motor-eight.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-30** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

A five-cluster periodic attractor is certified in all 23 state dimensions. The irregular comparison has six positive finite-time exponents over 20,000 units; it is not a certified chaotic attractor.

Physical coordinates (X0, X1, s0) of the 23D model; common scale.

[Report](motor-activity-continuation.md) · [Related result](motor-five-cycle-certificate.md)

### Three symmetric motor rows

![Three symmetric motor rows](assets/portraits/motor-three.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-30** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

Case 783 has an independently computed eight-dimensional spectrum with one positive and one near-zero rate. This image shows the three sine occupancies.

Sine occupancies (s0, s1, s2); a projection of the 8D flow.

[Report](motor-rows-gpu-survey.md)

### Motor-pair model

![Motor-pair model](assets/portraits/motor-pair.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-30** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

Case 444 retains positive growth in a long independent CPU calculation. The neutral rate still needs convergence; priority and asymptotic attraction remain unresolved.

Physical coordinates (h, Y, Z); final 2,000 saved time units.

[Report](motor-pair-gpu-survey.md)

### Synchronematic oscillators

![Synchronematic oscillators](assets/portraits/synchronematic.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-29** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

The four-particle candidate is displayed through gauge-invariant collective observables. These do not form an embedding of the seven-dimensional reduced angular dynamics.

Collective observables (Re Q, Im Q, C), not state coordinates.

[Report](synchronematic-research.md)

### Light-driven elastomer slider

![Light-driven elastomer slider](assets/portraits/lce-slider.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-29** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

The retained regime belongs to a three-state relay model. The report separates transverse growth estimates, switching assumptions and known related mechanisms.

Scaled full state (u/q, v/(sqrt(k) q), phi/q); scaling is explicitly physical.

[Report](lce-slider-research.md)

### Two-strain epidemic model

![Two-strain epidemic model](assets/portraits/two-strain.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-30** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

The numerical five-state regime contains unequal outbreaks and extremely small prevalence. Related epidemic chaos is established; physical interpretation requires the qualifications in the report.

Projection (S, 1000 I1, 1000 I2); explicit axis units, not an embedding.

[Report](two-strain-research.md)

### Coupled optical operators

![Coupled optical operators](assets/portraits/optical-relay.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-29** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

These six-state relay histories have persistent irregular motion in the recorded observation. Later continuations, captures and source-model assumptions are retained in the report.

Projection (x_minus, x_plus, y_minus - y_plus); all saved adaptive event/QR knots. The knots are nonuniform in time. Point density is sampling density, not an invariant measure.

[Report](coupled-operators-research.md)

### Smooth optical feedback

![Smooth optical feedback](assets/portraits/optical-smooth.png)

<!-- drawing-date -->
*First recorded portrait: **2026-10-29** · Current rendering: 2026-11-01.*
<!-- /drawing-date -->

At the same logistic-width parameters, saved starts produce regular and irregular histories. This is a finite coexistence observation, not an asymptotic attractor classification.

Projection (x_minus, x_plus, y_minus - y_plus); common scale.

[Report](coupled-operators-research.md)
