# An exactly straight geodesic on a surface with fifteen handles

A particle can follow a straight line around a periodic surface with fifteen handles. It returns to its starting point after winding (7, −7, −2) times around the three coordinate axes. An infinitesimal perturbation in the unstable direction grows by a factor between **4.45 and 4.47 million** per lap.

[Flow equations](#three-waves-and-a-unit-speed-flow) · [Surface topology](#why-the-surface-has-fifteen-handles) · [Hyperbolicity](#why-the-geodesic-dynamics-is-hyperbolic) · [Multiplier proof](#the-ray-and-its-multiplier) · [Reproduce](#reproduction)

![An exact straight ray and a generic geodesic on the fifteen-handle surface in one periodic cell](crystal.png)

<!-- drawing-date -->
*Drawn: **2026-11-22**.*
<!-- /drawing-date -->

## Three waves and a unit-speed flow

Position q = (x,y,z) is periodic modulo 2π in each coordinate. The figure uses n = 7 and k = 3/8. Write g and H for the gradient and Hessian of F, v for velocity, and u for its tangential part. Where g is nonzero, the flow is

```math
\begin{aligned}
 F(q)&=\cos x+\cos y+k\cos(x-y+nz),\qquad
 u=v-\frac{g\cdot v}{|g|^2}g,\\
 \dot q&=v,\qquad
 \dot v=-\frac{v^T H v+3g\cdot v+2F}{|g|^2}g+(1-|u|^2)u .
\end{aligned}
```

The figure shows the position surface F = 0. Its unit tangent vectors form the attractor in position–velocity space. On this set the acceleration is normal, so the particle follows a unit-speed geodesic.

To see attraction, differentiate the constraint along the flow. The feedback gives F″ + 3F′ + 2F = 0, so both F and F′ decay. Near the unit tangent bundle the squared speed S = |v|² obeys S′ = 2(1−S)S + O(|F| + |F′|), uniformly since |g| is bounded below. The equilibrium S = 1 is locally stable under this decaying forcing.

## Why the surface has fifteen handles

Set w = x − y + nz. The map (x,y,z) ↦ (x,y,w) is an n-fold torus cover. Project the basic level surface onto (x,y). For 0 < k < 2 its image |cos x + cos y| ≤ k is the torus minus two open disks. The excluded caps enclose the unique maximum and minimum, above and below the saddle value zero. Gluing two copies of the punctured torus along both boundary circles gives the basic surface, with Euler characteristic −4.

For k ≤ 1 the closed curve x = π/2, y = arccos(−k cos w) winds once in w and zero times in x and y. It cycles all n lifts in z, proving connectedness at k = 3/8. Regularity and compactness preserve the covering type throughout 0 < k < 2. The cover has Euler characteristic −4n and genus 2n + 1, hence fifteen handles for n = 7.

## Why the geodesic dynamics is hyperbolic

Write A = cos x, B = cos y and C = k cos w, so A+B+C = 0. In (x,y,w) the curvature numerator (gradient times adjugate Hessian times gradient) is k²AB + AC + BC. Changing back to (x,y,z) multiplies it by n². Dividing by |g|⁴ gives

```math
 K=-\frac{n^2}{|g|^4}\big(A^2+B^2+(2-k^2)AB\big).
```

For 0 < k < 2 the quadratic form is positive definite. Curvature is negative except at the finitely many points where all three cosines vanish. Every unit-time geodesic segment therefore contains an open interval of negative curvature.

In the Jacobi equation J″ = −KJ, the nonnegative (J,J′) cone is forward invariant. Negative curvature on part of each segment sends both boundary rays into its interior. Compactness makes the cone condition uniform. Time reversal gives the stable cone, proving that the unit tangent flow is Anosov. [Eberlein's criterion](https://doi.org/10.4310/jdg/1214431801) · [Xu on compact surfaces of nonpositive curvature](https://www.numdam.org/item/AFST_2015_6_24_3_625_0/)

## The ray and its multiplier

Launch from q₀ = (π/2, π/2, π/14) with unit velocity parallel to (1,−1,−2/7). Along the straight ray w = π/2 and y = π−x, so the first two cosines cancel and F stays zero. The acceleration vanishes. One primitive lap has length L = 2π√102. Its Jacobi coefficient repeats every half turn in x.

```math
 -K=\left(\frac{nk\cos x}{2\sin^2 x+(n^2+2)k^2}\right)^2 .
```

At 128-bit precision, bound this coefficient on each of 16,384 cells in a half turn. With nonnegative Jacobi initial data, the constant-coefficient solutions at the two bounds enclose each entry of the transfer matrix. Its determinant is one by the Wronskian identity. A trace bound t thus gives an eigenvalue bound (t + √(t² − 4))/2. One lap contains 14 half turns. Raising the two eigenvalue bounds to the fourteenth power gives the stated interval for infinitesimal transverse growth. [Interval certificate](ray-certificate.json) · [Verification code](../../../research/wow_search_20261122/lattice/certify.py)

## What was measured

A sweep of 72 winding frequencies and amplitudes selected the fifteen-handle case. For a deflection of 10⁻⁶ radians the particle first moved 0.2 distance units from the straight ray at time 47.454, before one lap (about 63.457). DOP853 and RK45 agreed within 10⁻⁶ time units for launch deflections of 10⁻⁴, 10⁻⁶ and 10⁻⁸ radians.

The gold path uses the first 337.5 time units of a second orbit, launched at a 0.321-radian deflection. In 280,000 positions sampled every 0.1 time units, it visited all 14 sheets (19,365–20,824 samples per sheet). The maximum absolute constraint residual was below 0.000026. Surface shading comes from the exact constraint. [Parameters, initial conditions and controls](control.json)

Kuznetsov studied three coupled rotators on a three-cosine surface with Anosov dynamics. Here the shear gives a fifteen-handle cover with an explicit straight orbit and a certified multiplier. Priority for the weighted cover is open. [Kuznetsov, *Chaos in three coupled rotators*](https://arxiv.org/abs/1708.04493)

A [classification of all straight closed geodesics](../../low_hanging_20261122/geodesics/README.md) counts 24 at these coefficients. Perturbing the second wave weight leaves eight, with the shortest straight lap more than seven times its former length.

## Reproduction

Run from the repository root with NumPy, SciPy, Numba, Matplotlib, Pillow and python-flint. Replace `RUN` with an external output directory.

```text
python research/wow_search_20261122/lattice/geodesic.py --out RUN/lattice --report RUN/lattice/control.json
python research/wow_search_20261122/lattice/certify.py --out RUN/lattice/ray-certificate.json
python research/wow_search_20261122/lattice/figure.py RUN/lattice/control.json --out RUN/lattice/crystal.png
```
