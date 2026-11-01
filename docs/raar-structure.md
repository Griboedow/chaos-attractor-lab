# Spatial-circle RAAR: bounds and fixed-point certificates

These are properties of an existing RAAR map, not a new algorithm. The numerical candidates and prior-work limits are described in [the report](raar-circles.md) and [priority review](raar-priority.md). General fixed-point characterizations already exist; the contribution here is an executable classification for six explicit spatial configurations.

## Map and absorbing bound

Let A be the unit circle in the xy plane, and let

\[
B(\gamma)=c+r(\cos\gamma,\cos\theta\sin\gamma,\sin\theta\sin\gamma).
\]

For a circle with center c, radius r and unit normal n, set v=(I−nnᵀ)(x−c). Away from its axis, P(x)=c+rv/||v||. The standard relaxed averaged alternating reflections map is

\[
F(x)=\tfrac\beta2(R_AR_Bx+x)+(1-\beta)P_Bx
=\beta x+\beta P_A(2P_Bx-x)+(1-2\beta)P_Bx.
\]

For 0<β<1, the triangle inequality gives

\[
\|F(x)\|\le\beta\|x\|+C,\quad
C=\beta+|1-2\beta|(\|c\|+r),
\]
\[
\|x_n\|\le\beta^n\|x_0\|+\frac{1-\beta^n}{1-\beta}C.
\]

Thus every ball of radius strictly greater than C/(1−β) absorbs every orbit for which the map is evaluated. The bound also holds for any nearest-point selection on the projection axes. This elementary bounded-set consequence is not claimed as a new general theorem. The implemented map stops within 10⁻¹² of an axis; no smooth extension is invented. A bounded orbit does not alone establish a continuous attractor or chaos.

## Derivative structure

At regular points, a circle projector has derivative (r/||v||)ttᵀ, where t is its unit tangent. Write DP_A(R_Bx)=a uuᵀ, DP_B(x)=b vvᵀ and d=uᵀv. Then

\[
J=\beta I-\beta a uu^T+(1-2\beta)bvv^T+2\beta ab d\,uv^T,
\]
\[
\det J=\beta^2[(1-a)(\beta+(1-2\beta)b)+ab d^2].
\]

The vector u×v is a **left** eigenvector with eigenvalue β when nonzero. The same eigenvalue persists when the tangents are parallel. The other eigenvalues solve λ²−Tλ+D=0, where

\[
T=2\beta-\beta a+(1-2\beta)b+2\beta ab d^2,
\quad D=\beta[(1-a)(\beta+(1-2\beta)b)+ab d^2].
\]

The symbolic audit checks the entire characteristic polynomial. Pointwise eigenvalues do not equal Lyapunov exponents of a product with changing eigendirections. In particular, this identity does not require one exponent to be log β. The geometry audit retains every negative determinant count, minimum observed axis distance and normalized identity residual. Two earlier failed audit versions are included, including the overly strict absolute residual threshold near a projection axis.

## Complete reduction of fixed points to a degree-eight polynomial

At a fixed point, let b=P_Bx and a=P_A(R_Bx), k=β/(1−β). The map equation gives x=b+k(a−b). Projection stationarity implies (a−b) is perpendicular to both circle tangents. Put

\[
X=c_x+r\cos\gamma,\quad Y=c_y+r\cos\theta\sin\gamma,\quad L=\sqrt{X^2+Y^2}.
\]

When L>0, a=σ(X,Y,0)/L, σ∈{−1,1}. The remaining scalar stationarity condition is f+σe/L=0, with

\[
f=-c_x\sin\gamma+(c_y\cos\theta+c_z\sin\theta)\cos\gamma,
\]
\[
e=c_x\sin\gamma-c_y\cos\theta\cos\gamma
+r\sin^2\theta\sin\gamma\cos\gamma.
\]

Set t=tan(γ/2), D₀=1+t², S=2t, Q=1−t². Define

\[
X_0=c_xD_0+rQ,\quad Y_0=c_yD_0+r\cos\theta S,
\]
\[
f_0=-c_xS+(c_y\cos\theta+c_z\sin\theta)Q,
\]
\[
e_0=c_xSD_0-c_y\cos\theta QD_0+r\sin^2\theta SQ.
\]

Every fixed point must arise from a real root of

\[
P(t)=f_0^2(X_0^2+Y_0^2)-e_0^2.
\]

The interval Sturm calculation certifies the total number of roots and isolates every one. It then rejects the spurious squared branch and tests the two necessary and sufficient radial inequalities

\[
(2b-x)^Ta>0,\qquad (x-c)^T(b-c)/r>0.
\]

Both signs are strict at every candidate, so no projection-axis fixed points remain. The omitted γ=π chart point is tested separately. B never intersects the A axis: the interval value of (c_x/r)²+(c_y/(r cos θ))²−1 is nonzero.

| Case | Real polynomial roots | Valid fixed points | Result |
|---|---:|---:|---|
| R75 | 6 | 0 | Fixed-point free |
| R820 | 4 | 0 | Fixed-point free |
| R1051 | 6 | 0 | Fixed-point free |
| R1311 | 8 | 0 | Fixed-point free |
| R1483 | 6 | 0 | Fixed-point free |
| R1679 | 6 | 2 | Both unstable |

For a valid point, the remaining quadratic is unstable if D>1, D<−1, q(1)<0 or q(−1)<0. Strict interval inequalities certify instability of both R1679 points. All six pairs of circles are disjoint, checked by intersecting A with the plane of B and testing both possible points against B's radius equation.

The mathematical parameters are **exact binary64 lifts**, not shortened table decimals; trigonometric values are enclosed at 80 decimal digits using mpmath interval arithmetic. Floating roots only propose rational intervals. Sturm variations independently certify completeness and all signs. The software certificate depends on the correctness of mpmath's outward-rounded arithmetic; it is not a proof-assistant formalization.

## Constraint linking, not attractor topology

The separate interval audit counts signed intersections of B with the oriented disk bounded by A. R1679 has linking number −1; the other five cases have linking number zero. The sign uses increasing γ on B and the +z normal on A's disk. R1679's two plane crossings have squared radii approximately 0.456549 and 4.755745, so exactly one crosses the interior disk. This classifies the **two input circles**, not the observed invariant cloud.

## Reproduce the certificates

```text
python research/certify_raar_fixed_points.py --input results/raar-circles/fixed-proof/input.json --out NEW_FIXED_PROOF_DIRECTORY
python research/audit_raar_algebra.py --input results/raar-circles/algebra/input.json --out NEW_ALGEBRA_DIRECTORY
```

The first requires mpmath 1.3.0; the second also requires SymPy 1.14.0. Fixed-point freedom, boundedness, constraint linking and numerical growth estimates are separate statements. None establishes historical uniqueness or excludes every long-period orbit.
