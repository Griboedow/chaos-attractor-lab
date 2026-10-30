# Connectedness and a dimension bound for the global motor attractor

The normalized motor model has a unique equilibrium. When $\rho>0$ and $r>\rho+1$, its unstable manifold has real dimension $2N-2$ and belongs to the compact global attractor. Consequently that attractor has both Hausdorff and topological dimension at least $2N-2$. The global attractor is also connected.

For the two exact eight-row source points considered below, this gives a connected global attractor of dimension **at least fourteen**, containing an equilibrium with **fourteen unstable and nine stable directions**. The [volume identity](motor-volume-contraction.md) simultaneously gives zero ordinary 23-dimensional Lebesgue measure. These statements do not determine an exact dimension, assert a noninteger dimension, or prove chaos.

## Unique rest state and its full linearization

Use the unchanged normalized field on $\sum_jX_j=0$:

```math
V_j=-\rho X_j+r(s_j-\bar s),\qquad
c_j'=1-c_j-V_js_j,\qquad
s_j'=b-s_j+V_jc_j,\qquad X_j'=V_j.
```

Here $N\geq1$ is finite, $r>0$, $\rho>0$, and $b$ is any real number. At equilibrium $X_j'=0$ makes every $V_j=0$, including the dependent final row because $\sum_jV_j=0$. The other equations then give $c_j=1$ and $s_j=b$. The mechanical equation finally gives $X_j=0$. Thus the unique equilibrium is

```math
e_*:\qquad c_j=1,\quad s_j=b,\quad X_j=0\quad(j=0,\ldots,N-1).
```

For perturbations $(u_j,v_j,x_j)$ about this state,

```math
u_j'=-u_j-b\,\delta V_j,\qquad
v_j'=-v_j+\delta V_j,\qquad
x_j'=\delta V_j,\qquad
\delta V_j=-\rho x_j+r(v_j-\bar v).
```

The two common modes $\bar u,\bar v$ each have eigenvalue $-1$. There is no common displacement mode because the displacement sum is constrained. On each of the $N-1$ independent zero-mean row directions, the remaining linearization is the same three-dimensional block

```math
B=\begin{pmatrix}
-1&-br&b\rho\\
0&r-1&-\rho\\
0&r&-\rho
\end{pmatrix},\qquad
\det(\lambda I-B)
=(\lambda+1)\big[\lambda^2+(\rho+1-r)\lambda+\rho\big].
```

The full constrained characteristic polynomial is therefore

```math
\det(\lambda I-Df(e_*))
=(\lambda+1)^{N+1}
\big[\lambda^2+(\rho+1-r)\lambda+\rho\big]^{N-1}.
```

The parameter $b$ affects some matrix entries and eigendirections but not these eigenvalues. Put $d=r-\rho-1$. If $d>0$ and $\rho>0$, the two quadratic roots either are both positive real numbers (positive sum $d$ and product $\rho$) or are a conjugate pair with real part $d/2>0$. Repeated positive roots are also allowed. None lies on the imaginary axis, so the equilibrium is hyperbolic, with

```math
\dim E^u=2(N-1),\qquad \dim E^s=N+1.
```

For $N=1$ the quadratic factor has multiplicity zero and the dimension lower bound is the trivial zero bound. The source applications here have $N=8$.

## Why the unstable manifold is inside the global attractor

The [global-attractor argument](motor-global-attractor.md) applies in normalized variables as well. For $z_j=c_j+is_j$, the forcing amplitude is $\sqrt{1+b^2}$ and

```math
|z_j(\tau)|\leq e^{-\tau}|z_j(0)|
+(1-e^{-\tau})\sqrt{1+b^2}.
```

With

```math
R=\sqrt{1+b^2}+1,\qquad
M=\frac{2rR}{\rho}+1,
```

the set

```math
D=\{(c,s,X):c_j^2+s_j^2\leq R^2,\ |X_j|\leq M,
\ \textstyle\sum_jX_j=0\}
```

is compact, convex, positively invariant and uniformly absorbs every bounded set. The proof is the same amplitude and stable-filter estimate as in the global-attractor note, with forcing $(1,b)$. In particular, every finite initial state has a globally defined bounded forward trajectory. Its compact global attractor $\mathcal A_f$ is exactly the set of points on complete bounded trajectories.

We now use the standard local unstable-manifold theorem for a smooth finite-dimensional ODE at a hyperbolic equilibrium. It supplies a smooth embedded local unstable manifold of real dimension $m=2N-2$. Points on a sufficiently small local patch have trajectories defined for all negative times, remaining near $e_*$ and converging to it as $\tau\to-\infty$.

For any such point, the negative-time trajectory is bounded by this local property, while the positive-time trajectory is bounded by the absorbing estimate. It is therefore a complete bounded trajectory and belongs to $\mathcal A_f$. This proves

```math
W^u_{\rm loc}(e_*)\subset\mathcal A_f
```

for a suitably chosen local patch. The same argument applies to its forward images. No numerical unstable-manifold integration is used to establish this inclusion.

Choose a smaller closed embedded $m$-disk inside that patch. It is compact and hence closed in $\mathcal A_f$. Such a disk has topological dimension $m$, and monotonicity for closed subsets gives the topological lower bound. Its local smooth coordinates can be restricted to a smaller compact chart with bi-Lipschitz bounds, so its Hausdorff dimension is also $m$. Inclusion then gives

```math
\dim_{\rm top}\mathcal A_f\geq2N-2,\qquad
\dim_{\rm H}\mathcal A_f\geq2N-2.
```

These are lower bounds obtained from an actual smooth disk within the global attractor. They do not assert equality or any fractal scaling law. The dynamical theorem invoked is the ordinary finite-dimensional unstable-manifold theorem; no computer-assisted manifold parameterization is claimed.

## Connectedness

The absorbing set $D$ above is convex and therefore connected. Every forward image $\phi_\tau(D)$ is compact and connected by continuity of the flow. Positive invariance makes these nonempty sets nested, and

```math
\mathcal A_f=\bigcap_{\tau\geq0}\phi_\tau(D).
```

An intersection of nested nonempty compact connected sets is connected. Hence $\mathcal A_f$ is connected. This does not establish path connectedness, identify a particular connecting orbit, or determine which locally attracting invariant sets have which basins.

## Exact eight-row source points

At

```math
N=8,\qquad r=1+\frac{19/20}{\sqrt{32001/500000}},\qquad b=1,
```

the [interval parameter controls](../results/motor-pde-lifting/global-geometry-algebra.json) give:

| Source point | Exact $\rho$ | $r-\rho-1$ (rounded) | Rest-state unstable/stable dimensions |
|---|---|---|---|
| Certified six-cluster saddle parameter | $817/250$ | $0.48714604775$ | $14/9$ |
| Eight-row activity-path candidate parameter | $1291/440$ | $0.82105513866$ | $14/9$ |

At both points the quadratic discriminant is strictly negative, so the fourteen unstable real dimensions form seven conjugate eigenvalue pairs. The record encloses their real parts and imaginary magnitudes using Arb and exact source constants. It also verifies the strict weighted-volume condition

```math
a_8=-16+7(r-\rho)<0.
```

Combining the results, the finite global attractor at each point is connected, has Hausdorff and topological dimension at least fourteen, and has zero 23-dimensional Lebesgue measure. Zero measure does not turn the dimension upper bound 23 into a strict upper bound.

The global attractor includes every bounded complete orbit, including this unstable rest state and periodic saddles. A locally attracting cycle or another local attracting set can coexist within it. The geometry proved here does not upgrade the activity-path candidate to a chaotic-attractor theorem.

Under the exact [density lifting assumptions](motor-pde-lifting.md), the normalized first harmonics are related to the physical coefficients by the fixed nonzero source scaling. The density lift is affine and bi-Lipschitz on this finite-dimensional phase space for each fixed $H^s$ norm. Connectedness and these dimension lower bounds therefore transfer to its PDE global attractor. The zero-Lebesgue-measure statement concerns the finite phase space; no infinite-dimensional Lebesgue measure is introduced.

## Reproducible controls

The [exact checker](../results/motor-pde-lifting/global-geometry-algebra.py) verifies the three-dimensional characteristic polynomial symbolically. For $N=1,2,3,8$, it differentiates the complete shifted nonlinear field, evaluates at rest, and checks the entire common-mode/difference-block decomposition with exact rational multivariate polynomials. The basis determinants are nonzero. These implementation controls accompany the general-$N$ decomposition proved above; a finite selection of $N$ is not presented as its proof.

```powershell
python results/motor-pde-lifting/global-geometry-algebra.py --output outputs/global-geometry-algebra-recheck.json
```

The script requires Python-FLINT and a new output path. Its record stores the exact source and supporting-note hashes, all polynomial checks, and 256-bit Arb parameter enclosures. No new trajectory simulations, exact attractor-dimension theorem, chaos certification or scientific-priority claim are produced.
