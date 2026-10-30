# A compact global attractor and its exact density lift

For the fixed-parameter motor model in the [PDE lifting note](motor-pde-lifting.md), the finite first-harmonic ODE has a unique compact global attractor. The full density PDE has exactly its first-harmonic lift as its compact global attractor in every periodic $H^s$ phase space with $s\geq0$. All complete bounded trajectories have zero omitted-mode residual.

This is an elementary consequence of the exact closure, dissipative estimates and translated residual formula already established in that note. It is not a mathematical-priority claim. A global attractor can contain equilibria, unstable periodic orbits and connecting trajectories; its existence proves neither chaos nor attraction to any particular periodic orbit. No exponential attraction rate for the finite global attractor is asserted.

## Fixed assumptions and phase spaces

Fix a finite positive integer $N$, real $A,B,\gamma$, and $\rho>0$. The source assumptions are a constant total switching rate, first-harmonic forcing and potential, and a spatially uniform row velocity. Time below is the normalized time $\tau$. The normalized equations, signs and density convention come from the [source derivation and exact lifting](motor-pde-lifting.md), based on [Alouges et al., Section 3.1](https://arxiv.org/pdf/2409.03506v2).

Put $z_j=C_j+iS_j$, $a=A+iB$ and $F=|a|=\sqrt{A^2+B^2}$. With real displacements constrained by $\sum_jX_j=0$, the finite equations are

```math
z_j'=(-1+iV_j)z_j+a,\qquad
X_j'=V_j=-\rho X_j-\gamma(S_j-\overline S),\qquad
\overline S=\frac1N\sum_jS_j.
```

Denote this $3N-1$ dimensional real phase space by $Y$ and its forward flow by $\phi_\tau$. The density equations are

```math
\partial_\tau q_j+V_j\partial_\theta q_j
=-q_j+\frac12+A\cos\theta+B\sin\theta,
\qquad \theta\in\mathbb T.
```

For a fixed $s\geq0$, let

```math
E_s=(H^s(\mathbb T;\mathbb R))^N
\times\{X\in\mathbb R^N:\textstyle\sum_jX_j=0\}.
```

Use a Fourier-defined $H^s$ norm and the product of density norms and the Euclidean displacement norm. Translations are isometries and are strongly continuous in these spaces. Write $\Pi(q,X)=(C,S,X)$ for first-harmonic projection, and include the displacement in the lift

```math
\mathcal L(y)=\left(\left[\frac12+C_j\cos\theta+S_j\sin\theta\right]_{j=1}^N,X\right).
```

The finite norm $\|y-\widetilde y\|_{f,s}$ is the norm induced by $\mathcal L(y)-\mathcal L(\widetilde y)$ in $E_s$. It is an ordinary finite-dimensional norm. The statements concern the specified sinusoidal model and do not extend automatically to spatially varying total rates or feedback from additional force harmonics.

## A compact positively invariant absorbing set

Variation of constants, with the mechanical phase having modulus one, gives for all forward times

```math
|z_j(\tau)|\leq e^{-\tau}|z_j(0)|+(1-e^{-\tau})F.
```

In particular, the first harmonics stay bounded. The mechanical variables are stable linear filters of these bounded amplitudes. This prevents finite-time forward escape and supplies a unique global finite flow.

Set

```math
R=F+1,\qquad M=\frac{2|\gamma|R}{\rho}+1,\qquad
D=\{y\in Y: |z_j|\leq R,\ |X_j|\leq M\text{ for all }j\}.
```

This is a nonempty compact set. On a harmonic boundary $|z_j|=R$,

```math
\frac12(|z_j|^2)'
=-|z_j|^2+\operatorname{Re}(\overline z_j a)
\leq -R^2+RF=-R<0.
```

Within the harmonic disks, $|S_j-\overline S|\leq2R$. On either displacement boundary $X_j=\pm M$, its outward derivative is bounded above by

```math
-\rho M+2|\gamma|R=-\rho<0.
```

Also $\sum_jX_j'=-\rho\sum_jX_j$, so the zero-sum constraint is preserved. These bounds make $D$ positively invariant, including at simultaneous boundary contacts.

Absorption is uniform for every bounded subset $B_0\subset Y$. To see this explicitly, choose bounds $Z_B\geq\max_j|z_j(0)|$ and $U_B\geq\max_j|X_j(0)|$ on $B_0$. After

```math
\tau_z=\log\max(1,Z_B)
```

all harmonics lie in the radius-$R$ disks. Before and after that time, the mechanical amplitude has the uniform bound

```math
\widetilde U_B=U_B+\frac{2|\gamma|}{\rho}\max(Z_B,F).
```

For $u\geq0$ measured from $\tau_z$, the stable-filter estimate gives

```math
\max_j|X_j(\tau_z+u)|
\leq\frac{2|\gamma|R}{\rho}+e^{-\rho u}\widetilde U_B.
```

Consequently every trajectory from $B_0$ enters $D$ by time

```math
\tau_z+\rho^{-1}\log\max(1,\widetilde U_B)
```

and stays there. These constants are convenient bounds, not optimized estimates.

## The finite global attractor and sharper complete-orbit bounds

Because $D$ is compact and positively invariant, its forward images are nonempty compact nested sets. Define

```math
\mathcal A_f=\bigcap_{\tau\geq0}\phi_\tau(D).
```

The intersection is nonempty and compact. It is strictly invariant: $\phi_u(\mathcal A_f)=\mathcal A_f$ for every $u\geq0$. For the reverse inclusion, preimages of a point in successively later images of $D$ have a convergent subsequence in $D$; continuity gives a preimage in every forward image, hence in $\mathcal A_f$.

The nested images approach their intersection uniformly in distance. Otherwise points in progressively later images staying a fixed distance from the intersection would have a subsequential limit in $D$ belonging to every image, a contradiction. Combining this with uniform absorption proves

```math
\sup_{y\in B_0}\operatorname{dist}_{f,s}(\phi_\tau y,\mathcal A_f)
\longrightarrow0
```

for every bounded $B_0\subset Y$. Thus $\mathcal A_f$ is the unique compact global attractor, with the usual definition requiring compactness, strict invariance and uniform attraction of bounded sets.

It is also the set of points on complete bounded finite trajectories. In one direction, strict invariance and compactness extend each trajectory through points of $\mathcal A_f$ to all real times. In the other direction, apply the uniform absorption estimate to the bounded range of a complete trajectory, starting arbitrarily far in the past, and then use the intersection formula.

For any complete bounded trajectory, the forward amplitude estimate applied from time $-T$ to time zero gives

```math
|z_j(0)|\leq e^{-T}|z_j(-T)|+(1-e^{-T})F\longrightarrow F.
```

Applying the mechanical filter estimate from the past then gives

```math
|z_j(\tau)|\leq F,\qquad
|X_j(\tau)|\leq\frac{2|\gamma|F}{\rho}
\quad\text{for every real }\tau.
```

These sharper bounds hold throughout $\mathcal A_f$. They do not assert that every forward trajectory enters the sharper box in finite time.

## Exact product coordinates for the PDE

Let $R_j=q_j-\mathcal L_j(\Pi(q,X))$. This residual contains the non-equilibrium mean and all higher harmonics, with zero first harmonic. The lifting calculation gives

```math
R_j(\theta,\tau)=e^{-\tau}
R_j(\theta-X_j(\tau)+X_j(0),0).
```

Introduce the moving residual

```math
W_j(\theta)=R_j(\theta+X_j).
```

The coordinate map $\Theta:E_s\to Y\times\mathcal R_s$, where $\mathcal R_s$ is the closed zero-first-harmonic residual subspace, is a homeomorphism. Its inverse is

```math
q_j(\theta)=\frac12+C_j\cos\theta+S_j\sin\theta+W_j(\theta-X_j).
```

Continuity follows from continuity of the finite projection and joint strong continuity of translations in $H^s$. Both maps take bounded sets to bounded sets because translations preserve norms and the first-harmonic lift is finite dimensional. No Fréchet differentiability of Eulerian translation with respect to $X$ on the same $H^s$ space is required.

In these coordinates the complete forward evolution is exactly

```math
\Theta\Psi_\tau\Theta^{-1}(y,W)=(\phi_\tau y,e^{-\tau}W).
```

This identity constructs the unique global continuous mild PDE semiflow $\Psi_\tau$ on $E_s$. The physical admissibility restriction discussed below is a forward-invariant subset of this space; its image under $\Theta$ need not be a Cartesian product.

## Global attractor, invariant sets and basins in the PDE

For every nonempty compact finite set $K\subset Y$, Fourier orthogonality gives the exact identity

```math
\operatorname{dist}_{E_s}(\Psi_\tau(q,X),\mathcal L K)^2
=\operatorname{dist}_{f,s}(\phi_\tau\Pi(q,X),K)^2
+e^{-2\tau}\|R(0)\|_{(H^s)^N}^2.
```

The first-harmonic difference is orthogonal to the residual, including its mean, in the Fourier $H^s$ norm. Taking the infimum over $K$ leaves the residual term unchanged. This identity, not continuous conjugacy alone, supplies uniform attraction for bounded sets: both the finite projection and the initial residual norms are uniformly bounded on every bounded subset of $E_s$.

It follows that

```math
\mathcal A_{\rm PDE}=\mathcal L(\mathcal A_f)
```

is compact, strictly invariant and uniformly attracts every bounded subset of $E_s$. It is the unique compact global attractor of the PDE. As a set of smooth first-harmonic densities with their displacements, it is the same lift for every $s\geq0$.

More generally, every compact **strictly invariant** PDE set has $R=W=0$. Indeed, if $\Psi_\tau K_*=K_*$ and the residual norms on $K_*$ are bounded by $C$, then every point of $K_*$ has a preimage in $K_*$ for every $\tau>0$, so its residual norm is at most $e^{-\tau}C$ and must vanish. Positive invariance alone does not give this conclusion. The same past-limit argument shows that every complete bounded PDE trajectory has zero residual and projects to a complete bounded finite trajectory.

Define a convergence basin by

```math
\mathcal B_f(K)=\{y:\operatorname{dist}_{f,s}(\phi_\tau y,K)\to0\},
\qquad
\mathcal B_{\rm PDE}(\mathcal L K)
=\{u:\operatorname{dist}_{E_s}(\Psi_\tau u,\mathcal L K)\to0\}.
```

For compact invariant attracting sets this is the usual basin of attraction. The distance identity also applies when $K$ is merely a nonempty compact target set and proves

```math
\mathcal B_{\rm PDE}(\mathcal L K)=\Pi^{-1}\mathcal B_f(K).
```

Only the residual term has an unconditional exponential decay rate. The finite attraction term can have a different, unspecified rate.

## Physical admissibility

When $F\leq1/2$, the forcing satisfies $0\leq f(\theta)\leq1$. The characteristic formula then preserves $0\leq q_j\leq1$, interpreted almost everywhere for rough data. This defines a closed positively invariant physical subset of $E_s$; no interior of the positive cone in $L^2$ is assumed.

The sharper complete-orbit bound makes every point of the global attractor physical:

```math
\frac12-F\leq\mathcal L_j(y)(\theta)\leq\frac12+F,
\qquad y\in\mathcal A_f.
```

Therefore the same $\mathcal L(\mathcal A_f)$ is the global attractor of the physically restricted semiflow. The basin identity restricts by intersecting its right-hand side with the physical phase space. Arbitrary transient finite projections of physical densities need not themselves have physical first-harmonic lifts; that issue does not affect the complete-orbit or attractor conclusion.

## Exact algebra controls and scope

The [boundedness controls](../results/motor-pde-lifting/global-attractor-algebra.py) verify the transport cancellation in the amplitude balance, the sum-of-squares form of the forcing bound, inward boundary expressions, scalar comparison equations and preservation of the displacement constraint. They use exact rational multivariate polynomials; the [record](../results/motor-pde-lifting/global-attractor-algebra.json) stores the source and prerequisite-note hashes.

```powershell
python results/motor-pde-lifting/global-attractor-algebra.py --output outputs/global-attractor-algebra-recheck.json
```

Choose a new output path. These identities supplement the analytic argument above; they are not a numerical PDE truncation or a replacement for the function-space, compactness and complete-trajectory arguments. No new trajectories, chaotic-attractor certificate or scientific-priority claim are produced.
