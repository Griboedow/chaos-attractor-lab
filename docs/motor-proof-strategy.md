# A rigorous route for the motor cluster cycle and its transverse instability

Status: original proof strategy for the follow-up task begun on 30 October 2026. The algebraic reductions below are exact. The earlier [cycle continuation](motor-cluster-transverse.md) supplies numerical starting guesses. The subsequently completed [cycle certificate](motor-cycle-certificate.md) establishes the local periodic-orbit and crossing-existence results, including a least-period extension; this note preserves the design and distinguishes uncompleted extensions. No chaotic attractor is rigorously established.

The first attainable theorem is smaller than a proof of hyperchaos: existence of the periodic four-plus-four orbit, its stability within the cluster subspace, and its sixfold transverse loss of stability. An unstable periodic orbit with six expanding directions does not itself prove a chaotic invariant set or an attractor.

## Fix the mathematical constants before validating

The source parameter prescription used in this repository gives

```math
a_0=c\sqrt{\Omega},\qquad a_1=b_1=d\Omega,\qquad
r=-\frac{\lambda d}{\ell c^2},\qquad
d=-\frac{c\ell}{2\pi\lambda}
\left(2\pi c+\frac{\zeta\ell}{\sqrt{\Omega_0}}\right),
\qquad \zeta=\frac{2\pi k}{\eta\ell}.
```

Consequently the constants involving the interaction potential, length and $\pi$ cancel:

```math
r=1+\frac{k}{\eta c\sqrt{\Omega_0}},\qquad
\rho=\frac{k}{\eta c\sqrt{\Omega_0(1+\delta)}}
=\frac{r-1}{\sqrt{1+\delta}}.
```

Using the declared decimal inputs as exact model constants,

```math
\Omega_0=15(0.0042668)=\frac{32001}{500000},\qquad
\frac{k}{\eta c}=\frac{19}{20},\qquad
(r-1)^2=\frac{451250}{32001},\qquad r\gt1.
```

Thus one can enclose the exact algebraic number $r$ by directed-rounding arithmetic, and recover activity by

```math
\delta=\frac{451250}{32001\rho^2}-1.
```

This is the unrounded parameter prescription implemented in `source_parameters`, not an assertion that a displayed rounded value of $d$ is exact. A theorem should state whether a quoted value of $\rho$ means a particular rational decimal, an interval, or the algebraic value derived from a specified activity. Treating source constants as exact defines a mathematical model; it does not establish the experimental accuracy of those constants. See the [source audit](motor-rows-source-audit.md) for the energy-unit convention.

## Every periodic equal-cluster orbit lies in an exact three-state manifold

Put $C=1-Z$. The core and its two extra cluster coordinates obey

```math
\begin{aligned}
h'&=(-\rho+rC)h-rY,&Y'&=hC-Y,&C'&=1-C-hY,\\
u'&=-u+hv,&v'&=-v-hu+b.
\end{aligned}
```

The identities

```math
u=bY,\qquad v=bC
```

define an invariant manifold. Indeed, for deviations $d=u-bY$ and $e=v-bC$,

```math
d'=-d+he,\qquad e'=-e-hd,\qquad
(d^2+e^2)'=-2(d^2+e^2).
```

If the full cluster orbit is periodic with positive period $T$, periodicity requires $d^2+e^2=e^{-2T}(d^2+e^2)$, hence $d=e=0$. Therefore **every periodic equal-cluster solution is determined by the three-state core alone**. No separate shooting of the two extra coordinates or numerical quadrature for a periodic slave solution is necessary.

The exact physical reconstruction is

```math
\begin{aligned}
c_A&=C+bY,&s_A&=bC-Y,&X_A&=(h-rY)/\rho,&V_A&=-h,\\
c_B&=C-bY,&s_B&=bC+Y,&X_B&=-X_A,&V_B&=h.
\end{aligned}
```

Both row radii satisfy

```math
c_A^2+s_A^2=c_B^2+s_B^2=(1+b^2)(C^2+Y^2).
```

In fact every bounded complete periodic core satisfies this disk bound analytically. Writing $q=C+iY$ gives $q'=1+(-1+ih)q$. Its bounded complete solution has the representation

```math
q(t)=\int_0^{\infty}e^{-s+i\int_{t-s}^{t}h(v)\,dv}\,ds,
\qquad |q(t)|\leq\int_0^{\infty}e^{-s}ds=1.
```

Thus existence of the periodic core automatically gives the stated row-disk bounds. Nonnegative transition rates and reconstructed density bounds still require the source rate-amplitude condition. This is a reduction of the equal-cluster periodic problem; it does not reduce general eight-row trajectories to three states.

There are also two analytically known passive multipliers. Since $X_A'=-h$, a periodic core has $\int_0^T h(t)\,dt=0$. In complex form $d+ie$ evolves by multiplication with $\exp[-T-i\int_0^T h]$. Its monodromy is exactly $e^{-T}I_2$. The five-state cluster spectrum consists of the core spectrum and two copies of $e^{-T}$.

## Validate a half-period boundary problem

In $q=(h,Y,Z)$ coordinates, define $R(h,Y,Z)=(-h,-Y,Z)$. Direct substitution gives $f(Rq)=Rf(q)$. This is a **commuting equivariance**, not a time-reversing symmetry.

For $q_0=(0,y,z)$, solve the three equations

```math
F(y,z,L;\rho)=\phi_L(0,y,z;\rho)-R(0,y,z)=0,
\qquad L\gt0,\quad y\lt0.
```

Equivariance and uniqueness then imply $\phi_{2L}(q_0)=q_0$. The inequalities exclude the equilibrium and give $h'(0)=-ry>0$. To claim **minimal** period $2L$, additionally exclude intermediate zeros of $h$ on each open half-orbit, using validated trajectory bounds and endpoint derivative signs. Existence of a nonconstant periodic orbit does not require that extra minimal-period claim.

Let $M_L=D_q\phi_L(q_0)$. The shooting Jacobian has the explicit columns

```math
D_{(y,z,L)}F=
\left[M_Le_Y+e_Y,\quad M_Le_Z-e_Z,\quad f(\phi_L(q_0))\right].
```

This requires only a validated core flow and its first state variational equation. A fixed unit integration interval can be used with $dq/d\sigma=L f(q)$, keeping uncertain $L$ as a parameter rather than asking an ordinary integrator to stop at an uncertain time.

A simple certificate uses a box $U$ around the numerical $(y,z,L)$, a center $u_0$, and a nonsingular approximate inverse $B$. Enclose

```math
K=u_0-BF(u_0)+(I-B[DF(U)])(U-u_0).
```

Strict inclusion $K\subset\mathrm{int}(U)$ together with a certified contraction bound $\sup\|I-BDF(U)\|\lt1$ proves existence and uniqueness in $U$ by the contraction mapping argument. Every operation and the ODE remainder must be outward enclosed. Ordinary high-precision shooting, a small residual, or two agreeing solvers do not supply this certificate. The CAPD review gives a primary worked periodic-orbit validation using an interval Newton operator [2].

For a family, validate this inclusion uniformly over a parameter interval. Adjacent parameter slabs must be linked by uniqueness domains containing the overlapping solution enclosures. Eighteen separately validated parameter points would not, by themselves, certify a continuous branch between them.

## Validate one transverse block and obtain six copies exactly

For one zero-sum contrast within cluster A,

```math
A(t)=
\begin{pmatrix}
-1&-V_A-rs_A&\rho s_A\\
V_A&-1+rc_A&-\rho c_A\\
0&r&-\rho
\end{pmatrix},\qquad
\Phi_A'=A(t)\Phi_A,\qquad \Phi_A(0)=I_3.
```

All coefficients are explicit functions of the validated three-state core. The full-period matrix $M_A=\Phi_A(T)$ therefore needs only one additional nine-component linear equation.

The certified half-period symmetry gives $A_B(t)=A_A(t+T/2)$. If $U_1$ and $U_2$ are the first- and second-half transfer matrices for A, then

```math
M_A=U_2U_1,\qquad M_B=U_1U_2=U_1M_AU_1^{-1}.
```

Fundamental matrices are invertible, so the transverse spectra coincide exactly. Each group of four has three independent zero-sum contrasts. Hence each eigenvalue of the one three-dimensional block occurs six times in the complete 23-dimensional tangent system. A numerical equality of six eigenvalues is unnecessary once this decomposition and symmetry have been certified.

There is a further determinant identity. Half-period symmetry gives $\int_0^T Y=0$. Since $\mathrm{tr}(A_A)=-\rho-2+r(C+bY)$ and $\mathrm{tr}(Df_{\rm core})=-\rho-2+rC$, Liouville's formula yields

```math
\det M_A=\det M_B=\det M_{\rm core}=d,\qquad
\det M_{23}=d^7e^{-2T}.
```

These provide independent consistency checks and allow the same rigorously enclosed determinant to be reused. The equality with the core determinant uses the verified half-period symmetry; it is not asserted for an arbitrary nonsymmetric core cycle.

The core's autonomous flow direction has multiplier exactly one. Write its other two multipliers as roots of

```math
\lambda^2-a\lambda+d=0,\qquad
a=\mathrm{tr}(M_{\rm core})-1,\qquad
d=\det M_{\rm core}
=\exp\left(\int_0^T[-\rho+r(1-Z)-2]dt\right).
```

Since $d>0$, strict interval inequalities

```math
1-a+d\gt0,\qquad 1+a+d\gt0,\qquad 1-d\gt0
```

certify both remaining core multipliers inside the unit circle. They avoid trying to distinguish a numerically rounded neutral multiplier from one. Together with $e^{-T}<1$, this proves orbital attraction **within the equal-cluster subspace**. The Poincare-map interpretation of the remaining multipliers is discussed in [3].

At the observed crossing the transverse block has three well-separated numerical roots near $-1$, $0.58815$ and $-0.04452$. For a rigorous enclosure, form its real cubic characteristic polynomial using interval matrix invariants. Opposite polynomial signs at the endpoints of three disjoint real intervals place at least one root in each. Degree three then gives exactly one simple root per interval. This can certify the two smaller roots remain inside the unit circle and isolate the critical real root, without a general interval eigenvalue library.

For the earlier stable control, where two transverse roots are complex, one can instead apply a cubic Schur test. If

```math
p(z)=z^3-tz^2+sz-d,
```

the Cayley transform gives

```math
(1-w)^3p\left(\frac{1+w}{1-w}\right)
=a_3w^3+a_2w^2+a_1w+a_0,
```

with $a_0=1-t+s-d$, $a_1=3-t-s+3d$, $a_2=3+t-s-3d$, and $a_3=1+t+s+d$. Positivity of all four coefficients and $a_2a_1>a_3a_0$ places the transformed roots in the left half-plane, hence the original roots inside the unit disk. This is an optional scalar-inequality alternative to certified complex root enclosures.

## A crossing theorem without an unjustified bifurcation theorem

On a rigorously continued local branch, define

```math
g(\rho)=\det(I+M_A(\rho)).
```

Opposite certified signs at two endpoints imply at least one parameter where a transverse multiplier is $-1$. If the remaining roots stay separated from the unit circle throughout the interval, this isolates the possible local stability change. It does **not** establish a unique critical parameter or a nonzero crossing speed.

To prove uniqueness and transversality, enclose $g'$ away from zero. The core family derivative follows from

```math
F_u u_\rho+F_\rho=0,\qquad f_\rho=(-h,0,0)^T
```

at fixed $r$. In a unit full-period coordinate $\sigma$, with $T=2L$,

```math
q_{\rho,\sigma}=T Df(q)q_\rho+T f_\rho+T_\rho f(q),
```

and the transverse matrix sensitivity obeys

```math
M_{\rho,\sigma}=T A M_\rho+(T_\rho A+T A_\rho)M.
```

The entries of $A_\rho$ are explicit using $C_\rho=-Z_\rho$, $c_{A,\rho}=C_\rho+bY_\rho$, $s_{A,\rho}=bC_\rho-Y_\rho$, and $V_{A,\rho}=-h_\rho$, including the direct occurrences of $\rho$ in A. This avoids a generic high-order derivative engine, although the resulting enlarged ODE still requires a validated enclosure. Differentiate the determinant polynomial directly near $g=0$; a formula using $(I+M_A)^{-1}$ is singular at the point being proved.

Even a certified simple critical root in one block becomes **sixfold** in the full system. The generic simple period-doubling theorem is therefore not directly applicable. Establishing doubled branches requires symmetry-specific nonlinear analysis. Reversible-system methods such as [4] provide useful methodological context, but their hypotheses do not automatically hold for this dissipative commuting symmetry.

## What must accompany a completed proof

The certificate should identify the exact model constants, parameter boxes, initial-state and half-period boxes, rounding backend, validated integration remainder bounds, Newton inclusion margins, uniqueness/contraction bounds, and spectral sign margins. All unsuccessful boxes and numerical limits should remain distinguishable from successful certificates. A small independent verifier should check the stored inequalities whenever feasible.

Proof of the local cycle and transverse crossing would still leave the observed [hyperchaotic candidate](motor-rows-activity-path.md) at numerical-evidence status. Proving chaos needs additional global information, such as suitable covering relations or a transverse homoclinic structure; proving an attracting hyperchaotic set requires still stronger trapping and expansion control. None follows solely from the six unstable directions of this periodic orbit.

## Primary references

1. F. Alouges et al., “Some mathematical models for flagellar activation mechanisms,” *Mathematical Models and Methods in Applied Sciences* **35**, 2395–2424 (2025). [DOI](https://doi.org/10.1142/S0218202525500423); [complete source v2](https://arxiv.org/html/2409.03506v2). Published physical model; the exact normalized field and source conventions are recorded in this repository's source audit.
2. T. Kapela, M. Mrozek, D. Wilczak, P. Zgliczynski, “CAPD::DynSys: A flexible C++ toolbox for rigorous numerical analysis of dynamical systems,” *Communications in Nonlinear Science and Numerical Simulation* **101**, 105578 (2021). [DOI](https://doi.org/10.1016/j.cnsns.2020.105578); [complete primary manuscript](https://arxiv.org/pdf/2010.07097), especially Section 4.2 and Theorem 4. The paper explicitly validates periodic points using flow/derivative enclosures and interval Newton inclusion; its implementation distinguishes rigorous and ordinary numerical solvers.
3. T. Kapela, D. Wilczak, P. Zgliczynski, “Recent advances in a rigorous computation of Poincare maps,” *Communications in Nonlinear Science and Numerical Simulation* **110**, 106366 (2022). [DOI](https://doi.org/10.1016/j.cnsns.2022.106366); [complete primary manuscript](https://arxiv.org/pdf/2104.08046), especially Sections 2–3. It treats return-time enclosures, section coordinates and the neutral flow multiplier.
4. I. Walawska, D. Wilczak, “Validated numerics for period-tupling and touch-and-go bifurcations of symmetric periodic orbits in reversible systems,” *Communications in Nonlinear Science and Numerical Simulation* **74**, 30–54 (2019). [DOI](https://doi.org/10.1016/j.cnsns.2019.03.005); [complete primary manuscript](https://arxiv.org/pdf/1903.01195). Its reversible/Hamiltonian setting is not assumed here; it illustrates why numerical detection and a rigorous nonlinear bifurcation proof are separate tasks.

The reductions in this note are direct algebraic deductions from the stated equations. No mathematical-priority claim is made for them.
