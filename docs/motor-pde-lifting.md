# Lifting motor-cycle results to the density PDE

Under the source's sinusoidal-rate assumptions, the closed first-harmonic ODE is an exact invariant part of the original motor-density PDE. Every omitted density component contracts by exactly $e^{-\tau}$ in a periodic Sobolev norm. On a periodic orbit its linearized return operator is exactly $e^{-T}I$. Consequently the PDE adds infinitely many stable directions, with no additional unstable or critical multipliers.

The argument below transfers the repository's finite-dimensional periodic-orbit results to the specified density PDE. It is an elementary deduction from the published equations, not a claim of mathematical priority. It proves neither a chaotic attractor nor a nonlinear period-doubling bifurcation.

## Assumptions and normalization

Fix chemical activity, $a_0>0$, positive $k,\eta,\ell$, and the source model with periodic spatial coordinate $\xi\in\mathbb R/\ell\mathbb Z$, potential difference $\Delta W=U\cos(2\pi\xi/\ell)$, constant total rate $\omega_A+\omega_B=a_0$, and

```math
\omega_B(\xi)=a_0/2+a_1\cos(2\pi\xi/\ell)+b_1\sin(2\pi\xi/\ell).
```

These are the numerical assumptions in Alouges et al., Section 3.1, printed PDF page 21. The broader theory permits additional odd forcing harmonics; the zero-tail statements here use their absence. The $N$-row density equations are Eq. (10) in the PDF, Eq. (15) in the HTML. The paper already discusses higher-mode relaxation in Proposition 2.2. [Primary PDF](https://arxiv.org/pdf/2409.03506v2), [primary HTML](https://arxiv.org/html/2409.03506v2).

Use angular coordinate $\theta=2\pi\xi/\ell$, time $\tau=a_0t$, and dimensionless density $q_j(\theta,\tau)=\ell Q_j(\xi,t)$. Set

```math
X_j=2\pi\Delta_j/\ell,\quad V_j=X_j',\quad
A=a_1/a_0,\quad B=b_1/a_0,\quad
\rho=k/(\eta a_0),\quad \gamma=2\pi^2U/(\eta\ell^2a_0).
```

The source equations become

```math
\partial_\tau q_j+V_j\partial_\theta q_j=-q_j+f(\theta),
\qquad f(\theta)=\frac12+A\cos\theta+B\sin\theta.
```

Define the mean and real Fourier coefficients by integrals, without a half-factor convention for the mean:

```math
m_j=\frac1{2\pi}\int_0^{2\pi}q_j\,d\theta,\qquad
C_{j,n}=\frac1\pi\int_0^{2\pi}q_j\cos(n\theta)\,d\theta,\qquad
S_{j,n}=\frac1\pi\int_0^{2\pi}q_j\sin(n\theta)\,d\theta.
```

Thus $q_j=m_j+\sum_{n\geq1}(C_{j,n}\cos n\theta+S_{j,n}\sin n\theta)$. The motor force is $-\pi U S_{j,1}/\ell$. Solving the force-difference equations with $\sum_jX_j=0$ gives

```math
V_j=-\rho X_j-\gamma(S_{j,1}-\overline S_1),
\qquad \overline S_1=\frac1N\sum_jS_{j,1}.
```

For the source's $A<0$, let $c_j=C_{j,1}/A$, $s_j=S_{j,1}/A$, $r=-\gamma A$, $b=B/A$. This gives exactly the [implemented first-harmonic field](motor-rows-source-audit.md), including the signs of its transport terms.

**Mean convention check.** On printed PDF page 12, Eq. (24) writes the constant term as $p_0/2$, while Eq. (25) evolves $p_0$ toward $1/(2\ell)$. Those two uses differ by a factor of two. We avoid that notation mismatch by integrating the PDE directly: the actual mean of $Q_j$ tends to $1/(2\ell)$, so $m_j$ tends to $1/2$. This concerns the displayed coefficient convention, not a change to the model's PDE or mechanical force.

## Exact mode equations and absence of feedback

Periodic integration by parts gives

```math
m_j'=\frac12-m_j,
```

```math
C_{j,n}'=-C_{j,n}-nV_jS_{j,n}+A\,\mathbf1_{n=1},\qquad
S_{j,n}'=-S_{j,n}+nV_jC_{j,n}+B\,\mathbf1_{n=1}.
```

In the convention $z_{j,n}=C_{j,n}+iS_{j,n}$,

```math
z_{j,n}'=(-1+inV_j)z_{j,n}+(A+iB)\mathbf1_{n=1}.
```

The usual complex Fourier coefficient with basis $e^{in\theta}$ is $(C_{j,n}-iS_{j,n})/2$ and therefore has the opposite phase sign. The two conventions describe the same translation.

For $n\geq2$, integration is exact:

```math
z_{j,n}(\tau)=e^{-\tau+in[X_j(\tau)-X_j(0)]}z_{j,n}(0).
```

Neither $m_j$ nor any $n\geq2$ coefficient enters $V_j$. The transport speed is independent of $\theta$, multiplication by $a_0$ mixes no Fourier modes, and the potential force selects only the first sine coefficient. These three facts give the exact closure; arbitrary mean or higher-mode initial data do not invalidate it.

## Function spaces and an exact nonlinear contraction

Take real periodic $H^s(\mathbb T)$ densities with $s\geq0$, including $L^2$ when $s=0$, and a real displacement vector with zero sum. Use the Fourier definition of the $H^s$ norm, so translations are isometries. Define the first-harmonic lift of the finite state $y=(C_{j,1},S_{j,1},X_j)_j$ by

```math
\mathcal L_j(y)(\theta)=\frac12+C_{j,1}\cos\theta+S_{j,1}\sin\theta,
\qquad R_j=q_j-\mathcal L_j(y).
```

The residual includes the non-equilibrium mean and all higher harmonics; its first harmonic is zero. Because $\mathcal L(y(\tau))$ solves the same forced PDE over the exact finite trajectory, subtraction gives

```math
\partial_\tau R_j+V_j\partial_\theta R_j=-R_j,
\qquad
R_j(\theta,\tau)=e^{-\tau}R_j(\theta-X_j(\tau)+X_j(0),0).
```

Hence, for the product norm over rows,

```math
\|R(\tau)\|_{(H^s)^N}=e^{-\tau}\|R(0)\|_{(H^s)^N}.
```

The same equality holds in translation-invariant $L^p$ norms, $1\leq p\leq\infty$, and standard periodic $C^k$ norms when defined. Two solutions with the same initial finite projection have identical mechanical trajectories, and their complete density difference obeys this equality.

This also constructs a unique global mild solution for the stated $H^s$ data. The finite first-harmonic amplitudes satisfy

```math
|z_{j,1}(\tau)|\leq e^{-\tau}|z_{j,1}(0)|+(1-e^{-\tau})\sqrt{A^2+B^2}.
```

The mechanical equation is a stable linear filter driven by these bounded amplitudes, so its finite-dimensional solution has no forward finite-time escape. The explicit translated residual then gives $q\in C([0,\infty);(H^s)^N)$, solving the PDE in the mild or distributional sense. For periodic $C^1$ initial densities it gives the classical transport solution. This is a direct construction under the assumptions above, rather than an invocation of a general semilinear-PDE theorem.

The affine manifold $R=0$ has dimension $3N-1$, is invariant, and attracts the residual at rate one. That statement does not imply attraction to any particular periodic orbit. A projected low-harmonic shadow of an arbitrary physical initial density need not itself be pointwise nonnegative during its initial transient.

## Periodic solutions and their complete linearized spectrum

**Periodic lifting theorem.** Every finite first-harmonic periodic orbit lifts to a periodic density solution by $q=\mathcal L(y)$. Conversely, any full-PDE periodic solution in the stated spaces has $R=0$: if its period is $T>0$, periodicity and the norm identity imply $\|R(0)\|=e^{-T}\|R(0)\|$. The displacement and first harmonics can be recovered from the full state, so the correspondence preserves the least positive period of every nonconstant orbit.

Linearize about a lifted periodic solution. Its derivative with respect to $\theta$ has only first harmonics. Therefore the term involving a speed perturbation, $\delta V_j\partial_\theta q_j^*$, contributes only to the finite first-harmonic block. The residual linearization is exactly the homogeneous transport equation above. Since $X_j^*(T)=X_j^*(0)$ as a real displacement,

```math
M_{\rm PDE}=M_{\rm finite}\oplus e^{-T}I_{\rm residual},\qquad
\sigma(M_{\rm PDE})=\sigma(M_{\rm finite})\cup\{e^{-T}\}.
```

This is the bounded monodromy operator of the linearized PDE on the complexification of $(H^s)^N$ with constrained displacements. The residual subspace includes $N$ mean directions and infinitely many higher-mode directions. Its multiplier $e^{-T}$ has infinite multiplicity, and all its normalized Floquet growth exponents are $-1$. In dimensional time the decay rate is $a_0$ and the multiplier is $e^{-a_0T_{\rm physical}}$.

In particular, the [certified eight-row target orbit](motor-cycle-certificate.md) has exactly six unstable directions and one neutral flow direction in the density PDE; its sixteen stable finite-part directions are supplemented by infinitely many stable density directions. At the [certified transverse crossing](motor-crossing-certificate.md), exactly six critical multipliers equal $-1$, the flow multiplier equals $1$, and all other directions are stable. The omitted modes cannot change that multiplicity or the certified crossing speed.

**Regularity qualification.** The spectrum above is the spectrum of the linearized transport evolution. Translations are strongly continuous on $H^s$ but not continuous in operator norm; dependence on a translation parameter can cost one spatial derivative. We do not assume that the nonlinear Eulerian solution map is Fréchet differentiable on the same $H^s$ space. One can instead use the continuous invertible moving-coordinate change $W_j(\theta)=R_j(\theta+X_j)$, for which $W_j'=-W_j$ exactly. Nonlinear stability below follows directly from norms, without a smooth-semigroup or normal-hyperbolicity assertion.

## Nonlinear orbital stability and physical densities

Let $\Gamma$ be a periodic finite-dimensional orbit, measured with the norm induced by its density lift and the displacement norm. Fourier orthogonality makes the following identity exact:

```math
\operatorname{dist}_{s}\big((q(\tau),X(\tau)),\mathcal L\Gamma\big)^2
=\operatorname{dist}_{s}\big(y(\tau),\Gamma\big)^2
+e^{-2\tau}\|R(0)\|_{(H^s)^N}^2.
```

The notation $\mathcal L\Gamma$ includes the displacement component. Thus finite-dimensional orbital stability, asymptotic orbital stability, or instability transfers to the density PDE. If the finite cycle has a local exponential orbital-attraction estimate with rate $\alpha>0$, the full estimate has rate $\min(\alpha,1)$. The subsequently [certified attracting four-pair orbit](motor-paired-certificate.md) meets this full-system stability condition and therefore lifts to an attracting PDE orbit. Stability only inside a cluster restriction would be insufficient.

The [global-attractor extension](motor-global-attractor.md) combines this exact splitting with explicit absorbing bounds. It proves that the entire compact global PDE attractor is the first-harmonic lift of the finite global attractor, and that convergence basins correspond under first-harmonic projection. This concerns global dissipativity, not a chaos or scientific-priority assertion.

Nonnegative transition rates require $0\leq f\leq1$, equivalently $\sqrt{A^2+B^2}\leq1/2$. Along characteristics the full density obeys a scalar relaxation equation, so $0\leq q_j(\theta,0)\leq1$ implies $0\leq q_j(\theta,\tau)\leq1$ for all forward times (almost everywhere for rough data). These are the probabilities of state A; the complementary state has density $(1-q_j)/\ell$.

Every finite periodic orbit automatically gives a physical lift under this rate condition. Indeed,

```math
\frac12(|z_{j,1}|^2)'=-|z_{j,1}|^2+AC_{j,1}+BS_{j,1}.
```

At a maximum over a period this implies $|z_{j,1}|\leq\sqrt{A^2+B^2}$, hence

```math
\frac12-\sqrt{A^2+B^2}\leq q_j^*(\theta,\tau)
\leq\frac12+\sqrt{A^2+B^2}.
```

For a strict rate inequality the periodic density is uniformly separated from both physical boundaries. In $H^s$ with $s>1/2$, sufficiently small density perturbations are also pointwise admissible by Sobolev embedding. In $L^2$ the stability statement can instead be restricted to admissible initial densities; no interior of the positive cone is assumed.

The conclusions depend on constant total rate, first-harmonic potential, fixed parameters and the stated spatially uniform transport speeds. Spatially varying total rates or additional force harmonics can restore feedback from omitted modes. With higher forcing harmonics, the tails are driven rather than zero; the simple zero-tail invariant manifold used here no longer applies. A bending filament coupled to fluid dynamics adds other states and is outside this theorem.

The [source and algebra record](../results/motor-pde-lifting/audit.json) identifies the inspected PDF and the normalization checks. The [exact polynomial controls](../results/motor-pde-lifting/check_lifting_algebra.py) verify 13 sign, scaling and amplitude identities; their [output](../results/motor-pde-lifting/algebra.json) records the executed-source hash, with the original script preserved alongside it. Run them with `python results/motor-pde-lifting/check_lifting_algebra.py --output outputs/pde-lifting-algebra.json` in the proof environment. Choose a new output path. These controls supplement the function-space argument; no new trajectory integration or infinite-dimensional numerical truncation is needed for the lifting theorem.
