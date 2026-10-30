# Follow-up priority audit and an even-row corollary

Audit date: 30 October 2026. This supplements the [earlier priority audit](motor-rows-priority-audit.md) and the [periodic-orbit certificate](motor-cycle-certificate.md).

**Priority remains unresolved.** This follow-up strengthens the comparison with a 2026 citing paper and identifies a second relevant author thesis that has not been inspected in full. It also derives a consequence of the existing certificate: the same periodic saddle exists in every even-row mathematical system, with exactly $N-2$ unstable Floquet directions at the certified target parameter. This is an algebraic extension of the repository's result, without a claim of historical novelty or of a new chaotic attractor.

## What the additional source checks establish

**The exact source.** Alouges et al., arXiv v2, Section 3.3 and Fig. 9, use eight motor rows to represent the axonemal configuration with nine doublets and a fixed bridge. Their examples display periodic two-group motion; general cluster selection and the higher-dimensional bifurcation problem are left open. The general equations precede those examples. Thus the source itself supports the mathematical $N$-row family, while its biological eight-row interpretation is specific. [Complete primary source](https://arxiv.org/html/2409.03506v2).

**Agostinelli et al. (2026).** The official repository's indexed PDF excerpt for Section 4.1, printed page 11, explicitly replaces the microscopic motor kinetics by a cubic active-force law. Its Eq. (35) is

```math
\dot f_a=-\frac1\tau\left(f_a-\lambda\dot u+\mu\tau\dot u^3\right).
```

The projected cilium equation, Eq. (40), is a second-order oscillator with cubic velocity damping. These are approximations to mechanochemical dynamics, not an identification with the exact $3N-1$ first-harmonic field. This is now a coefficient-level distinction rather than an inference from an abstract. The directly downloaded authoritative PDF remains inaccessible (HTTP 403). A public transcription of the authors' pre-proof also exposes Sections 4.1–4.2 and the discussion; it agrees on this model choice, but is supplementary evidence and was not verified against every page of the published PDF. Its machine-extracted signs and formatting are not used where ambiguous. [Official record](https://www.iris.sssup.it/handle/11382/587254), [indexed primary PDF](https://www.iris.sssup.it/retrieve/1bc1e5d7-7032-41ad-ade0-cd5e0efeb690/1-s2.0-S0022509626000293-main-4.pdf), [supplementary pre-proof transcription](https://www.scribd.com/document/990072681/1-s2-0-S0022509626000293-main).

**Anello's 2025 thesis remains a gap.** The institutional record and indexed abstract are available, but the complete thesis could not be read. An abstract describing Hopf results and periodic examples cannot exclude additional results elsewhere in a thesis. [Author's institutional page](https://www.math.sissa.it/users/irene-anello), [thesis record](https://iris.sissa.it/handle/20.500.11767/147991).

**Levillain's 2024 thesis is another relevant gap.** The author's website and official HAL metadata identify *Filament models for swimming at the microscopic scale*, defended 27 September 2024 and deposited 5 November 2024. Its abstract explicitly places the two-row analysis and the all-row numerical model in its final chapter. Ordinary document requests returned an anti-automation HTML page, not a PDF; the full chapter was not inspected. This is a particularly relevant missing source because Levillain coauthored the exact model. [Author's page](https://jesnial.github.io/), [official thesis record](https://theses.hal.science/tel-04766775v1), [public HAL metadata query](https://api.archives-ouvertes.fr/search/?q=halId_s%3Atel-04766775&fl=title_s,uri_s,fileMain_s,files_s,abstract_s,docType_s,nntId_s,doiId_s,defenseDate_s,submittedDate_s&wt=json).

**Cluster breakup is established prior art.** Haugland and Krischer, *Physical Review E* **103**, L060201 (2021), analyze four globally coupled Stuart–Landau oscillators. A half-period-exchanging $2+2$ state loses stability by equivariant period doubling; $2+1+1$ and fully asymmetric branches undergo cascades to chaos, connected by symmetry-breaking bifurcations. Their Eq. (1) uses complex cubic oscillators and nonlinear mean coupling, with harmonic mean motion; it is not the motor displacement/Fourier field. This does not identify our exact equations with theirs, but prevents presenting cluster-splitting period doubling and subsequent chaos as a newly invented general mechanism. A matching local symmetry pattern would support comparison, not establish global equivalence or priority. [Complete primary preprint, Eq. (1), Fig. 1 and associated discussion](https://arxiv.org/html/2102.10138v1), [journal DOI](https://doi.org/10.1103/PhysRevE.103.L060201).

The earlier pair-model equivalence to Oriola et al. and related sarcomere-chaos precedents remain in force. These checks do not convert a failure to find an exact-match paper into evidence that none exists. The [retrieval record](../results/motor-priority-followup/source-audit.json) separates primary full text, indexed excerpts, metadata and unread documents.

## Exact replication to every even row count

This part is a deduction from the normalized equations and the existing certificate, separate from the literature assessment. Fix

```math
b=1,\qquad r=1+\frac{19/20}{\sqrt{32001/500000}},\qquad
\rho=\frac{1291}{440}.
```

For every integer $m\geq1$, take $N=2m$ and partition the rows into two groups of size $m$. Let $(h,Y,C)$ be the certified nonconstant periodic core orbit:

```math
h'=(-\rho+rC)h-rY,\qquad Y'=hC-Y,\qquad C'=1-C-hY.
```

Assign every A row and every B row the respective states

```math
(c_A,s_A,X_A)=(C+Y,C-Y,(h-rY)/\rho),
```

```math
(c_B,s_B,X_B)=(C-Y,C+Y,-X_A).
```

Equal group sizes give $\bar s=C$ and $\sum_jX_j=0$, independently of $m$. Substitution gives $V_A=-h$, $V_B=h$ and exactly the three core equations. Thus the orbit and its period are unchanged by replication. Conversely, $C=(c_A+c_B)/2$, $Y=(c_A-c_B)/2$ and $h=\rho X_A+rY$ recover the core from any embedded state. Therefore the embedding preserves the least positive period certified for that core.

The tangent space splits as

```math
\mathcal T=\mathcal T_{\rm cluster}\oplus
(\mathcal U_A\otimes\mathbb R^3)\oplus
(\mathcal U_B\otimes\mathbb R^3),
\qquad \dim\mathcal U_A=\dim\mathcal U_B=m-1.
```

Here $\mathcal U_A$ consists of coefficients summing to zero within A, and similarly for B. The cluster component has dimension five because two identical row states have six coordinates and one mechanical constraint. Hence $5+6(m-1)=3N-1$.

Each contrast has zero change in $\bar s$. Along its own row state it obeys

```math
J_\perp=
\begin{pmatrix}
-1&-V-rs&\rho s\\
V&-1+rc&-\rho c\\
0&r&-\rho
\end{pmatrix}.
```

There are $m-1$ copies of this block in each group. The certified half-period exchange makes the A and B monodromies similar, so every block multiplier occurs $2(m-1)=N-2$ times. At the target the three simple roots are isolated in

```math
(-3.71,-3.70),\qquad(-0.012,-0.010),\qquad(0.420,0.423).
```

The cluster block retains one neutral and four stable directions. Two stable directions are the core's non-flow modes; the two additional cluster deviations obey $d'=-d+he$, $e'=-e-hd$, so their norm decays exactly as $e^{-t}$. Consequently the full target orbit has

```math
\boxed{N-2\ \text{unstable},\quad 1\ \text{neutral},\quad 2N\ \text{stable directions}.}
```

For $N=2$ there are no contrasts, so this same orbit is orbitally attracting in its five-state system. For $N\geq4$ it is a saddle. These counts concern a periodic orbit; they do not count the positive exponents of a chaotic trajectory, establish a chaotic attractor for any $N$, or explain the odd-row $N=3$ candidate.

The separate certified family over $I=[3.370218399569,3.370218399589]$ also replicates exactly. The inherited endpoint counts become $(N-2,1,2N)$ at the lower endpoint and $(0,1,3N-2)$ at the upper endpoint. For every even $N\geq4$, an interior parameter has multiplier $-1$ with multiplicity at least $N-2$, directly from the inherited continuity argument. The subsequently completed [crossing certificate](motor-crossing-certificate.md) sharpens this: in its smaller interval the crossing is unique, has nonzero speed and multiplicity **exactly $N-2$**. This requires neither a nonlinear flip theorem nor a claim that the long parameter interval joining this family to the target has been validated.

## General balanced-cluster replication identity

The same calculation applies to any $k$-row solution, with $k\geq2$, replicated $m$ times per row into $N=mk$ rows. Write the row states as $q_a=(c_a,s_a,X_a)$, $a=1,\ldots,k$, with $\sum_aX_a=0$. The replicated mean is

```math
\frac1{mk}\sum_{a=1}^k\sum_{i=1}^m s_a
=\frac1k\sum_{a=1}^k s_a.
```

Consequently the balanced $k$-cluster quotient is exactly the original $3k-1$ dimensional $k$-row field, independent of $m$. Along a periodic quotient orbit, let $M_k$ be its constrained quotient monodromy and $M_a$ the three-dimensional zero-mean contrast monodromy for cluster $a$. The invariant tangent splitting gives the characteristic-polynomial identity

```math
\det(\zeta I-M_N)=\det(\zeta I-M_k)
\prod_{a=1}^{k}\det(\zeta I-M_a)^{m-1}.
```

The dimensions match: $(3k-1)+3k(m-1)=3N-1$. This factorization holds for the specified balanced clustering without assuming that the periodic orbit has a rotating symmetry.

If a time shift cyclically permutes all $k$ cluster states, their variational coefficients are time shifts of each other. The monodromies at different time origins are similar, so all $M_a$ have the same characteristic polynomial. The factorization reduces to

```math
\det(\zeta I-M_N)=\det(\zeta I-M_k)
\det(\zeta I-M_\perp)^{N-k}.
```

Thus each simple transverse root contributes $N-k$ copies; coincidences with quotient roots would add to its full-system multiplicity. If the quotient has one flow multiplier and all other multipliers strictly inside the unit circle, and a transverse block has exactly one simple root outside, the full cycle has $N-k$ unstable, one neutral and $2N+k-2$ stable directions.

This explains the dimension counts relevant to the numerical follow-up: for $N=8$, a $4+4$ cycle has $k=2$ and six possible expanding contrast copies; a rotating $2+2+2+2$ cycle has $k=4$ and four. The latter count applies only when the stated quotient stability, rotating symmetry and transverse spectrum actually hold. Existence and bifurcation properties of that four-cluster cycle require their own numerical or interval evidence. The identity alone proves neither the transition between the two cycles nor a cascade to chaos. A symmetry-forced multiple root at $-1$ still does not satisfy the hypothesis of a generic simple flip theorem.

For an **unequal partition**, with $m_a\geq1$ rows in group $a$ and $N=\sum_a m_a$, the quotient instead has weighted mean $\bar s=\sum_a(m_a/N)s_a$ and constraint $\sum_a m_a X_a=0$. These weights cannot be discarded. It is still a $3k-1$ dimensional invariant system, and its periodic monodromy satisfies

```math
\det(\zeta I-M_N)=\det(\zeta I-M_{\rm weighted})
\prod_{a=1}^{k}\det(\zeta I-M_a)^{m_a-1}.
```

To see this, split each cluster's perturbation into its row average and a zero-sum remainder. A remainder leaves the weighted mean and mechanical constraint unchanged; all its independent coefficients obey that cluster's $J_\perp$. The average component evolves by the weighted quotient. The dimension identity is $(3k-1)+3\sum_a(m_a-1)=3N-1$. No equality of different cluster spectra follows without an additional verified symmetry. For example, the eight-row partition $2+2+1+1+1+1$ has a 17-state quotient and two three-state contrast blocks. Quotient attraction alone cannot establish attraction in the full 23-state system. This algebra permits partial cluster splitting and does not predict an all-rows splitting cascade.

## Physical closure and scope

No $N$-dependent rescaling of the row's chemical rates or mechanical constants is needed for this mathematical replication. The coupling already averages the $N$ sine coefficients. For the periodic core, $(C^2+Y^2)'=2C-2(C^2+Y^2)$ implies $C^2+Y^2\leq1$ by considering its maximum over a period. Every replicated row therefore satisfies $c_j^2+s_j^2=2(C^2+Y^2)\leq2$.

The source conditions $|A|\sqrt2\leq1/2$ for nonnegative transition rates and $\ell Q_j\in[1/2-|A|\sqrt{c_j^2+s_j^2},\,1/2+|A|\sqrt{c_j^2+s_j^2}]$ for the reconstructed first-harmonic density are unchanged with $N$. Mechanical endpoint closure is also unchanged: set $\Delta_j=\ell X_j/(2\pi)$ and reconstruct $x_i=\sum_{j\leq i}\Delta_j$; then $x_0=x_N=0$. This does not preserve every absolute filament displacement as $N$ grows, nor does it assert an axoneme with arbitrary row count. It transfers the theorem within the published mathematical family.

The [exact algebra check](../results/motor-priority-followup/check_even_rows.py) verifies 22 polynomial identities over rational coefficients, including the embedding, both transverse blocks and the physical-radius relations. Its [output](../results/motor-priority-followup/even-rows-algebra.json) records the inherited certificate hash. It does not repeat the interval integrations or replace the direct-sum argument above. Reproduce with:

```bash
python results/motor-priority-followup/check_even_rows.py --output outputs/even-rows-algebra.json
```

The useful claim at present is an exact-model periodic-orbit stability theorem with a dimension-dependent consequence. Establishing priority, a nonlinear route to chaos and a chaotic attracting invariant set remain distinct tasks.
