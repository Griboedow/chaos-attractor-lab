# What three sinks force in the driven cubic

Let $`P`$ be the time-$`2\pi`$ map at forcing phase zero. At gain 0.34 the [moving-disk certificate](../shape_search_20261121/frontiers/cubic/finite-gains.json) gives three sinks with chronological braid $`(\sigma_1^{-1}\sigma_2)^3`$. We use the [original gradient equation](README.md#system) throughout.

## Entropy in the original flow

The radial component points strictly inward for $`r=|z|\ge3`$, since $`r^6-4r^4+3r^2-(29/4)(r^3+r)>0`$ there. Extend to the radius-four disk by multiplying the field by a smooth radial function equal to one up to radius three, positive inside radius four and zero on its boundary. The resulting flow is complete on the larger disk. Its time map fixes the boundary. Every point of the open outer annulus is wandering by strict radial decrease. The fixed boundary carries zero entropy. Thus the entropy is supported on the original disk, by the nonwandering-set characterization of topological entropy.

In slope coordinates the elementary twists and the one-period braid act by

```math
S_1=\begin{pmatrix}1&1\\0&1\end{pmatrix},\quad
S_2=\begin{pmatrix}1&0\\-1&1\end{pmatrix},\quad
A=(S_2 S_1^{-1})^3=\begin{pmatrix}5&-8\\-8&13\end{pmatrix}.
```

These are the reduced Burau matrices at $`t=-1`$ in a homology basis of the double cover. The braid entropy bound for disk homeomorphisms gives

```math
h_{\rm top}(P)\ \ge\ \log\rho(A)=\log(9+4\sqrt5).
```

The unit is one forcing period. Recurrence places any ergodic positive-entropy measure outside the three sink basins up to a null set. Whether a chaotic set attracts an open set remains open. [Band and Boyland, *The Burau estimate for the entropy of a braid*](https://arxiv.org/abs/math/0612716) · [Periodic orbits as moving obstacles](https://arxiv.org/abs/nlin/0510075) · [Braids in driven oscillators](https://doi.org/10.1142/S0218127493001100)

## Counting periodic points

Collapse the outer boundary to a repelling fixed point. The resulting sphere has four marked fixed points, the three sinks and the repeller. Its double cover branched there is a torus with involution $`x\mapsto-x`$. The two lifts of $`P^n`$ act on torus homology by $`B=A^n`$ and $`B=-A^n`$. Write $`T_n=\mathrm{tr}(A^n)`$.

Choose the origin at a branch point. Each universal lift $`F:\mathbb R^2\to\mathbb R^2`$ fixing zero is odd and satisfies $`F(x+k)=F(x)+Bk`$ for integer $`k`$. For any integer vector $`v`$, the map $`x\mapsto F(x)-x-v`$ differs from the invertible linear map $`B-I`$ by a bounded function. Its degree is the sign of $`\det(B-I)`$. Each residue class of $`v`$ modulo $`(B-I)\mathbb Z^2`$ therefore supplies a fixed point of the torus map, giving at least $`|\det(B-I)|`$ points.

At the four half-lattice points $`s`$, oddness gives $`F(s)=Bs`$. The four classes of $`(B-I)s`$ are distinct. For $`B=A^n`$ the class degree is −1. Its branch point has local index +1 since it is a sink or a repeller also upstairs. There must be another zero in each such class. The involution pairs it with a distinct nonbranch zero in that same class since $`2(B-I)s\in(B-I)\mathbb Z^2`$. This adds eight torus points for this lift.

The two elementary bounds are $`T_n-2`$ and $`T_n+2`$. Each fixed sphere point contributes two to the combined count upstairs, including each branch point, counted once by each lift. Removing the point at infinity gives

```math
\#\mathrm{Fix}(P^n)\ \ge\ \tfrac12[(T_n-2)+8+(T_n+2)]-1=T_n+3.
```

The bounds for one, two and three returns are 21, 325 and 5,781 section points of period dividing $`n`$. [Exact matrix and counts](bounds.json) · [Boyland on braid and fixed-point methods](https://people.clas.ufl.edu/boyland/files/topmethods.pdf)

## Increasing entropy at a fixed fundamental period

For any positive integer $`m`$ use the same equation with

```math
b_m(\theta)=\tfrac14\sin\theta-4\sin(3m\theta)+3i\sin(6m\theta),\qquad \kappa\ge1200m.
```

The amplitude stays below 7.25 and the fundamental period stays $`2\pi`$. Gain and bandwidth increase with $`m`$. The three attracting cycles force at least $`m\log(9+4\sqrt5)`$ entropy per period.

For clearance from the critical values $`\pm2`$, put $`x=|\sin(3m\theta)|`$. If $`x<1/4`$ or $`x>3/4`$, the absolute real part is respectively below 1.25 or above 2.75. Otherwise the absolute imaginary part is $`6x\sqrt{1-x^2}\ge3\sqrt7/8>1/2`$. These bounds persist while removing the first harmonic. The remaining drive traverses the same root braid $`m`$ times, with spectral radius $`\rho(A)^m`$. Indeed its slit crossings alternate between the two outer rays with opposite twist signs, as in the [certified primitive drive](../shape_search_20261121/frontiers/cubic/#recovering-the-link-from-the-coefficients).

Use the three components of $`|p|\le1/32`$. They lie inside radius three. The [gradient tube estimate](../shape_search_20261121/frontiers/cubic/#rigorous-gain-regimes-for-the-gradient-flow) gives $`|p'|^2\ge9(1/2-1/32)/5=27/32`$. This exceeds $`6|z|/32`$, giving contraction. Also $`\kappa(27/32)/32>30m+1/4\ge\max|b'_m|`$, giving strict trapping. The unique cycles in these tubes have the root braid. The disk and entropy arguments apply for every $`m`$.
