# Exact structure of scalar-quadratic Adamax

These elementary deductions concern a known optimizer, not a newly proposed map. They explain the folded sheets in the [numerical geometry study](adamax-geometry.md). The statements below establish identities and boundedness under specified assumptions; they do not establish chaotic attraction or historical priority.

Write the autonomous normalized update as

$$
M=bm+(1-b)x,\qquad U=\max(cu,|x|+e),\qquad X=x-M/U,
$$

with state `(x,m,u)`, output `(X,M,U)`, `0<b,c<1`, and `e>0`. The half-space `u>=0` maps into `U>=e`, so division is well-defined. For the actual cold-start optimizer, the last equation is `X=x-M/[(1-b^t)U]` at update `t`. The autonomous map is its large-time limit.

## Reset surfaces and exact loss of one input direction

On the reset branch `|x|+e>=cu`, set `sigma=sign(x)`. Eliminating the old position gives

$$
\boxed{M=\sigma U(U-e)-XU.}
$$

At `x=0`, either sign gives the same relation because `U=e`. Thus every reset image lies on one of two quadratic surfaces, each a graph over `(X,U)`. This is an exact identity for the autonomous map. For the time-dependent bias-corrected update, the right-hand side is multiplied by `1-b^t`; its surfaces approach the displayed ones.

Within a smooth reset branch, the derivative's third column is exactly zero. Its rank is two for `b>0` and `x!=0`. A derivative product containing such a reset has a zero singular value, so its third product exponent is minus infinity. A large finite negative value from three-column floating-point QR would obscure this exact rank loss. Our numerical checks therefore evolve the two surviving tangent directions. Derivatives at the switching boundary and at the absolute-value kink require separate treatment; they are not certified by this branch calculation.

On the memory branch, the map is

$$
G(x,m,u)=\left(x-\frac{bm+(1-b)x}{cu},\;bm+(1-b)x,\;cu\right),
\qquad \det DG=bc.
$$

The determinant follows by first replacing `(m,u)` with `(M,U)` and then applying the shear `X=x-M/U`. Subsequent memory steps transport a reset surface through `G`, producing the observed layers. On the two-million-state preview, the largest observed number of consecutive memory steps is four for case 847 and one for case 4099. This finite observation is not a global return-time bound.

The symmetry `R(x,m,u)=(-x,-m,u)` satisfies `F R=R F`. The two sign-related halves are not separate discoveries.

## Every bounded autonomous limit set lies on finitely many sheets

**Statement.** For `e>0` and `0<c<1`, the omega-limit set of any bounded forward orbit of this autonomous map has Hausdorff dimension at most two.

**Proof.** Choose a finite bound `Umax` for its memory coordinate and an integer `L` satisfying `c^L Umax<e`. There cannot be `L` consecutive memory-branch updates: they would produce `u<e`, whereas every update has `u>=e`. After the initial `L` steps, every state is therefore at most `L-1` memory updates past a reset.

Choose a compact box containing the orbit and its reset outputs. Let `S+` and `S-` be the portions of the two reset surfaces in that box. Every sufficiently late state belongs to

$$
\bigcup_{\sigma\in\{+,-\}}\;\bigcup_{k=0}^{L-1}G^k(S_\sigma).
$$

The extension of `G^k` to each compact patch is smooth: its denominator stays positive, at least `c^k e`. Each image is compact and is a Lipschitz image of a two-dimensional surface patch. The finite union is closed and has Hausdorff dimension at most two; it contains the omega-limit set. This proves the statement.

The statement is conditional on orbit boundedness. It does not prove that an arbitrary initial state is bounded, that the limit set is chaotic, or that its dimension equals two. It also does not rule out a topological conjugacy with familiar attractors: Hausdorff dimension is not invariant under arbitrary homeomorphisms. When `e=0`, the finite-memory-age argument no longer applies.

## A bound for cold starts when the first memory decays faster

**Statement.** Suppose `0<b<c<1`, `m0=u0=0`, and `e>0`. For both the autonomous update and the continuously bias-corrected Adamax update,

$$
|x_t|\le B,\quad |m_t|\le B,\quad e\le u_t\le B+e\quad(t\ge1),
\qquad B=\max\left(|x_0|,\frac{c}{c-b}\right).
$$

**Proof.** The maximum-memory recurrence gives `u_t>=c^k |x_(t-1-k)|`. Consequently, for `mhat_t=m_t/(1-b^t)`,

$$
\frac{|\widehat m_t|}{u_t}
\le\frac{1-b}{1-b^t}\sum_{k=0}^{t-1}(b/c)^k
\le K:=\frac{1-b}{1-b/c}.
$$

The last inequality uses `(b/c)^t>=b^t`. The bias-corrected step has magnitude at most `K`; the autonomous step is no larger. Hence positions `k` updates apart differ by at most `kK`.

The normalized first moment is a truncated geometric average of past positions. Its mean lag is at most `b/(1-b)`, giving

$$
|\widehat m_t-x_{t-1}|\le D:=\frac{bK}{1-b}.
$$

When `|x_(t-1)|>D`, the update points toward zero. Its magnitude is at most `K`, so an overshoot has magnitude at most `K`. Inside `[-D,D]`, the next position has magnitude at most `D+K=c/(c-b)`. Induction proves the position bound. The first-moment and maximum-memory bounds follow directly from their recurrences.

This proves boundedness of the selected cold starts for 847, 239 and the standard default pair 4096. It does not cover case 4099, whose `b>c`, or arbitrary nonzero initial memories. It is not a claim that the rectangular box itself is forward invariant for every point in it. Combined with the preceding statement, it bounds the Hausdorff dimension of the corresponding autonomous cold-start limit sets by two.

## A closed period-two branch

Let `kappa=(1-b)/(1+b)` and `a=kappa/2-e>0`. Direct substitution gives the orbit

$$
(a,-\kappa a,\kappa/2)\ \longleftrightarrow\ (-a,\kappa a,\kappa/2).
$$

Both points strictly select the reset branch whenever `c<1`. The nonzero two-step multipliers are the eigenvalues of `B2^2`, where

$$
B_2=\begin{pmatrix}
1-2b-4e/\kappa&-2b/\kappa\\
1-b&b
\end{pmatrix},
\qquad
\operatorname{tr}B_2=1-b-4e/\kappa,\quad
\det B_2=b(3-4e/\kappa).
$$

The third multiplier is zero. Local stability of this branch is independent of `c`. For `e=5e-6`, its lower unit-determinant crossing is the smaller root of

$$
(3+4e)b^2-(4-4e)b+1=0.
$$

This identifies a linear unit-modulus crossing. We have not computed the nonlinear coefficient needed for a Neimark–Sacker bifurcation theorem, nor established a route from this branch to the displayed chaotic candidates.

The [structural controls](../results/adamax-geometry/structure/analysis.json) check 256 period-two identities and full two-step derivatives, with maximum absolute discrepancy below `8.6e-14`. The reset-sheet residuals are also checked directly on actual states. These numerical controls support the implementation; the algebra above supplies the arguments. A separately closed attracting four-cycle at case 239 illustrates coexistence and is retained in the same packet.
