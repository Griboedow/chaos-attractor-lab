# Structure of the existing smooth Yogi map

This note studies the pre-existing TensorFlow Addons `Yogi(activation="tanh")` option on a scalar quadratic. The algorithm and the occurrence of optimizer oscillations are prior work. The exact statements below describe this map; they do not establish a new attractor class or historical priority. The [experiment report](yogi-tanh.md) separates numerical evidence from these statements.

## Map and normalization

For `0<=b<1`, `0<c<1`, `e>0`, `k>0` and state `(x,m,v)` with `v>=0`, write

```math
T=\tanh(k(v-x^2)),\quad M=bm+(1-b)x,\quad
V=v-(1-c)Tx^2,\quad F(x,m,v)=(x-M/(\sqrt V+e),M,V).
```

In the [existing library implementation](https://github.com/tensorflow/addons/blob/d208d752e98c310280938efa939117bf635a60a8/tensorflow_addons/optimizers/yogi.py), the coefficient inside `tanh` is 10 in physical units. On `f(theta)=h*theta^2/2`, with `h>0`, the substitutions

```math
\theta=\alpha x,\quad m_{\rm phys}=\alpha h m,\quad
v_{\rm phys}=(\alpha h)^2v,\quad
e=\epsilon/(\alpha h),\quad k=10(\alpha h)^2
```

give the displayed map after the bias factors have decayed. Changing curvature can therefore compensate a change of learning rate. These are equivalent normalizations, not different discoveries. During the library transient, the first coordinate uses `r_t*M` instead of `M`, where `r_t=sqrt(1-c^t)/(1-b^t)`. Epsilon remains outside the uncorrected square root.

## All forward orbits are bounded

The legal half-space is invariant. If `v>=x^2`, then `V>=v-(1-c)x^2>=c*v`. If `v<x^2`, then `V>=v`. Moreover `V>0` whenever `(x,v)!=(0,0)`. Thus the map is continuous on the closed half-space and real analytic where `v>0`.

Set

```math
\sigma=\min\{e\sqrt{k},\sqrt c,1/\sqrt2,
\sqrt{(1-c)\tanh(1/2)}\}>0.
```

The denominator `d=sqrt(V)+e` satisfies `d>=e` and `d>=sigma*|x|`. To see the second inequality, for `|x|<1/sqrt(k)` use `d>=e`. For larger `|x|`, split the memory into three ranges. If `v>=x^2`, then `V>=c*x^2`. If `x^2/2<=v<x^2`, then `V>=v>=x^2/2`. If `v<x^2/2`, then `k(x^2-v)>=1/2`, so `V>=(1-c)*tanh(1/2)*x^2`.

The [scalar first-exit theorem](scalar-optimizer-boundedness.md) now bounds both `x` and `m` from every finite legal initial state, for every `0<=b<1`. If `|x|<=B`, then

$$v_t\le\max\{v_0,(2-c)B^2\}.$$

Indeed, above `B^2` the second memory cannot increase, and below it the next value is at most `(2-c)B^2`. The same argument covers the time-dependent bias factors: `0<r_t<=1/(1-b)`, so the effective denominator `d/r_t` has floors `(1-b)e` and `(1-b)sigma*|x|`.

Every autonomous orbit consequently has a nonempty compact invariant omega-limit set. Continuity and compactness prove invariance by taking convergent subsequences of both forward images and predecessors. This does not identify that set as chaotic or prove attraction of a neighborhood. The fixed line `(0,0,v)` is unbounded, so there is no compact global attractor that attracts all finite legal states.

## Folds and noninvertibility

For `v>0`, the derivative determinant is

```math
\det DF=b[1-(1-c)kx^2\operatorname{sech}^2(k(v-x^2))].
```

For `b>0`, its critical sheets satisfy

```math
v=x^2\ \pm\ \frac{\operatorname{arcosh}\sqrt{(1-c)kx^2}}k,
\qquad (1-c)kx^2\ge1,\quad v\ge0,
```

with arbitrary `m`. The map is globally noninjective for every legal parameter tuple. Choose `x` with `(1-c)kx^2>1`. The derivative of `v -> V` is negative at `v=x^2` and positive for sufficiently large `v`; hence two legal values of `v` have the same image. With `m=-(1-b)x/b`, both have `M=0` and `X=x`, so the full images agree. For `b=0`, the map already forgets `m`.

This excludes global conjugacy to an invertible map on the entire legal domain. It does not exclude equivalence on a particular invariant subset. Positive average log absolute determinant along a bounded orbit is possible for a noninvertible map; it is not a contradiction and is not a license to apply an invertible-map dimension formula without further hypotheses.

## A credited symmetric cycle and its exact local stability interval

Symmetric two-cycles of Adam on quadratics are already described by [Bock and Weiss](https://arxiv.org/abs/2210.02070), and in the epsilon-outside convention by [Fong and Yang, Proposition 3.1](https://arxiv.org/abs/2608.20638). The construction below transfers that family to the existing smooth Yogi update. Bock and Weiss put positive epsilon inside the square root; that positive-epsilon convention must not be identified with ours without changing equations.

Assume `0<b<1`, define `A=(1-b)/(1+b)`, `d=A/2`, and take `0<e<d`. Then `a=d-e>0` and

$$z_+=(a,-Aa,a^2),\qquad z_-=(-a,Aa,a^2)$$

form a two-cycle, since `T=0`, `V=a^2` and the position update subtracts `2a`. This asserts the existence of this symmetric family, not a classification of every two-cycle. Put

$$r=a/d\in(0,1),\qquad q=(1-c)ka^2>0.$$

For the reflection `R=diag(-1,-1,1)`, let `L=R*DF(z_+)`. Symmetry gives `DF^2(z_+)=L^2`. Direct differentiation yields

$$L=\begin{pmatrix}
1+2b-2rq & 2b/A & -2(1-q)/A\\
-(1-b)&-b&0\\
qrA&0&1-q
\end{pmatrix}.$$

Its characteristic polynomial is `p(z)=z^3+a1*z^2+a2*z+a3`, with

$$a_1=-b-2+q(1+2r),\quad
a_2=1+2b+q(2br-b-1),\quad a_3=b(q-1).$$

The cycle is locally exponentially stable **if and only if**

```math
\max(0,q_-)<q<q_+,\qquad
q_-={4br-(1-b)^2\over b(1-b+2r)},\quad
q_+={2(1+b)\over1+b+r(1-b)}.
```

Here “exponentially” excludes undecided nonhyperbolic boundary cases. The interval is nonempty, because

$$q_+-q_-={ (1-b)(1+r)((1+b)^2-4br)\over
b(1-b+2r)(1+b+r(1-b))}>0.$$

### Proof of the root criterion used here

For a real monic cubic `p`, set `p*(z)=z^3*p(1/z)`. When `|a3|<1`, its Schur transform

$$Q(z)=[p(z)-a_3p^*(z)]/z=Uz^2+Wz+H$$

has `U=1-a3^2`, `W=a1-a3*a2`, `H=a2-a3*a1`. All roots of `p` are in the open unit disk if and only if all roots of `Q` are there. In one direction apply Rouche's theorem to `p-a3*p*` on the unit circle; the converse follows from `U*p=z*Q+a3*Q*` and the same theorem. A stable monic cubic necessarily has `|a3|<1` from its root product.

For the quadratic, the substitution `z=(1+s)/(1-s)` gives coefficients `U-W+H`, `2(U-H)`, `U+W+H`. The elementary real quadratic left-half-plane criterion therefore gives `U-H>0` and `U+H>|W|`. Consequently the cubic conditions are

$$|a_3|<1,\quad p(1)>0,\quad -p(-1)>0,\quad
1-a_2+a_1a_3-a_3^2>0.$$

In this case they factor as

$$p(1)=2qr(1+b),$$
$$-p(-1)=4(1+b)-2q[(1+b)+r(1-b)],$$
$$1-a_2+a_1a_3-a_3^2
=q[(1-b)^2-4br+b(1-b+2r)q].$$

They give the interval above; `q_+<=2` makes `|b(q-1)|<1` automatic inside it. Eigenvalues of the full two-step derivative are squares of those of `L`, proving the exponential stability statement.

### What is already equivalent to Adam

At this cycle, `V_x=2qa` and `V_v=1-q`. These equal the derivatives of the Adam memory update `V=(1-q)v+q*x^2`. With the same epsilon-outside position update, the **entire local Jacobian is identical**. For `0<q<1` this is ordinary Adam with second-memory coefficient `1-q`; this part must not be promoted as a new mechanism. For `q>1` the equivalent linearization has a negative memory coefficient and is outside standard Adam's parameter range. The lower boundary exceeds one exactly when `b>1/(1+2r)`.

At `q=q_+`, `L` has eigenvalue `-1`, but the full two-step derivative has eigenvalue `+1`. Calling this a period-doubling of the two-cycle would be unjustified. At a positive `q_-`, the unit-modulus pair alone does not prove a nonlinear torus bifurcation. No nonlinear boundary classification is claimed here.

The [structure controls](../results/yogi-tanh/structure/analysis.json) verify the identities symbolically, the Adam Jacobian equivalence, and 1,024 cycle/root tests. They supplement the proof; they are not a novelty certificate.
