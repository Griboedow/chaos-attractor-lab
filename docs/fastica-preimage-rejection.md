# A rejected intrinsic-collision candidate in FastICA

**Result, 2026-11-01.** A proposed collision in the structured FastICA map is
real, but its extra preimage is outside the attracting set. Exact rational
checks certify two distinct points in the trapping torus with the same image
and exclude a whole neighborhood of one of them from the image of that torus.
This candidate therefore supplies no evidence that the map restricted to the
attractor is noninjective.

[Construction and portraits](fastica-tripling.md) ·
[Ambient local-noninjectivity theorem](fastica-tripling-structure.md) ·
[Certificate and reproduction](../results/fastica-preimage-rejection/README.md)

## Statement

Use the exact model at $`\epsilon=1/25`$, with homogeneous cubic numerator
$`h`$, normalization $`N(v)=v/\|v\|`$ and $`F(w)=N(h(w))`$. The
[proved trapping region](fastica-tripling.md#strict-trapping) is

```math
T=\{(A,B,C,D)\in S^3:C^2+D^2\leq\tfrac4{25}(A^2+B^2)\},
\qquad K=\bigcap_{n\geq0}F^n(T).
```

The following entries specify an **exact decimal rational vector**, obtained
by rationalizing one finite numerical sample:

```text
v = (-0.38715568266959405, -0.9020845179494923,
     -0.19040449850479912,  0.01000633807187778)
p = v / ||v||
```

Define a closed box $`X`$ in the chart $`q=N(a,-1,c,d)`$:

| Coordinate | Exact decimal center | Half-width |
| --- | --- | --- |
| $`a`$ | -0.43931892 | $`10^{-7}`$ |
| $`c`$ | -0.29018779 | $`10^{-7}`$ |
| $`d`$ | -0.03294978 | $`10^{-7}`$ |

There is exactly one $`x_*=(a_*,c_*,d_*)\in X`$ for which

```math
q_*=N(a_*,-1,c_*,d_*),\qquad F(q_*)=F(p).
```

Both points lie strictly inside $`T`$, and $`q_*\neq p`$. Moreover, the
entire normalized box $`N(X_1,-1,X_2,X_3)`$ is disjoint from $`F(T)`$.
Consequently $`q_*\notin K`$. Membership of $`p`$ in $`K`$ is not claimed:
a rationalized finite trajectory sample is not an exact invariant-set witness.

## The collision exists

Write $`h_0=h(v)`$ and $`h_x=h(a,-1,c,d)`$. The three polynomial equations are

```math
f_i(x)=(h_x)_i(h_0)_B-(h_x)_B(h_0)_i=0,
\qquad i\in\{A,C,D\}.
```

At the exact center $`x_0`$ of $`X`$, invert $`Df(x_0)`$ over the rational
numbers and call the result $`C_0`$. The fixed points of
$`G(x)=x-C_0f(x)`$ are exactly the zeros of $`f`$.
Exact interval evaluation establishes

```math
\sup_{x\in X}\|I-C_0Df(x)\|_\infty
<0.000226<1,\qquad G(X)\subset\mathrm{int}(X).
```

The smallest coordinate margin in the self-map inclusion exceeds
$`9.70\times10^{-8}`$. These displayed bounds are weakened decimal bounds;
the certificate retains the deciding rational inequalities. The contraction
mapping theorem gives existence and uniqueness in the closed box.

The cubic component $`(h_0)_B`$ and the entire interval $`(h_x)_B`$ are
negative. Thus the common projective image has the same sign and proves
$`F(q_*)=F(p)`$ on the sphere, not just after identifying antipodes.

Exact substitutions also give strict tube margins above $`0.1178`$ for $`v`$
and $`0.1055`$ for every raw vector in the proposed box. In the chart with
$`B=-1`$, the $`C`$ coordinate of $`p`$ lies more than $`0.0791`$ above
the upper face of $`X_2`$. Hence the two normalized points are distinct and
are not antipodal.

## The extra preimage cannot have a predecessor in the torus

Every point of $`T`$ has a nonzero planar component. After rescaling and,
if necessary, changing its overall sign, it has one of two forms:

```math
(1,t,u,v)\quad\text{or}\quad(t,1,u,v),\qquad
|t|\leq1,\quad |u|,|v|\leq\tfrac35.
```

These boxes cover the torus because
$`u^2+v^2\leq(4/25)(1+t^2)\leq8/25<9/25`$. The antipodal change is harmless:
the homogeneous odd numerator changes sign, while the equations for
parallelism test both target signs.

For a target $`(a,-1,c,d)`$, a necessary condition for a preimage is

```math
h_A+a h_B=0,\qquad h_C+c h_B=0,\qquad h_D+d h_B=0.
```

The checker bisects the two covering boxes using exact rational endpoints.
It rejects a leaf only if it lies strictly outside the torus inequality or
one of these equations has an interval range excluding zero for **every**
target in $`X`$. Both children of every split are processed.

| Exhaustive accounting | Count |
| --- | --- |
| Initial charts | 2 |
| Binary splits | 43 |
| Processed boxes | 88 |
| Leaves outside the torus | 4 |
| Leaves excluded by an equation | 41 |
| Maximum split depth | 7 |
| Unresolved leaves | 0 |

All 45 leaf addresses and their exact exclusion witnesses are retained.
An address audit reconstructs the complete binary partitions independently
of the search order. Since $`F`$ is defined throughout $`T`$, the absence of
parallel numerators excludes the normalized target box from $`F(T)`$.
The already proved inclusion $`K\subset F(T)`$ then excludes $`q_*`$ from $`K`$.

## Search scope and consequence

The exploratory search used 16 targets from one binary64 trajectory after
10,000 warm-up steps, spaced by 137 iterations. Each target received 27
least-squares starts from three approximate angular branches and a 3 by 3
transverse grid, plus its known predecessor as a recovery control: 448 solves
in total. Fifteen targets returned one accepted preimage each; the target
with zero-based index 13 returned two. These solver outcomes are not a
complete root count. The extra root at index 13 motivated the certificate.

The proof uses the exact rational proposal above, not a claim that numerical
iteration or least squares is exact. Its interval operations and polynomial
identity checks are reproducible, but the argument is not formalized in a
proof assistant. The arithmetic replay and partition audit reuse the stated
interval implementation; they are not independent implementations of interval
analysis.

The previous ambient local-noninjectivity theorem remains valid. This result
rejects one tempting route from that theorem to an intrinsic solenoid
obstruction. It proves neither injectivity nor noninjectivity of $`F|_K`$.
The known circle-expansion mechanism and the unresolved historical priority
of the finite-data construction retain their
[earlier assessment](fastica-tripling.md#prior-art-and-scope). No new attractor
class or uniqueness probability follows from this check.
