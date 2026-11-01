# A Whitney fold in the attracting spatial-circle RAAR set

## Exact setting and conclusion

Use the exact binary64 parameter lifts of the R1051 instance:

`[beta, r, theta, cx, cy, cz] = [0.9708823579074815, 1.9886915469542146, 0.9028818255996501, -0.1435409113764763, 0.5745496526360512, -0.7745124697685242]`.

Let A be the unit circle in the xy plane and B the circle
`c + r(cos(t), cos(theta)sin(t), sin(theta)sin(t))`. Both projectors are nearest
projections onto the circumferences. On their regular domain,

`F(x) = beta*x + beta*P_A(2*P_B(x)-x) + (1-2*beta)*P_B(x)`.

The earlier frozen trapping certificate constructs compact Q with
`F(Q) + [-2^-35,2^-35]^3 ⊂ Q` and strictly positive projector-axis distances.
Thus F is smooth on a neighborhood of Q and
`K = intersection_{n>=0} F^n(Q)` is a compact invariant set attracting a
neighborhood of Q. Its already certified topological entropy is at least
`log(2)/49`. This note uses that completed certificate as an explicit dependency.

**Theorem.** K contains a corank-one critical point of F at which the
one-dimensional kernel of DF is transverse to `det DF = 0`. Consequently F
has a smooth Whitney fold at a point of K. At least one such point lies in

`[0.20984,0.22744] × [-2.86973,-2.84328] × [3.86195,3.88395]`.

There is no conjugacy by an ambient homeomorphism, with the conjugacy identity
holding on an open neighborhood of K, to a map locally injective at every
point of the image of K. This excludes diffeomorphisms and local
diffeomorphisms under that precise notion of equivalence. It does not establish
noninjectivity of the restriction `F|K`, exclude conjugacies defined only on K,
or distinguish the system from all previously described noninvertible maps.

## A thin self-covering window inside Q

Write `N = c0 + S[-1,1]^3`, where `S = B0 diag(r0)`. The center and basis are
unchanged from the earlier window N0; only its two stable widths are reduced:

```
c0 = [0.005907460086282232, -2.444586561951587, 4.139526201112402]
B0 = [[-0.2323036470597354, -0.09754348589931008, 0.9602881626052217],
      [ 0.8662592186016529,  0.7569218185352454,  0.2673319366440987],
      [ 0.44230078199162143,-0.6461846709601738,  0.07987665748141355]]
r0 = [0.0009880210333907218, 0.00001, 0.000001]
```

All listed scalars denote exact binary64 lifts, and their products in S are
mathematical products. Exact comparison with the old widths gives
`N ⊂ N0 ⊂ Q`. The matrix is invertible. Let `G=F^7` and express G in N's
normalized coordinates as `g=(g1,gs)`. The computations establish:

1. On the left face `u=-1`, `g1>1`; on the right face `u=1`, `g1<-1`.
2. For every point of `[-1,1]^3`, either `g1<-1`, `g1>1`, or both stable
   coordinates of g lie strictly between -1 and 1.

The vector-projector implementation at 70 digits covers the volume with 60
leaves; the scalar automatic-derivative implementation at 100 digits uses 41.
Each exit face passes as one leaf in each implementation. These are the actual
inequalities used below; no cone condition or smooth unstable-manifold theorem
is assumed.

## Elementary crossing-continuum lemma

Call a continuous path in N a crossing path when its endpoints lie on the two
opposite u faces. The two conditions above imply that the G-image of any
crossing path contains a crossing subpath in N. Indeed, its unstable coordinate
starts beyond one face and ends beyond the other. Select a subinterval from the
last crossing of the starting face before the first crossing of the opposite
face. Its unstable coordinate stays in [-1,1]. Condition 2 puts both stable
coordinates strictly inside their intervals throughout that subpath, including
its endpoints. Reversing the parameter when needed gives the desired path.

Start with the centerline crossing path P0. Inductively choose a crossing
subpath `P_(n+1) ⊂ G(P_n) ∩ N`. Then each Pn is nonempty, compact and connected,
meets both closed faces, and belongs to `G^n(N) ⊂ F^(7n)(Q)`. The hyperspace of
nonempty compact subsets of compact N is compact in the Hausdorff metric, so a
subsequence with indices tending to infinity converges to a compact nonempty
set C. A Hausdorff limit of compact connected sets is connected: a separation
of the limit into two disjoint compact pieces would have disjoint small
neighborhoods which both meet, and together contain, all sufficiently late Pn,
contradicting their connectedness. Since each Pn meets both closed faces, C
does also.

For every integer m, all sufficiently late Pn lie in `F^m(Q)`, because
`F^(7n)(Q) ⊂ F^m(Q)` whenever `7n>=m`. These sets are compact and hence closed.
Thus `C ⊂ intersection_m F^m(Q) = K`. This proves the existence of a connected
subset of `K ∩ N` crossing both faces. It does not assert that C is a smooth
curve or that the numerical centerline lies in K.

## A critical point is forced inside K

Define the continuous real function `h(x)=det DF(F^42(x))` on N. Regularity
follows from `N ⊂ Q`, forward invariance of Q and its positive axis bounds.
The whole-face interval calculations give

`h(left face) < -0.0311`, and `h(right face) > 0.1744`.

These conservative displayed bounds follow from the 70-digit reconditioned
enclosure. The 100-digit scalar/automatic-derivative realization independently
confirms the signs. An unreconditioned mean-value calculation also passes,
using four left-face leaves and one right-face leaf, with weaker margins.

The connected set C meets both faces, so h(C) is an interval containing values
of opposite signs. There is x in C with h(x)=0. Then `z=F^42(x)` belongs to K
by invariance, and `det DF(z)=0`. This proves membership, rather than merely
proximity to a finite numerical cloud.

## Certifying fold type and location

Subdivide the entire normalized cube of N. For each source cell, enclose its
42-step image X. Accept the cell if either:

* `det DF(X)` excludes zero; or
* two selected rows of DF have a cross product v with a component bounded
  away from zero, and `gradient(det DF)(X) dot v(X)` excludes zero.

At a critical point in the second case, the nonzero minor gives rank at least
two, and the zero determinant gives rank exactly two. The cross product spans
the kernel. Its nonzero contraction with the determinant gradient proves both
regularity of the critical hypersurface and kernel transversality. The standard
smooth-fold criterion therefore applies. The source and target changes of
coordinates in the fold normal form are different; they are not a dynamical
conjugacy.

The primary calculation accepts 478 leaves: 452 regular and 26 fold-capable.
The second calculation accepts 283: 270 regular and 13 fold-capable. It uses
scalar circle coordinates, nested reflections and forward first/second
automatic derivatives; the determinant gradient is differentiated by its six
permutation products. All source cells are covered, not just samples near a
guessed zero. The earlier forced point z must lie in the union of the
fold-capable image boxes. Their primary coordinate hull lies strictly inside
the outward decimal box in the theorem. We do not locate one unique fold or
claim that every point of this hull is critical.

## Enclosure and audit assumptions

The reconditioned mean-value method maintains
`X ⊂ C + A*delta + B*r`, with interval vectors and point matrices. Choose
arbitrary finite point matrices Anew and Bnew, the latter verified invertible
by interval arithmetic. The mean-value theorem gives

```
F(X) ⊂ F(C) + Anew*delta + E,
E = (DF(X)*A-Anew)*delta + DF(X)*B*r.
rnew = enclosure(Bnew^-1 * E).
```

The next representation is `F(C)+Anew*delta+Bnew*rnew`. A numerical QR
factorization selects Bnew; neither exact orthogonality nor numerical accuracy
of that selection is assumed. The interval inverse and mismatch term account
for the actual chosen matrices. All segment domains are enclosed and positive
projector-axis bounds are checked. The proof relies on the stated mpmath
interval arithmetic and inspected formulas, not on sampled error estimates.
It is not a proof-assistant formalization.

The second map/derivative realization shares the enclosure driver and interval
library. It is therefore not an entirely independent software stack. The
separate accounting auditor verifies exact window containment, dependency
manifests, full binary subdivision trees, strict inequalities and source
identities. It does not itself recompute every enclosure. Failed and incomplete
development attempts are retained and excluded from the proof.

## Local noninjectivity and the equivalence boundary

At a Whitney fold, smooth source and target coordinates give the local form
`(s,t,w) -> (s,t,w^2)`, up to signs and permutations. Hence F is not injective
on any open neighborhood of z. If an ambient homeomorphism H conjugated F to
T on an open neighborhood of K and T were locally injective at H(z), choose
a sufficiently small open neighborhood of z on which the conjugacy identity
holds and whose H-image lies in that injectivity neighborhood. Injectivity of
H and T would force F to be injective there, a contradiction.

The nearby pairs witnessing a fold need not both belong to K. Consequently
this argument does not show `F|K` is noninjective and does not prohibit an
abstract topological conjugacy of restrictions to invariant sets. A smooth
fold, noninvertible chaos and computer-assisted covering proofs are all
established concepts. This theorem provides a specific obstruction to ambient
local equivalence, not a classification of historically new attractors.
