# Full permutation isotropy of the certified motor saddle

The certified six-cluster saddle has **exactly eight orbital symmetries in $S_8$**. It therefore has **5,040 distinct labeled periodic-orbit images**, counted up to time shift. These are symmetry copies of the same saddle, with the same least period and Floquet classification, not newly discovered attractors.

Extending the two numerically observed departure rules to all row permutations gives a formal directed graph with **5,040 nodes, 10,080 edges and 630 components of eight nodes each**. The orbit count uses the interval existence, spectrum and least-period certificates. The edge count is exact permutation combinatorics; it does not certify any infinite-time connecting trajectory or the structure of the entire dynamical network.

## Certified inputs

The [saddle certificate](motor-saddle-certificate.md) concerns the unchanged eight-row model at exact $\rho=817/250$, $r=1+(19/20)/\sqrt{32001/500000}$ and $b=1$. The clusters, using zero-based row labels, are

```text
(0), (1), (2,3), (4), (5), (6,7).
```

The [combined existence and spectrum record](../results/motor-saddle-proof/similarity-certificate.json) separates every pair of the six initial row triples and certifies a nonconstant orbit. Its full-period contrast blocks satisfy:

| Pair split by a contrast | Certified full-period spectrum |
|---|---|
| $(2,3)$ | All three multiplier moduli are below $0.9$ |
| $(6,7)$ | One real multiplier lies in $(1.12,1.13)$; the other two moduli are below $0.9$ |

The [least-period certificate](../results/motor-saddle-proof/least-period.json) revalidates the same shooting root inside the original uniqueness box and proves

```math
T=13.6331457678306460232747\;\mathbin{\pm}\;9.71\times10^{-24}
```

is the least orbit period. It also proves that the scalar signal $X_2$ has least period $T/2$. Since the complete row-2 state is fixed by the half-period row permutation, its least period is likewise $T/2$.

The certified half-period relation is represented on full rows by

```math
h=(0\ 1)(2\ 3)(4\ 5)(6\ 7).
```

Swapping the two members of either synchronized pair has no effect on the orbit. The reduced six-cluster permutation is $(1,0,2,4,3,5)$, exactly the one used in the shooting certificate.

## The full orbital stabilizer

The full-row ODE is $S_8$-equivariant: every row obeys the same equations, and row permutations preserve the arithmetic mean and the zero-sum displacement constraint. Use the pullback action

```math
(q\cdot a)_j=q_{a_j},\qquad
(a\circ b)_j=a_{b_j}.
```

Here $q_j$ denotes the finite physical row triple; all permutations act on all three entries of a row together. Let $\Gamma$ be the certified orbit, and let

```math
S=\{a\in S_8:\Gamma\cdot a=\Gamma\}
```

be its full orbital stabilizer. The distinct cluster triples give the pointwise stabilizer

```math
H_0=\langle(2\ 3),(6\ 7)\rangle,\qquad |H_0|=4.
```

By uniqueness and equivariance, a permutation fixing one point of the orbit fixes its full trajectory. Conversely, distinct row states cannot merge at a later phase. Thus the same pointwise stabilizer applies at every phase.

Each $a\in S$ gives a unique phase shift $\vartheta(a)\in\mathbb R/T\mathbb Z$. This is a homomorphism because the action commutes with the flow, and its kernel is exactly $H_0$. In particular, $H_0$ is normal in $S$.

An orbital symmetry must therefore permute the two two-row orbits of $H_0$, namely the paired sets $\{2,3\}$ and $\{6,7\}$. It could initially appear to exchange them. Such an exchange would carry the variational contrast equation of one pair to that of the other, with a time shift and an invertible change of contrast coordinates. Full-period monodromy matrices at different phases are similar. Exchanging the pairs would consequently force the two full-period contrast spectra to coincide.

The certified spectra exclude this: one block is strictly stable, while the other has a multiplier larger than $1.12$. Every $a\in S$ therefore preserves each paired set individually.

In particular, the row-2 scalar signal must obey

```math
X_2(\tau+\vartheta(a))=X_2(\tau).
```

Its least period is $T/2$, so the only possible orbital phases modulo $T$ are zero and $T/2$. The zero-phase fibre has four elements, precisely $H_0$. The known half-period permutation $h\notin H_0$ supplies the other fibre, $H_0h$, also with four elements. Therefore

```math
S=H=\langle H_0,h\rangle,\qquad |S|=8.
```

This establishes the full stabilizer, strengthening the earlier [normalizer argument within the smaller generated group](motor-network-combinatorics.md). The orbit-stabilizer formula now gives

```math
\#\{\text{distinct labeled orbit images under }S_8\}
=[S_8:H]=\frac{8!}{8}=5040.
```

“Distinct” here means distinct subsets of the labeled phase space after allowing time shifts. They form one orbit class after quotienting by row permutations, and their dynamics are conjugate by symmetry. No equivalence of different spatial patterns in an actual axoneme is asserted.

## The formal full-permutation departure graph

The two observed finite-time destination permutations are

```text
g_plus  = [5,4,6,7,2,3,0,1]
g_minus = [5,4,7,6,2,3,1,0].
```

Their interpretation and numerical limitations are recorded in the [eight-node graph report](motor-network-combinatorics.md). Set

```math
G=\langle H,g_+,g_-\rangle\cong C_2^4\rtimes C_4,
\qquad |G|=64.
```

For the stated row action, orbit nodes are **left cosets** $Ha$ with $a\in S_8$. Transporting a departure by $a$ gives the two edges

```math
Ha\longrightarrow H(g_+\circ a),\qquad
Ha\longrightarrow H(g_-\circ a).
```

The unordered destination set is independent of the representative chosen in $Ha$. A global plus/minus label is not assumed: changing a node representative can exchange the local sign labels.

Every edge from a node inside a left coset $Ga$ stays inside $Ga$, because its multiplier belongs to $G$. The nodes there are the eight $H$-cosets within $Ga$. Right translation by $a$ is a graph isomorphism from the reference eight-node graph to this component. The reference graph is strongly connected, so these sets are exactly the strongly connected components of the specified formal graph. There are

```math
[S_8:G]=\frac{40320}{64}=630
```

components. Each has eight nodes and sixteen edges; every node has in-degree two and out-degree two. Hence the complete graph has 5,040 nodes and 10,080 directed edges.

This decomposition describes only the edge family generated from the two supplied numerical departures. It does not rule out other dynamical connections, establish convergence to these targets, prove attraction of a component, or show that an actual global invariant network decomposes into these 630 components. It makes no claim of 630 new attractors.

## Exact enumeration and provenance

The [extension script](../results/motor-network-combinatorics/full-permutation-extension.py) enumerates all 40,320 row permutations using integer tuples. It partitions them into all 5,040 $H$-cosets and 630 $G$-cosets, checks all 40,320 representative choices for the two-target rule, verifies every in-degree and out-degree, and checks reachability and reference-graph isomorphism in every component.

The [complete result](../results/motor-network-combinatorics/full-permutation-extension.json) stores every node representative, its two destinations, component membership and reference-node label, together with every component representative. It also checks the archived least-period and spectral theorem gates, exact dyadic sign/margin bounds, certificate linkage and source-snapshot hashes. These checks consume the existing interval certificates; they do not repeat their integrations or constitute a separate formal verification of the interval engine.

```powershell
python results/motor-network-combinatorics/full-permutation-extension.py --output outputs/full-permutation-extension-recheck.json
```

Only the Python standard library is needed for this enumeration and certificate-record check. Choose a new output path. The record includes the script, report and input-certificate SHA256 hashes. No new trajectory simulation is used.
