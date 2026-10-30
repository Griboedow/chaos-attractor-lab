# Eight-node candidate graph for the motor saddle

The two observed departure permutations generate a formal graph with **eight nodes and sixteen directed edges**. There are two nodes at each of four levels, and every node has an edge to both nodes at the next level. The previously displayed four-saddle circuit is one of sixteen simple four-cycles in this graph.

These are exact finite-group computations. Identifying all eight nodes with distinct periodic orbits requires the explicit dynamical assumptions in the lemma below. The edges encode the finite-time numerical destinations in the [postcritical study](motor-postcritical.md); their existence as infinite-time heteroclinic connections and attraction of the graph remain unproved. Graph cycles are sequences of candidate connections, not additional periodic trajectories.

![Eight-node candidate graph](../results/motor-network-combinatorics/figure-layout-review/candidate-graph.png)

## Row action and cosets

All row indices are zero-based. Write the physical row state as an eight-entry list of triples and use the observed convention

\[
(q\mathbin{\cdot}a)_j=q_{a_j},\qquad
(a\circ b)_j=a_{b_j}.
\]

Thus \((q\cdot a)\cdot b=q\cdot(a\circ b)\). At the reference saddle the paired rows are \((2,3)\) and \((6,7)\), while rows \(0,1,4,5\) are singletons. Define

\[
H_0=\langle(2\ 3),(6\ 7)\rangle,\qquad
h=(0\ 1)(2\ 3)(4\ 5)(6\ 7),\qquad
H=\langle H_0,h\rangle.
\]

The known half-period relation supplies the orbital symmetry \(h\). The numerical outgoing permutations are

```text
g_plus  = [5,4,6,7,2,3,0,1]
g_minus = [5,4,7,6,2,3,1,0]
```

If \(H\) is the orbital stabilizer within the generated group, orbit labels are **left cosets** \(Ha\). Transporting the observed connection by \(a\) gives

\[
Ha\longrightarrow H(g_+\circ a),\qquad
Ha\longrightarrow H(g_-\circ a).
\]

The computation checks all eight representatives of every coset: the unordered two-target set is unchanged. A plus/minus name can change with the representative, so the figure gives unlabelled edges rather than a global sign convention.

## Exact group and graph

Let \(s_j=(2j\ 2j+1)\). The pair-flip kernel

\[
K=\langle s_0,s_1,s_2,s_3\rangle\cong C_2^4
\]

has order sixteen. An order-four lift of the pair rotation is

```text
R = [4,5,6,7,2,3,0,1]
pair action of R = [2,3,1,0]
```

It permutes the four flip coordinates cyclically. Exact enumeration verifies

\[
G=\langle H,g_+,g_-\rangle
 =K\rtimes\langle R\rangle\cong C_2^4\rtimes C_4,
\quad |G|=64,\quad |H_0|=4,\quad |H|=8.
\]

The set \(\{kR^j:k\in K,\ 0\leq j<4\}\) is exactly \(G\). Moreover, \(g_+=R\circ s_0\) and both \(g_+^4\) and \(g_-^4\) equal \(h\); each generator has permutation order eight and orbit-coset order four.

| Level | Nodes | Expanding pair | Contracting pair | Targets from either node |
|---|---|---|---|---|
| 0 | L0a, L0b | 6,7 | 2,3 | L1a and L1b |
| 1 | L1a, L1b | 2,3 | 4,5 | L2a and L2b |
| 2 | L2a, L2b | 4,5 | 0,1 | L3a and L3b |
| 3 | L3a, L3b | 0,1 | 6,7 | L0a and L0b |

The pair labels come from transporting the reference numerical Floquet modes. The two nodes at a level have the same cluster memberships but different orbit cosets. Each node has in-degree two and out-degree two. Exhaustive enumeration gives sixteen simple directed cycles of length four and eight of length eight, with cycles identified up to cyclic rotation. There are no other simple-cycle lengths.

The previously selected positive-generator circuit is

```text
L0a -> L1b -> L2b -> L3b -> L0a.
```

The negative-generator powers give

```text
L0a -> L1a -> L2a -> L3b -> L0a.
```

These statements concern the two supplied destination permutations and their symmetry images. They do not establish the entire dynamical network under all of \(S_8\).

## Conditional lemma: eight distinct images without full isotropy

Assume a nonconstant periodic orbit exists in the permutation-equivariant eight-row ODE, its six cluster triples are pairwise distinct at one phase with precisely the two pairs above, and \(h\) maps this orbit to itself with a time shift. No least-period certificate is assumed here; a nonconstant periodic orbit has a mathematical least positive period.

Let \(S\) be its full orbital stabilizer in \(S_8\). Distinct cluster triples make the point stabilizer exactly \(H_0\). By uniqueness of ODE solutions and equivariance, a permutation fixing one point fixes the entire trajectory. Therefore the kernel of the orbital phase-shift homomorphism is \(H_0\), and \(S/H_0\) is a finite subgroup of phase rotations, hence cyclic. In particular, \(S\) normalizes \(H_0\).

The exact permutation calculation gives

\[
N_G(H_0)=K,\qquad K/H_0\cong C_2^2.
\]

One can also see the normalizer directly: only the identity power of the pair rotation preserves the set of the two paired base pairs; every flip in \(K\) does. Consequently \((S\cap G)/H_0\) is a cyclic subgroup of an exponent-two group and has order at most two. The known \(h\notin H_0\) supplies a subgroup of order two. It follows that

\[
S\cap G=H,\qquad [G:S\cap G]=8.
\]

Thus the stated dynamical assumptions suffice for eight distinct periodic-orbit images **within \(G\)**, even if additional orbital symmetries exist outside \(G\). This does not identify the full \(S_8\) stabilizer and does not prove any connecting orbit. The computer calculation proves the finite-group ingredients; an orbit certificate must supply the dynamical assumptions.

## Numerical phase check

The separate phase audit uses the frozen [six-cluster cycle](../results/motor-postcritical/six-cluster-cycle-rho-3p268/analysis.json) at \(\rho=3.268\), with period approximately \(13.633145767830683\). It compares all 64 permuted initial physical row states against the stored orbit, using periodic cubic interpolation of 1,001 samples. For each permutation it refines every local minimum found on the 1,000-phase grid and also retains the original grid candidates.

The eight elements of \(H\) have maximum Euclidean residual \(1.3325\times10^{-13}\). The smallest residual among the other 56 permutations is **0.0714879262** in the 24 physical coordinates \((c_j,s_j,X_j)\). This is numerical support for the orbit identification, not an interval exclusion over all phases. Interpolation forces the final stored sample equal to the first; the unmodified endpoint discrepancy is recorded. No fresh ODE integration or GPU work was performed.

## Reproduction and evidence

```powershell
../chua-attractor/.venv/Scripts/python.exe -B research/motor_network_combinatorics.py --output results/motor-network-combinatorics-recheck
```

The output directory must not already exist. The complete [exact group and graph](../results/motor-network-combinatorics/exact-graph.json) includes all 64 permutations, eight cosets, sixteen edges, all simple cycles and normalizer checks. The [phase audit](../results/motor-network-combinatorics/phase-audit.json) retains all 64 cases. Each run records source and input SHA256 hashes and archives its executed source.

The initial figure placed some arrows through intermediate node glyphs. Its original bundle is retained. The [figure-layout review](../results/motor-network-combinatorics/figure-layout-review/analysis.json) changes only the drawing layout and is the figure displayed here; its exact graph and phase audit are byte-identical to the first run. Both figures were inspected, and the second resolves the routing ambiguity. The two runs took approximately 1.20 and 1.09 CPU wall seconds, respectively.
