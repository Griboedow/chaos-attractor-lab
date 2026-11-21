# Higher-order cyclic attractors in the three-axis twist map

Nine parameter vectors have localized attracting sets with an exact number of nondegenerate connected components. The counts cover every integer from 8 through 14. These vectors were screened before commit `6610acc`. The strict trapping and noncollapse proofs are new. All nine original-coordinate point clouds contain 240,000 consecutive float32 positions under the key `points`. The [manifest](manifest.json) gives their cache-relative and absolute paths, SHA-256 hashes and independent numerical exponents.

| Components | Case IDs |
|---:|---|
| 8 | 2200235, 3600210 |
| 9 | 3503658 |
| 10 | 2300491 |
| 11 | 3701267 |
| 12 | 2302642 |
| 13 | 2102331, 2700304 |
| 14 | 2000254 |

For the 14-component example, the canonical start enters the trapping box after 57 steps and its independent finite-time largest exponent is 0.0116554 per iterate. The 12-component example gives 0.0177788. The separate [trajectory checks](checks/) include the full QR spectrum, exponential rates over time windows, additional starts and nearby coefficients. These finite trajectories do not give an asymptotic exponent or an entropy theorem. The certified localized set may coexist with other attractors.

## Basis of the component counts

For each case, the two implementations prove that the convex compact return box T misses every earlier forward image up to iterate m−1 and that the m-th image is strictly inside T. Each twist preserves its axis coordinate and has a global inverse obtained by negating its angle. The affine contraction has positive, nonzero ratio. Consequently the full map F is invertible. The decreasing intersection of the connected compact sets F^(mn)(T) is a nonempty connected set A₀. It is invariant under F^m and attracts a neighborhood under this return. Its first m iterates are pairwise disjoint by invertibility and the earlier-image exclusions. Their union therefore has exactly m connected components and attracts an open neighborhood under F.

The directed verifier evaluates the derivative product with the exact telescoping identity obtained by writing each actual Jacobian Jᵢ as a point matrix Cᵢ plus its enclosed remainder. For Dᵢ=Jᵢ…J₁M, the contribution of remainder i to a final product is propagated by Cⱼ…Cᵢ₊₁ on the left. A mean-value form bounds every point of a dyadic input cell. The saved split tree covers the complete box. A second implementation uses 150-bit Arb matrices, fresh trigonometric enclosures and extra subdivisions where needed. The [proof packets](traps/) identify the exact binary64 coefficients and box frames.

For every case a [384-bit interval Newton check](noncollapse/) proves existence of a solution of F^(2m)(x)=x in the interior of T. Its m-th image is disjoint from its isolating cube. The preconditioned residual maps that closed cube strictly inside itself with Lipschitz constant below one. The preconditioner is an invertible interval enclosure of the inverse point Jacobian. The two orbit points x and F^m(x) lie in every F^(mn)(T) by periodicity and invariance of T under the return. Thus every component contains at least two points and, by connectedness, uncountably many. A numeric canonical start is verified in the trap exactly when the manifest supplies an entry step. The entry for 3600210 remains unverified. The 1,792-bit point propagation proves canonical entry at step 284 for case 2302642 and step 281 for case 3701267.

## Coexisting periodic sinks

The second search screened 1,920 distinct new coefficient vectors in three neighborhoods. The audit regenerated the 175,104 distinct earlier designs and checked full-vector equality. The two searches have zero shared tuples. Of 449 numerically selected candidates, 26 retained representatives have verified attracting periodic orbits and 14 have two distinct attracting periodic orbits at the same coefficients. Examples include (4,348), (29,232), (4,145), (4,116), (2,96) and (10,24). The [index](cycles.json) and [portable certificates](cycle-certificates/) give the exact coefficients, proposed centers and linear frames.

For each phase the verifier encloses the one-step Jacobian on the complete ellipsoid, checks positive definiteness of q²I−JᵀJ in the successor frame by Sylvester's criterion, and verifies that the normalized center shift is below (1−q) times its radius. All neighboring maps are strict contractions in the successive charts. Every pair of phase ellipsoids has a strict native-coordinate separator. The contraction theorem gives a unique attracting orbit of the stated least period. Distinct least periods guarantee different attracting orbits for the same parameters. The source uses exact dyadic centers for strongly anisotropic charts.

## Sharper lower dimension bound on the earlier five-component set

Over the same closed ten-dimensional coefficient block of radius 2⁻²⁴ around case 2503155, the [150-bit Arb computation](dimension-finest.json) gives **dimH A > 1.041063**. The existing upper bound is 1.45. The [earlier proof](../../maps/dimension-proof.md) supplies the symbolic slabs, exact slice and arbitrary past, uniform slice content and native separation. Only the inverse Lipschitz estimate changes. Intersect the two previously certified native hulls for each phase. They are convex and contain both endpoints and the segment joining them. Bound each average inverse Jacobian in arbitrary invertible phase-specific bases on 32³ subcells per phase, compose in the inverse order 0,4,3,2,1, and bound the product in a positive weighted maximum norm. The output gives L < 21,425,228. The fixed start/end chart makes this the same norm at each return. Equivalence of this norm and native Euclidean distance contributes a single constant independent of past length. Thus the slice argument gives dimH A ≥ 1+log(2)/log(L). The frame optimization only proposes bases. Every inverse, enclosure and final rational inequality is checked by Arb.

## Reproduction

Run from the repository root with Python, NumPy, SciPy, Numba and python-flint. Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and `SHAPE_SEARCH_THREADS=1`. The following PowerShell commands assume the shared cache location used in this packet.

```powershell
$s='research/shape_search_20261121/extension/maps'
$r='results/shape_search_20261121/extension/maps'
$c='D:/Codex/chaos-attractor-cache/shape-search-20261121/extension/maps'
$old='D:/Codex/chaos-attractor-cache/shape-search-20261121/maps'
python "$s/multistable.py" $c --count 640
python "$s/cycle_tubes.py" $c $r
python "$s/noncollapse.py" "$r/traps" "$r/noncollapse"
python "$s/adapted_dimension.py" "$r/dimension-finest.json" --frames "$r/dimension.json" --grid 5
python "$s/../../maps/validate.py" --cache $c --results $r --preview "$old/review-cyclic/preview.json" --proposals "$r/proposals.json"
python "$s/../../maps/validate.py" --cache $c --results $r --preview "$old/review-cyclic/preview.json" --proposals "$r/proposals-supplement.json"
python "$s/package.py" $c $r
```

To independently reproduce the trapping checks directly from the tracked packets, run the following from the repository root. Each replay also checks the stored input digest.

```powershell
python -c "import sys,json,pathlib; sys.path.insert(0,'research/shape_search_20261121/extension/maps'); from replay_traps import replay; from telescoped_return import run; ps=[json.loads(f.read_text()) for f in pathlib.Path('results/shape_search_20261121/extension/maps/traps').glob('*.json')]; print([(p['box']['case'],run(p['box'],9000)['nodes'],replay(p['box'],p['directed'])['accepted_leaves']) for p in ps])"
```
