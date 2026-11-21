# Three chaotic attractors linked only as a triple

We constructed 80 maps with three coexisting connected solenoidal attractors in Euclidean three-space. Any two of the three compact sets can be separated by an embedded sphere after the third is removed. No sphere in the complement of all three separates any component from another. At least **68 of the triples are mutually inequivalent as subsets of three-space under arbitrary ambient homeomorphisms**, including reflections. Their three internal dynamical systems are the same throughout the family. Here the distinction is a property of the embedding of the limiting attracting sets.

An exhaustive search of 217,690 cyclically reduced balanced pure-braid words with a prescribed positive first generator and at most sixteen blocks supplied witnesses for absolute triple invariants from 1 to 16. Integer lattice polygons supplied the additional 64 cases. Each of the 80 cloud files contains 40,000 consecutive iterates from each of three initial states. Regenerating all 80 in a separate directory reproduced their SHA-256 checksums byte for byte. [Words and map parameters](designs.json) · [Trajectories and SHA-256 digests](clouds.json) · [Independent checks](checks.json)

## Local maps and disjoint trapping tubes

Start with three disk centers at horizontal positions −1, 0 and 1 in a transverse plane. A positive or negative braid generator exchanges an adjacent pair by half a rotation around their midpoint in that direction. Its time parameter on each sector is the quintic 10t³ − 15t⁴ + 6t⁵ for t from zero to one. The two exchanged centers are at distance one throughout, and the third center is at least one from either. The word uses full twists on adjacent pairs, so the three center paths close separately. Write their transverse coordinates as uᵢ(θ), vᵢ(θ).

With a = 0.27 and R = 3, embed the disk of radius one around component i using the chart below. The cylindrical radius always stays positive, so the ambient point determines θ. Disks in different charts at that θ are separated by at least 1 − 2a in the transverse plane. Each chart is a solid torus.

```math
H_i(\theta,x+iy)=\left((R+u_i(\theta)+ax)\cos\theta,\ (R+u_i(\theta)+ax)\sin\theta,\ v_i(\theta)+ay\right).
```

Use the following local map in each chart, with respective degrees 2, 3 and 5. Angles are taken modulo two pi.

```math
F_i(\theta,w)=\left(d_i\theta+0.12\sin(d_i\theta)+0.137,\quad 0.105w+0.55e^{i\theta}\right).
```

The image lies strictly inside the disk since 0.105 + 0.55 < 1. The derivative of the angular map is at least 0.88 times the degree. Two input angles with the same output differ by a multiple of 2π/dᵢ, because the common outer circle map has strictly positive derivative. Distinct branch disks have center distance at least 1.1 sin(π/5), exceeding their combined radii 0.21. Thus each map embeds its trapping torus into its interior. Intersecting its nested images gives a connected attractor. Its dynamics are the inverse limit of the expanding degree-dᵢ circle map: a full backward angular history determines a unique fiber point by contraction, and injectivity recovers that history. Its topological entropy is log dᵢ. These three formulas define one map on a neighborhood of the three compact sets. If a global map of three-space is wanted, extend the charts a little, then multiply the Cartesian formula on each by a smooth cutoff equal to one on the closed trapping torus and constant outside the enlarged tubes. The local attractors are unchanged. The charts and the local formulas are twice continuously differentiable.

## Exact Hausdorff dimensions

Each of the three component dimensions is independent of its surrounding braid. For common fiber contraction ρ = 0.105 and degree d, the full component has Hausdorff dimension 1 + log(d)/log(1/ρ). This gives approximately 1.307547, 1.487450 and 1.714101. Three target angles for each degree were checked through seven inverse levels, reaching 78,125 disjoint disks per checked degree-five slice. [Fiber counts and observed gaps](fiber-sections.json)

On a sufficiently short angular interval, label the d inverse branches continuously and code full backward histories by one-sided sequences in d symbols. Give the sequences distance ρⁿ when their first difference follows n common inverse choices. At fixed angle, two corresponding fiber points are at least gρⁿ and at most 2ρⁿ apart, where g = 1.1 sin(π/d) − 2ρ > 0. Indeed, every common prefix applies a fiber similarity of ratio ρ; the next distinct branches lie in the separated first-level disks. For a fixed history, each past angle at depth j varies with Lipschitz constant at most (0.88d)⁻ʲ on that interval. The weighted sum of its complex exponentials giving the fiber point consequently varies Lipschitz-continuously with the present angle, uniformly in the history.

These bounds make the chart on the attractor locally bi-Lipschitz for the product of angular distance and symbolic distance, with the maximum metric. For the inverse bound, compare two histories first at a common angle and then subtract the uniform Lipschitz change due to their angular difference. A ball of radius comparable to ρⁿ in the sequence space fixes n digits, so its uniform Bernoulli measure is comparable to d⁻ⁿ. Multiplying that measure by ordinary length on the angular interval gives both dimension bounds and the stated value. The physical tube chart and its local inverse have bounded derivatives on compact subcharts.

## Verifying the core links

The two block generators A and B are full twists of strands 1–2 and 2–3 respectively. Every retained word has zero signed exponent sum for each generator. Deleting the middle strand kills both generators. Deleting either outside strand leaves a pure two-strand braid with zero total exponent. Closing these braids therefore gives the two-component unlink in each case. All three pairwise linking numbers vanish.

The first calculation tracks conjugators of meridians through every half-twist, retaining the Magnus coefficients up to degree two after setting the third variable to zero. The coefficient of the ordered term X₁X₂ in the third longitude is Milnor's integer μ₁₂₃ for zero pairwise linking. Factors in the meridian of component three do not affect this coefficient. A second implementation conjugates pairs of meridians through entire full twists using the integral Heisenberg group. It checks the longitude separately for every chosen word. The absolute value is insensitive to the conventions used here for orientation and order. [Magnus implementation](../../../../research/shape_search_20261121/extension/brunnian/search.py) · [Full-twist verifier](../../../../research/shape_search_20261121/extension/brunnian/verify.py) · [Milnor invariants from longitudes](https://msp.org/agt/2014/14-2/agt-v14-n2-p20-s.pdf)

## Passing from tubes to the compact sets

For three zero-pairwise-linking cores with invariant μ, replace core i by a curve in its own solid torus of signed winding wᵢ. The new triple has zero pairwise linking and invariant μw₁w₂w₃. To see this, homotope each new component inside its tube to a wᵢ-fold traversal of its core. Self-intersections are permitted in a link homotopy and the three tubes are disjoint. The generalized Gauss map of the three curves then factors through the product of the three degree-wᵢ circle maps. The primary degrees vanish. Its secondary invariant can be computed by the integral of a primitive of the pulled-back normalized area form wedged with that form. Pulling back this integral multiplies it by the degree of the product, w₁w₂w₃. The secondary invariant is twice Milnor's invariant, giving the asserted law. [DeTurck et al., Theorems A and B](https://arxiv.org/abs/1207.1793)

At common depth n the three solid tori have winding 2ⁿ, 3ⁿ and 5ⁿ in their original tubes. Their core invariant consequently has absolute value 30ⁿ times its initial value. This is nonzero at every depth. Any separating sphere avoiding the three compact attractors would already avoid all three tori at a sufficiently large finite depth, by compactness of the sphere and nesting. Each torus is connected, so this sphere would split that depth's core link, contradicting its nonzero triple invariant. Conversely, each pair of outer cores forms the unlink. Its two disjoint tubular neighborhoods can be shrunk by ambient isotopy toward their cores, so a sphere also separates the original pair of tubes and the two compact sets they contain.

## Distinguishing 68 ambient types

For a positive integer m, remove all factors 2, 3 and 5 and denote the remaining integer by q(m). Suppose an ambient homeomorphism takes a compact triple with initial absolute invariant m to another with invariant m′. Its three connected components can be permuted. For a sufficiently large common n, the images of the three depth-n source tori lie inside their corresponding outer target tori. Their core link has absolute triple invariant m·30ⁿ. Write wᵢ for the signed winding of each image torus in its target tube. By the winding law above, this same core link has absolute invariant m′·|w₁w₂w₃|.

For each target component, a deeper target torus lies inside the corresponding image torus because the target compact component lies in the interior of the image torus. The homomorphisms on first integral homology from this deeper target torus through the image torus into the outer target torus compose to multiplication by a power of its degree dᵢ. Thus wᵢ divides a nonzero power of 2, 3 or 5. Removing all these prime factors from the equality gives q(m) = q(m′). This is a necessary condition even when the homeomorphism permutes components or reverses orientation. The 80 verified examples yield 68 distinct values of q. For instance, Cases 30300001 and 30300007 have q equal to 1 and 7. This distinguishes the two full triples.

The same argument gives infinitely many examples at the fixed internal degrees 2, 3 and 5. The pure-braid commutator of Aᵏ and B has absolute invariant k: its two full-twist blocks act on the class-two longitude by the commutator of the first two meridians to power ±k. Taking k over positive integers coprime to 30 gives distinct q values. The finite dataset above provides the sampled and independently verified cases. The topology of embedded solenoids and the triple-linking invariants are established prior work. Historical priority for this specific family of simultaneous attracting sets remains open. [Jiang et al. on embedded solenoids](https://arxiv.org/abs/math/0611900) · [DeTurck et al. on three-component links](https://arxiv.org/abs/1101.3374)

## Reproduction

Use a fresh external directory named RUN. The first command enumerates the words. The next two sample the three local maps in Cartesian coordinates, save a component label for each state and check the resulting files.

```text
python research/shape_search_20261121/extension/brunnian/search.py --out RUN/designs.json
python research/shape_search_20261121/extension/brunnian/sample.py RUN/designs.json --out RUN/extension/brunnian --manifest RUN/clouds.json
python research/shape_search_20261121/extension/brunnian/verify.py --packet RUN --cache RUN --out RUN/checks.json
python research/shape_search_20261121/extension/brunnian/figure.py RUN/designs.json --cache RUN/extension/brunnian --out RUN/brunnian.png
```
