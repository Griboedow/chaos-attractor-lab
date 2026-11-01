# RAAR candidates: prior art and equivalence audit

**Verdict as of 1 November 2026:** this is a reproducible study of particular spatial-circle configurations in an existing algorithm. We have not located an exact published match for the selected cases. Historical uniqueness, a new topological attractor class and inequivalence to every previously described map remain **unestablished**.

The working labels “serrated arch” and “layered loop” describe views of actual states. They are not mathematical names or priority claims. Changing a parameter, camera or coordinate system does not by itself create a new discovery.

## Closest known work

| Primary source | Already known | Difference and review limit |
|---|---|---|
| [Luke, Relaxed averaged alternating reflections for diffraction imaging (2004 preprint / 2005)](https://arxiv.org/abs/math/0405208) | The exact RAAR iteration, its relaxation and connections to related projection algorithms. | Our algorithm is this published algorithm. Abstract, introduction, Sections 2.1–2.3 and equations 10, 11, 14 inspected. |
| [Luke and Martins, Convergence Analysis of the Relaxed Douglas–Rachford Algorithm (2020)](https://doi.org/10.1137/18M1229638) | Explicit RAAR examples with **two planar circles**, including intersecting, separated, nested and concentric cases; general fixed-point and local convergence theory. Example 3.17 already supplies a fixed-point-free circle/point case. | Our selected circles are noncoplanar in R³. Section 5 and selected earlier examples inspected; general proofs not fully re-audited. Fixed-point freedom in RAAR is not new in general. |
| [Borwein et al., Dynamics of the Douglas–Rachford method for ellipses and p-spheres (2017 online; 2018 volume)](https://arxiv.org/abs/1610.03975) | Rich periodic domains, long transients and earlier circle-line singular-set phenomena. | β=1 Douglas–Rachford and different constraints; our bounded β<1 spatial examples require separate analysis. Most descriptive text and five figure pages inspected, not every proof/figure. |
| [Elser, Rankenburg and Thibault, Searching with iterated maps (2007)](https://doi.org/10.1073/pnas.0606359104) | Strange-attractor search dynamics in projection/difference-map methods explicitly discussed. | Excludes any claim to discover chaos in constraint-search algorithms. Publisher-indexed abstract and body excerpts inspected; full paper remains a review gap. |
| [Marchesini, A unified evaluation of iterative projection algorithms for phase retrieval (2007), institutional figure 5](https://digital.library.unt.edu/ark:/67531/metadc877846/m1/10/) | Elementary circle examples and algorithmic failure modes. | Two circles form a **union in one constraint**, with a line as the other; not our two separate spatial circles. Figure caption and indexed text inspected; full paper review incomplete. |
| [Dinh, Jansen and Luke, Cyclic Relaxed Douglas–Rachford Splitting for Inconsistent Nonconvex Feasibility (2025 preprint / 2026)](https://arxiv.org/abs/2502.12285) | Expanded RAAR formula, fixed-point linear combinations, shadow relations and local convergence theory. | Section 3.1 and selected definitions inspected. Our expansion and general fixed-point reasoning are not new algorithmic principles. |
| [Dinh et al., Algorithmic approaches to avoiding bad local minima in nonconvex inconsistent feasibility](https://arxiv.org/html/2502.19052v1) | Product-space relaxed Douglas–Rachford can escape bad minima in simulated and real orbital tomography. Section 4.1 reports nonvanishing iterate changes with convergent gaps. | Actual imaging application differs from our toy benchmark. Relevant methods, fixed-point characterizations and numerical sections inspected; figure captions, not all pixels. No inference that the reported behavior is our strange attractor. |
| [Chen, Local saddle structure in relaxed averaged alternating reflections Algorithms on phase retrieval (2021)](https://arxiv.org/abs/2101.11145) | ADMM and primal-dual interpretations; parameter-dependent basins and stationary-point screening. | Primary abstract only. These equivalent terms expand the necessary future search. |
| [Elser, Solving Problems with Projections, Origins chapter (2025)](https://doi.org/10.1017/9781009475518.003) | Established context of nonconvex projection dynamics. | Publisher chapter excerpt only; full book not reviewed. |

## What the present packet supports

1. Two visually selected noncoplanar-circle regimes, with actual-state density portraits and rigid rotations, survive long fresh-start, high-precision and original-library checks.
2. A third setting has two positive mean growth estimates, with an explicitly retained negative second-rate short window in the original-library run.
3. Complete fixed-point classifications and circle-disjointness certificates at six exact settings; one pair of input circles is certified linked.
4. All 2,051 screened settings, all 53 flagged follow-ups, periodic neighbors, 771 parameter-slice settings and the full declared initial-state experiments are retained.

These are concrete case-study results. Their combination may merit further study, but it does not establish priority. The interval certificates prove finite-dimensional geometric statements; the numerical evidence remains finite-time evidence.

## Equivalent descriptions and unresolved checks

Search terms include RAAR, relaxed Douglas–Rachford, averaged alternating reflections, projection-reflection dynamics, difference maps, circle feasibility, spatial/noncoplanar/linked circles, strange attractors and hyperchaos. The β=1 relation to AAR/HPR and a difference-map instance is explicit in Luke's original paper; the relaxed family should not be silently identified with that same difference-map formula. Chen's ADMM/primal-dual interpretation concerns its phase-retrieval setting and is not automatically an explicit conjugacy for this toy map.

Translation, rigid rotation and a common spatial scale conjugate geometrically equivalent circle instances. We normalize A to the unit xy circle and orient the B-plane tilt; relabeling A/B changes the ordered iteration and needs its own check. The selected input geometries being noncoplanar does not rule out a nonlinear conjugacy or a known symbolic dynamics on their attracting sets. A global topological/symbolic comparison has not been completed.

The open requirements for a stronger claim are: a full review of older projection-chaos papers and supplements, an explicit inequivalence argument or useful invariant, a certified invariant region avoiding projection singularities, and an asymptotic chaos/attraction proof. We do not replace these with a lack of search hits.

Source hashes, search terms and inspection scopes are frozen in [the priority packet](../results/raar-circles/priority/). Third-party publication PDFs remain local and are not redistributed. The pinned [official ProxToolbox evaluator](https://gitlab.gwdg.de/nam/ProxPython/-/blob/702154fc33a0b05747ef7f2980c88613236e3626/proxtoolbox/algorithms/RAAR.py) is retained with its repository license and attribution.
