# All 18 certified one-period saddles meet the three basins

Every neighborhood of each [one-period saddle](../cubic-forced-orbits/) contains starts converging to each of the three sinks. For the 16 moving saddles, 48 interval integrations send three portions of each unstable manifold into three different trapping disks. The two [stationary saddles](../cubic-saddle-basins/) supply the remaining cases. The inclination lemma then puts each saddle's entire stable manifold on the common boundary. [Exit certificates, phases and starting offsets](boundary-orbits.json)

The system is the [forced cubic gradient flow](../shape_search_20261121/frontiers/cubic/#a-single-gradient-field) at gain 0.34. The smallest starting displacements are around $`10^{-35}`$. A logarithmic scan in both unstable directions selected offsets reaching different sinks. The checker evolves the displacement relative to the certified moving orbit until ordinary coordinates resolve it. Each endpoint lies within 0.00181 of its phase-zero trap center, inside the certified radius of 0.0055. Repeating the weakest-expanding saddle and orbit 14 at 2,560 cells and Taylor degree 12 confirmed all six exits.

<details>
<summary>Why the starting disks meet the unstable manifold</summary>

Center the return map $`P`$ on the periodic orbit. Validated short flow segments give its expanding eigenpair $`(\lambda,v)`$, a Lipschitz bound $`H`$, and a quadratic error bound $`|P(w)-DP(0)w|\le C|w|^2`$. The verifier integrates the Hessian remainder along the enclosed variational flow and checks $`H<|\lambda|^2`$. For each signed offset $`s`$, iterating $`s\lambda^{-n}v`$ for $`n`$ periods has a limit within $`C s^2/(|\lambda|^2-H)`$ of $`sv`$. Its backward iterates tend to the saddle. The entire resulting disk is propagated with an interval Taylor defect and a Euclidean logarithmic-norm bound. [Local-limit argument](../cubic-saddle-basins/stable-crossing.md#building-the-arc-in-reverse-time) · [Checker](../../research/shape_search_20261121/frontiers/cubic/boundary_orbits.py)

</details>

```sh
python research/shape_search_20261121/frontiers/cubic/boundary_orbits.py --workers 6
```
