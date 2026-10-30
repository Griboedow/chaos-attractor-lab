# Destinations of all four physical perturbations of lane 49

All four declared perturbations (cases 16–19) were continued for 500 original tau in the complete, unprojected 23-dimensional system. Upward X0 sections and close-return lags 1–16 were retained. A fixed close-return threshold selected the smallest lag for periodic shooting. Grouping at 1e-10 was applied only after full-space sampling; all projection errors and shooting outcomes remain recorded.

| Case | Original perturbation | Sampled endpoint groups | Shooting period | Expanding / neutral multipliers |
|---|---|---|---:|---:|
| 16 | +1e-7 | 2+2+1+2+1 | 25.001684375896957 | 0 / 1 |
| 17 | -1e-7 | 2+1+1+2+1+1 | 15.657983546061825 | 3 / 1 |
| 18 | +1e-9 | 2+1+1+2+1+1 | 15.657983546061750 | 3 / 1 |
| 19 | -1e-9 | 1+2+1+2+1+1 | 15.657983546061738 | 3 / 1 |

These are floating-point conclusions. Cases 17–19 shoot back to four-pair roots: their temporarily unresolved row differences should not be interpreted as six distinct periodic clusters. Phase and optimal row-permutation fits identify them with the previously found lane49 periodic saddle to full-row RMS errors below 1.7e-13. Their three expanding multipliers are near -1.7278395673. A near-zero final growth window therefore did not establish attraction.

Case 16 instead resolves a distinct five-cluster periodic candidate with groups `[[0,3],[1,2],[4],[5,6],[7]]`. Its complete DOP853 monodromy has one flow multiplier and twenty-two multipliers strictly within the unit circle numerically; the largest nonflow modulus is 0.127482288890853. The next conjugate pair has modulus 0.075980015243. Full-state closure is below 8.2e-14.

An independent full 23D Radau state integration and all 23 Radau variational columns confirm this spectrum. Matrix disagreement with DOP853 is below 1.73e-11, matched eigenvalue disagreement below 7.41e-13, state closure below 5.59e-13, and flow-identity residual below 2.95e-12. Radau gives maximum nonflow modulus 0.127482288890874 and exactly one neutral direction. The sparse column-major Jacobian convention is directly checked. The Radau state interpolant and variational calculation are numerical, not interval enclosures.

The five-cluster candidate has numerical weighted permutation closures at T/2, T/3 and T/6. In the stored group ordering, the T/6 permutation is `[3,0,4,1,2]`, with order six, and the maximum closure residual is below 1.3e-13. The T/3 permutation `[1,3,2,0,4]` cycles the three paired clusters and fixes both singletons. These are numerical clues for a separate validated shooting calculation, not exact-symmetry or least-period proofs.

Files `protocol.json`, `summary.json`, every `case-*/analysis.json`, all full-state/section samples and the complete periodic analyses preserve the original source hashes, settings, selection and measurements. `confirmation/` contains the independently computed matrices and refined phase/permutation comparisons. `five-cluster-symmetry.json` records all tested time divisors, including failures. Executed sources are copied beside these outputs.

The cases were selected for this follow-up after observing their near-zero final GPU windows; this is not random basin sampling and does not change the original cohort. The evidence supports a distinct strongly attracting periodic candidate for case 16 and return to symmetry-related saddles for the other three cases. It does not prove interval stability, asymptotic basin fractions, historical novelty, or a chaotic attractor. A separate interval theorem, if completed, belongs to a separate proof package.
