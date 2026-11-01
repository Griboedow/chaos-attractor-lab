# Certified symbolic dynamics in spatial-circle RAAR

The eight edges of a branching covering graph pass two interval map implementations. Standard covering theory gives a compact invariant set with topological entropy at least log(2)/49 per original iteration. This is not a certificate of an attracting chaotic set or historical uniqueness.

All runs are preserved under their original directory names in `runs`; source paths inside reports record provenance. The primary certificates are `raar-covering-proof-pilot01` and `raar-covering-proof-edge1-01` through `edge7-01`. The second implementation is `raar-covering-ad-edge0-02` through `edge7-02`. Both use the identical exact input from `raar-covering-refined01/input.json`.

The initial 36 numerical window proposals all failed their nonlinear sample tests. The refined proposal succeeds numerically and is then certified. The interrupted AD pilot `raar-covering-ad-edge0-01` is incomplete; its gzip stream may lack a footer and is never counted as a certificate. The partial structural audit is retained separately from both complete audits.

The two map implementations use analytic vector projectors versus scalar nested reflectors with automatic derivatives, and different matrix inversion methods. They share the subdivision driver and mpmath interval library. These are computer-assisted interval proofs, not proof-assistant formalizations.

The located least-period-49 saddle has a separate 200-digit existence/least-period certificate, a certified window itinerary, and exactly one unstable multiplier. Longer coded multiple-shooting examples are numerical illustrations. Their binary64 monodromy products are ill-conditioned: small eigenvalues must not be used as stability evidence.

The return scout's legacy `max_excursion` field measures distance from one selected phase of the seven-cycle, not distance from the whole cycle. It is not used as a geometric excursion estimate. The final itinerary figure uses nearest distance to all seven cycle points.

Final figure files are in `runs/raar-symbolic-render02`. Earlier renderings are retained. Run `python research/verify_raar_symbolic.py` from the repository root. See [the report](../../docs/raar-symbolic-dynamics.md) and [manuscript](../../papers/raar-symbolic-dynamics/manuscript.pdf).
