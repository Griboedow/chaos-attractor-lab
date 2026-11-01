# Local equilibrium certificates

These files extend the frozen unexpected-search packet. For each of tetrahedral
cases 421 and 565, 33 distinct equilibrium boxes pass outward interval Krawczyk
inclusion at 60 and 90 decimal digits. Routh first-column signs classify their
unstable eigenvalue counts. The root search is not exhaustive, and this is not a
chaos or historical-priority certificate.

Only the corrected certificates using the actual interval displacement X-x0 are
included. An earlier scratch attempt used a nominal displacement radius; it was
replaced before publication. The portable script additionally replaces the
floating-point box-separation check with an outward interval comparison.

Input parameters are the exact binary64 values stored in the earlier packet,
converted to intervals by integer ratios. Reports retain full interval bounds.
See ../../docs/tetrahedral-equilibria.md for scope and reproduction.
