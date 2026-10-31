# Soft pneumatic limbs: a source-code gate before a chaos search

Comoretto, Schomaker and Overvelde's [2025 Science paper](https://doi.org/10.1126/science.adr3661) describes autonomous soft limbs driven by a constant air supply. Its [accepted manuscript and supplement](https://ir.amolf.nl/pub/11152/17082acceptedVersion.pdf) provide a discrete mechanical model; [the authors' version-2 data deposit](https://zenodo.org/records/15399470) contains its code. This is an interesting physical setting, but the retrieved implementation is not yet a suitable basis for claiming a new attractor.

The model has 17 planar nodes, elastic forces and delayed stiffening associated with a moving kink. Its discrete memory is essential to the implemented dynamics. We read the model notebooks without executing them and probed the unchanged `helper.py`. Its SHA-256 is `c021a5b02b12b290f5320fdd705e2b98c61fb33d0e1be07e2fdc48a3c2930862`.

## What the reproducible probes show

Our [audit script](../research/audit_pneumatic_author_code.py) finds four specific problems relevant to numerical dynamical classification:

1. The differential-equation callback modifies its input state in place. In the declared 17-node probe, the largest modification is 4 coordinate units.
2. It also changes an external kink-history array during derivative evaluation. Repeating the same time and explicit state after this update changes the derivative by `5.55e-10` in the probe. The value is small in this example, but the callback is history dependent and does not roll back rejected solver stages. That is a specification problem, not an estimated chaos effect.
3. The torsional-force routine overwrites its sigmoid-midpoint argument `x0` with a neighboring node's horizontal coordinate. Changing the intended midpoint from 1.5 to 100 has zero effect; translating the same shape horizontally by 10 changes the computed torsional force by as much as `0.1823`.
4. At positive PCHIP knots `0, 0.3, 2π/3, 2π`, supplement Eq. 16 lists values `0, 1.2, 0.1, 0.2`, whereas the notebook supplies `0, 0.2, 0.1, 1.2`. The supplement was checked visually on PDF page 38.

The [frozen audit record](../results/pneumatic-source-audit/) contains the exact probe output, source-file hashes and archive paths. Original author code is not redistributed in this packet. Public HTTP byte-range requests retrieved only the relevant ZIP entries; ZIP CRC checks verified extraction.

## Research decision

Do not interpret a positive numerical exponent from the unchanged helper as a property of a well-defined autonomous ODE. A follow-up would first need an explicit hybrid state, an accepted-event rule for the memory, a resolved force convention, and convergence checks under changes of solver and step size. Any corrected variant must be identified as such.

These findings concern this retrieved implementation and its suitability for our proposed chaos search. They do not refute the experimental oscillations, establish that the physical device cannot be chaotic, or settle whether its dynamics has prior equivalents. No chaotic attractor is claimed or simulated here.
