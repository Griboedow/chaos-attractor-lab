# FastICA evidence packet

[Report](../../docs/fastica-dynamics.md) · [Reproduction](../../research/fastica/README.md)

The packet records a 576-case screen, five 200,000-step follow-ups, five
five-million-step tangent spectra, 320 fresh starts, 96 input perturbations,
library execution, stopping-rule continuations and independent precision
controls. Case 40 is rejected as a long transient. Cases 15 and 351 have
regular competing outcomes. No historical priority or asymptotic chaos
certificate is claimed.

| Files | Scope |
| --- | --- |
| `pilot-01/` | Declared design, all cases and all scout outcomes |
| `follow-01/` | Five exact input matrices and initial vectors; intermediate spectra |
| `long-01/report.json` | All long spectra, including the captured case |
| `basin-01/` | 320 exact starts and all finite outcomes |
| `perturb-01/` | Perturbation seeds, hashes and all 96 outcomes |
| `precision-70.json`, `precision-100.json` | Independent mpmath growth windows |
| `library-controls-v2.json` | Library receipt with explicit full/partial window counts |
| `controls.json`, `stopping.json` | Step/Jacobian checks, stable capture and stopping-state continuations |
| `histograms/` | Unthinned portrait counts and exact cameras |
| `source-original/`, `portability.json` | Executed scripts and path-only reproduction mapping |
| `external-data.json` | Hashes of 149,783,441 bytes of generated arrays kept on D: |
| `manifest.json` | Hashes of the committed packet and reproduction sources |

The earlier `library-controls.json` receipt retains the initial observation
run; its partial-window accumulations must not be read as full-window rates.
The second run records the stopping weight and explicit counts. Both runs
agree on the first-unit exit of Case 246.

The long and intermediate runs save temporally decimated clouds only in the
external run directory. Those clouds are not portrait inputs. The separate
portrait run uses consecutive states and no sign folding. No PDF or large
trajectory array is needed to read or audit this packet.
