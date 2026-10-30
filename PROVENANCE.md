# Provenance

## Public-clone verification

On 30 October 2026, an independent public clone of commit `a5145092d86ff5ad51d42bc5e1a7a4eb64d2d802` was installed in a new Python environment. Four GPU/CPU checks and the short optional DCR control passed, as did `pip check`. Execution used isolated Python outside both repositories and left the clone's tracked files unchanged. [Numerical packaging record](provenance/public-clone-validation.json) and [718-file comparison](provenance/publication-copy-integrity.json) record the exact scope. The OpenCL binding emitted a nonempty-compiler-output warning during its first build; all requested controls completed successfully. This pass did not inspect that build output, and it does not certify every archived model.

The reusable [verification driver](tools/verify_public_clone.py) requires an explicit clean clone and separate output directory. Historical staging checks below retain their actual earlier snapshots; they have not been relabeled as the public-clone run. The `prepublication-manifest.json` is the immutable staging inventory before publication wording changes.

## Source extraction

Subsequent exact copies are recorded individually under `provenance/increments/`: the optical cycle fold from `f0010d1`, molecular-motor core audit from `643d8eb`, N-row source controls from `1584ba3`, and optical saddle departures from `ffd3892`. These records preserve source revisions and file hashes; copying did not rerun the experiments. The motor GPU surveys are still separate work in progress.

This project was extracted from [chua-gpu-atlas](https://github.com/Griboedow/chua-gpu-atlas) at commit `b60418fc2fbdd8a34ab873c5581e31f1b9c0727a` for publication as [chaos-attractor-lab](https://github.com/Griboedow/chaos-attractor-lab). The original MIT license and Copyright (c) 2026 Griboedow are retained verbatim. The old repository's history is not rewritten.

`provenance/source-manifest.json` records original Git blob IDs and SHA-256 hashes. Four research modules replace an inherited Chua runtime with a generic OpenCL device/context/buffer helper; their vector fields and paired kernels are unchanged. No numerical experiment has been rerun as part of this preparation. Historical repository commit URLs in the notebooks are intentionally retained. The historical requirements lock includes dependencies used by the Chua project and is evidence of the original environment, not an independently verified lock for this project.

All findings retain their source-model distinctions, physical-domain warnings, finite-time limitations and negative outcomes. Referenced papers, cached PDFs, source archives and downloaded browser captures are not redistributed here or relicensed under MIT. The microfluidic source controls, stiff-gel endpoint follow-up, mineral-protocell audit, DCR source controls, optical-order/priority audits, optical-cycle continuation, optical lifetimes through 100,000 thermal-time units, and heat-pipe/thermal-flutter source audits are included at this snapshot. Later science still requires a final provenance refresh. All Chua parameter-scaling work remains in the Chua repository. Historical rendering-audit links to excluded Chua reports are changed only in this staged copy to their immutable original-repository URLs.

The sanitized `provenance/packaging-validation.json` preserves the parent-run fresh-install, pip-check and four GPU/CPU controls at `8a3564c`, with a byte-identity audit for the same tested source files at the refreshed snapshot. Those computations were not rerun. `requirements-lock.txt` records the fresh Windows environment including the separately installed optional DCR dependency; its provenance is in `provenance/optional-dcr-packaging.json`. The older original-project lock remains separately labeled as history.

The optional DCR record preserves a separate parent-run isolated CPU control at `3cf58d0`, including its original timestamp and duration 0.01. The model, optional requirements and archived comparison array are byte-identical at this refreshed snapshot. Absolute workstation paths were removed from the public provenance copy. Neither that computation nor the older four-GPU controls were rerun.

`docs/research-findings.md` is a new synthesis prepared for this standalone project, not a copied document from the source Git tree. Its evidence comes from the individually cited reports. The source manifest labels it as derived and records its content hash.
