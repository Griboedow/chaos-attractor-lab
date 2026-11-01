# Ice-cloud microphysics: no retained chaotic candidate

All 2,048 trajectories in the declared coefficient-family screen complete. None passes the primary chaos filter. The two weaker signals, and their companion starts, close to numerical period-one orbits under an independent implicit solver. This branch is therefore deprioritized for the search for an unusual strange attractor.

## Model and scope

We implement the full three-state positive-interior model in Eqs. 55–57 of [Bergner and Spichtinger, *Ice clouds as nonlinear oscillators*, arXiv:2507.03475v1](https://arxiv.org/html/2507.03475v1), subsequently published in [Chaos 36, 043115 (2026)](https://doi.org/10.1063/5.0297531). Its states are crystal number mixing ratio, ice mass mixing ratio, and saturation over ice. Nucleation, vapor growth/evaporation and sedimentation are retained, including the small nucleation-mass term omitted in the source's later analytical reduction.

This is the published autonomous approximation with fixed pressure and temperature while retaining an updraft source. It does not simulate the changing temperature and pressure of an actually ascending parcel. We declare `cp=1005 J/(kg K)` because that value was not recovered from the source. The final printed sedimentation coefficients are used literally; an earlier air-density correction is not silently inserted. These choices prevent a claim of exact reproduction of a measured cloud or a particular source figure.

The source already establishes periodic cloud oscillations and maps their environmental parameter dependence. The exploratory extension here varies aerosol abundance and the mass-distribution moment ratio as well. No state, forcing term or feedback law is added. [Executed model](../research/ice_cloud.py).

## Source controls and a preserved failed check

At 230 K, 300 hPa, sedimentation factor 1, aerosol abundance 300 per cubic centimetre and moment ratio 3, independent 100-unit integrations approach stable equilibria at updrafts .01 and 1 m/s. At .1 m/s they approach a cycle with period 1,823.851 seconds. The maximum difference across its final 100 full-state section returns is `1.68e-9`.

A separately coded dimensional field integrated with Radau agrees with normalized logarithmic-coordinate DOP853 within scaled error `2.57e-10` on the short controls. The final five-point Jacobian assertion then fails at extreme exponential-dominated random states: subtraction loses the much smaller derivative components. This failure is preserved alongside all three completed histories. A separate audit replays those same 48 smooth-side states using complex-step differentiation of the independent dimensional field. Maximum scaled Jacobian disagreement is `5.36e-15`; field disagreement is `2.67e-14`. No differentiability at the evaporation switch is asserted. [Original failure](../results/ice-cloud-screen/ice-cloud-controls/failure.json) · [Replacement derivative audit](../results/ice-cloud-screen/ice-cloud-algebra/validation.json).

Before searching, five physical tuples at two tolerances pass independent Radau/GPU state-and-tangent controls. The largest discrepancy among these short state, tangent and growth-rate comparisons is below `1.2e-8`. [GPU controls](../results/ice-cloud-screen/ice-cloud-gpu-controls/validation.json).

## Search and independent returns

The design uses 1,024 scrambled Sobol tuples, each with two declared starts: a 1% perturbation of the full equilibrium, and a separate low-number/high-mass state. Ranges are 190–240 K, 200–300 hPa, updraft .0005–2 m/s, sedimentation factor .01–1, aerosol abundance 3–3,000 per cubic centimetre, and moment ratio 1.05–30. The broad width range is coefficient exploration, not a validated family of atmospheric size distributions.

Time is rescaled by the local equilibrium's linear period; this is a fixed positive scaling for each tuple. After 200 warm-up units, each lane is observed for 800 units. The primary rule requires all four 200-unit leading estimates above .005; the weaker rule uses .001. There are **zero primary and two weak** cases. [Complete design and conventions](../results/ice-cloud-screen/ice-cloud-pilot/protocol.json) · [Screen outcome](../results/ice-cloud-screen/ice-cloud-pilot/analysis.json).

Independent dimensional-field Radau continues both starts at each weak tuple for 400 rescaled units. The declared return rule uses the final 50 consecutive upward equilibrium-saturation crossings, comparing all three logarithmic state coordinates with tolerance `1e-7`.

| Cases | Crossings per run | Largest final return difference | Period in seconds |
| --- | --- | --- | --- |
| 1872 / 1873 | 148 / 147 | `6.03e-12` | 4,032.731 |
| 1930 / 1931 | 92 / 93 | `2.23e-11` | 195,355.519 |

Both pairs meet the numerical period-one rule. The weak positive estimates therefore do not survive the return check. The long physical period in the second pair is a result of the frozen-environment approximation, not evidence that a real parcel maintains those conditions for days. [All four histories and return results](../results/ice-cloud-screen/ice-cloud-returns/analysis.json).

The result excludes these weak signals from our candidate list. It is not a proof that chaos is absent from the whole model family. No new cloud oscillator or historical priority is claimed. All saved runs, original and corrected controls, inputs, statuses and source snapshots are retained in the [manifest](../results/ice-cloud-screen/manifest.json); the [scope record](../results/ice-cloud-screen/scope-and-corrections.json) identifies the exact archive contents.
