# Icosahedral animation inputs

`icosa-animation.npz` stores 320,000 three-dimensional states for each of cases
177 and 183 (`points_177`, `points_183`). These are every fourth state of the
original 1,280,000-state dense trajectories used for the 16,000-unit portraits.
This packet adds a denser rendering sample; it does not contain a new simulation.

The dense source hashes, dimensions and selection are recorded in
[`source.json`](source.json). Before selection, every sixteenth dense state was
checked bitwise against the corresponding trajectory in the frozen
[`icosa-long/sample.npz`](../icosahedral-flow/icosa-long/sample.npz). Equivalently,
every fourth state of each new array equals that frozen 80,000-state sample.
`python tools/check_portraits.py` repeats the latter check.

[`tools/animate_portraits.py`](../../tools/animate_portraits.py) combines these
states with the archived initial camera and rotates the view at a fixed scale.
The animation record includes the input hash, selected count and every frame's
retained point count. The scientific equations, integration settings and
validation remain in the [original packet](../icosahedral-flow/README.md).
