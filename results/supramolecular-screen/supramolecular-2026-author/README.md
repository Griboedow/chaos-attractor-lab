# Reaction-diffusion simulation of coupled cooperative supramolecular polymer system

The code is written to run numerical integration of a two-dimensional reaction-diffusion system using Euler time integration and finite-difference five point stencil for spatial discretization. 

---

## Model variables
d    : deactivated monomer concentration </br>
a1   : activated monomer concentration </br>
m1   : polymer mass concentration </br>
m0   : polymer number concentration </br>

## Model parameters
ka   :  monomer activation rate constant </br>
kd2  :  polymer end deactivation rate constant </br>
kd3  :  polymer chain deactivation rate constant </br>
kp   :  rate constant for polymer assembly </br>
ctot :  total monomer concentration 

## Numerical method
- Explicit fourth-order Runge-Kutta time scheme
- Second-order central differences for diffusion
- No-flux (Neumann) boundary condition

## Initial condition
We fix the initial condition depending on the type of perturbation we apply on the system.
- Random perturbation
- Gaussian perturbation

## Reference
If you find the code useful please cite:</br>
```
@misc{singh2026emergencespatiotemporalpatternsfueldriven,
      title={Emergence of spatiotemporal patterns in a fuel-driven coupled cooperative supramolecular system}, 
      author={Akta Singh and Nayana Mukherjee and Jagannath Mondal and Pushpita Ghosh},
      year={2026},
      eprint={2601.15662},
      archivePrefix={arXiv},
      primaryClass={cond-mat.soft},
      url={https://arxiv.org/abs/2601.15662}, 
}
```
