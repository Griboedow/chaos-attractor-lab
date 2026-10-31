# Minimal biofilm oscillator: an exact prior-art exclusion

The three-variable metabolic model of Garde, Ibrahim, Kovacs and Schuster (2020) is already covered by the global analysis of the Wilhelm–Heinrich chemical oscillator. It is not a new autonomous-chaos target. This decision follows from an exact equation match, not a failed numerical search.

With constant external glutamate `GE`, the biofilm equations are

```math
\dot G_p=(k_1G_E-k_4)G_p-k_2AG_p,\qquad
\dot A=k_5G_i-k_3A,\qquad
\dot G_i=k_4G_p-k_5G_i.
```

These coincide with Smith's equation (1.1) under `(x,y,z)=(Gp,A,Gi)` and `k=k1*GE-k4`. Alternatively, the invertible positive scaling

```math
x=\frac{k_2k_4}{k_3}G_p,\qquad y=k_2A,\qquad
z=\frac{k_2k_5}{k_3}G_i
```

gives his equation (1.2):

```math
\dot x=x(k-y),\qquad \dot y=k_3(z-y),\qquad
\dot z=k_5(x-z).
```

Smith identifies this as a competitive and monotone cyclic feedback system and proves the relevant global convergence results. For `k>k3+k5`, positive solutions outside the equilibrium's stable manifold approach periodic orbits. His no-chaos conclusion therefore applies to the literal minimal biofilm model. The biochemical interpretation does not change this equivalence.

This is a bounded exclusion: it does **not** automatically settle reversible reactions, Michaelis–Menten uptake, six-compartment extensions, delay, spatial transport, or time-dependent forcing. Each changes the equation-level question. No numerical scan was performed for the excluded base model.

## Primary sources checked

1. Garde, R., Ibrahim, B., Kovacs, A. T., and Schuster, S. (2020). *Differential equation-based minimal model describing metabolic oscillations in Bacillus subtilis biofilms*. Royal Society Open Science 7, 190810. [DOI](https://doi.org/10.1098/rsos.190810), [institutional full text](https://backend.orbit.dtu.dk/ws/files/205732981/75_Garde_RoyalSocOpenSci_2020.pdf). Equations and rate table checked.
2. Smith, H. L. (2011). *Global Dynamics of the smallest Chemical Reaction System with Hopf Bifurcation*. [Primary full text](https://arxiv.org/html/1109.5204v1), especially equations (1.1)–(1.2) and Theorem 1.1. This explicitly addresses the earlier open possibility of chaos.
3. Garde, R., Ibrahim, B., and Schuster, S. (2020). *Extending the minimal model of metabolic oscillations in Bacillus subtilis biofilms*. Scientific Reports 10. [DOI / full text](https://doi.org/10.1038/s41598-020-62526-6). Its extensions are distinct targets; the base-model theorem is not silently transferred to them.

Audit date: 31 October 2026. The new contribution here is only the screening cross-reference, not the theorem or the oscillator.
