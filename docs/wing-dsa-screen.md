# Delayed stretch activation: a finite negative screen

The positive-feedback four-state wing model printed by Lynch et al. (2022) yielded no retained chaos candidate in the declared survey. Ten weak signals also failed retention after 40,000 additional observed time units. This result does not classify every regime of the model, and does not cover the opposite feedback sign in later insect-flight formulations.

## Source and implementation

The source is Lynch, Gau, Sponberg and Gravish, *Autonomous Actuation of Flapping Wing Robots Inspired by Asynchronous Insect Muscle*, ICRA (2022), [author-hosted full paper](https://sponberg.gatech.edu/wp-content/uploads/2022/08/Gau-et-al-2022_02.pdf). We implement its literal dimensional Eq. (9), with the coefficients of Eq. (4):

\[
I\ddot\theta+k\theta+\Gamma|\dot\theta|\dot\theta-\mu f=0,
\]
\[
\ddot f+r_3(1+\kappa)\dot f+\kappa r_3^2 f
=r_3(1-\kappa)\dot\theta.
\]

All physical coefficients are positive and \(0<\kappa=r_4/r_3<1\). The printed first-order Eq. (10) omits a division of drag by inertia relative to Eq. (9); the implementation follows Eq. (9) explicitly. The prose stiffness 0.416 is used, not the inconsistent figure-caption value. The quadratic-drag value needed for quantitative angular amplitudes was not recovered.

With \(\omega_n=\sqrt{k/I}\), \(\tau=\omega_nt\), \(\theta_*=I/\Gamma\), \(x=\theta/\theta_*\) and \(h=\mu f/(k\theta_*)\), the exact normalized equations are

\[
x'=v,\quad v'=-x-v|v|+h,\quad h'=w,
\]
\[
w'=\eta r(1-\kappa)v-\kappa r^2h-r(1+\kappa)w,
\quad r=r_3/\omega_n,\quad\eta=\mu/k.
\]

Quadratic drag sets the coordinate unit. Choosing a drag value for an algebra test does not fit an experimental amplitude. Three 1,600-unit controls compare DOP853 against an independent transformation of the dimensional equation integrated with Radau. Maximum scaled differences are \(2.60\times10^{-11}\), \(1.12\times10^{-9}\) and \(4.34\times10^{-10}\). The cases include a stable origin, a large oscillation and the reported mechanical/gain tuple. Random field, Jacobian and characteristic-polynomial checks also pass.

## Declared search and outcome

The 8,192-proposal Sobol design uses logarithmic bounds \(r\in[0.03,30]\), \(\kappa\in[0.01,0.99]\), \(\eta\in[0.001,100]\), seed 2026110103. Selecting an origin growth rate above 0.001 leaves 4,098 tuples. Two distinct fresh starts produce 8,196 trajectories. Hidden attractors around a stable origin are outside this design.

| Stage | Complete | Other outcome | Retained |
| --- | ---: | --- | ---: |
| Warm 500, observe 2,000; four windows above 0.005 | 7,530 | 666 magnitude-cutoff exits | 0 |
| Same pilot, exploratory four-window threshold 0.001 | 7,530 | Same exits, not another survey | 10 weak |
| All ten weak endpoints: warm 2,000, observe 40,000; four windows above 0.001 | 10 | No failure or budget exit | 0 |

The 666 cutoff exits are unresolved, not evidence of either bounded chaos or asymptotic divergence. The cutoff is \(10^7\) in normalized coordinates. No global absorbing-set proof has been supplied for the whole parameter box.

Independent Radau continuations of all ten weak pilot endpoints were also attempted over 800 units after a 500-unit warm-up. They have only one to six upward zero-displacement crossings each, too few for the declared return criterion. They must not be described as verified cycles. The subsequent longer GPU observation lowers every case below the retention rule, but does not settle arbitrarily weak asymptotic growth. Some normalized amplitudes are very large; these are mathematical parameter explorations, not predictions for feasible insect wing angles.

## Equivalent terminology and remaining scope

Delayed stretch activation is already an established self-excitation mechanism. The [2023 Nature paper by Gau et al.](https://sponberg.gatech.edu/wp-content/uploads/2024/06/s41586-023-06606-3.pdf) expresses it as a convolution with a difference of exponentials; its Eqs. (10–11) include a **negative** strain-rate feedback sign. The [supplement](https://sponberg.gatech.edu/wp-content/uploads/2024/12/41586_2023_6606_MOESM1_ESM.pdf) derives the normalized kernel and work output. Their reported simulation uses a sampled, truncated FIR filter. Its continuous, untruncated kernel has a finite-state realization, but that realization is not identical to the sampled implementation.

The [2025 Wold et al. preprint](https://sponberg.gatech.edu/wp-content/uploads/2025/06/2025.05.07.652680.full_.pdf), *Supra-resonant wingbeats in insects*, also prints negative velocity feedback in its force equation. Its complete supplementary coefficient specification has not been read. Accordingly, this positive-sign 2022 screen cannot be relabeled as an exclusion for the later model, nor does a different notation establish a new mechanism. The 2026 *Supra-resonant wing oscillations in insect flight* article is a separate paper and was not substituted for this preprint.

The [frozen packet](../results/wing-dsa-screen/manifest.json) records the numerical settings, exact executed code, all pilot endpoints, all long endpoints and the independent return events. Compact histories preserve their final 50 units; complete local histories are hash-indexed. The appropriate current decision is **no retained candidate in this declared positive-sign search**, with physical calibration and the later opposite-sign formulation explicitly outside its conclusion.
