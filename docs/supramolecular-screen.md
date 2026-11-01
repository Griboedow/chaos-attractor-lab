# Cooperative polymerization: no retained chaos candidate

The published approximate polymer model produced strong short-run signals in ordinary concentration coordinates, but none survived the completed long follow-ups in positive logarithmic coordinates. This is a useful numerical exclusion from the current search queue, not a theorem that the model cannot be chaotic.

The [frozen evidence packet](../results/supramolecular-screen/manifest.json) includes the declared inputs, executed source, failed-domain outcomes, corrected calculations, independent return checks and hashes of the full local histories. No parameter tuple from this screen is promoted as a new attractor.

## Source and prior work

The equations come from Sharko, Livitz, De Piccoli, Bishop and Hermans, *Insights into Chemically Fueled Supramolecular Polymers*, Chemical Reviews 122 (2022), 11759–11777, [DOI 10.1021/acs.chemrev.1c00958](https://doi.org/10.1021/acs.chemrev.1c00958). The public [author dissertation](https://publication-theses.unistra.fr/public/theses_doctorat/2022/SHARKO_Anastasiia_2022_ED222.pdf), printed pages 48–49, supplies the approximate moment equations, correction coefficient and periodic reference. These are a long-polymer closure, not the full polymer-size reaction network.

The same homogeneous closure reappears in Singh et al., [arXiv:2601.15662v1](https://arxiv.org/abs/2601.15662), and the later article *Reaction–Transport Coupling Drives Spatiotemporal Organization in Fuel-Driven Supramolecular Polymerization*, [DOI 10.1002/anie.5055917](https://doi.org/10.1002/anie.5055917). The [author code](https://github.com/JMLab-tifrh/supramolecular_polymer_patterns/tree/cf2a8f464350b14f56ddc5fbe05f7de5c7304406) is pinned and its Git blob hashes checked. It was inspected, not executed. Its spatial transport model is a separate system. The final article's full equations were not recovered; the equation comparison uses the full preprint and pinned code.

Thus neither this reaction mechanism nor its regular oscillation is a new 2026 discovery. A future chaos claim would need to distinguish the homogeneous closure from spatial patterns and to address the closure's physical limits.

## Exact reduced equations used

With the source choice of elongation-rate unit, write

\[
x=a/k_{d,2},\quad y=d/k_{d,2},\quad z=n/k_{d,3},\quad
\tau=k_{d,3}t,
\]

and parameters \(\alpha=k_a/k_{d,3}\), \(\beta=k_{d,2}/k_{d,3}\), \(T=c_{\rm tot}/k_{d,2}\), \(C=4.6\). Conservation removes polymer mass, giving

\[
x'=\alpha y-2xz,\qquad y'=-\alpha y+2z,
\]
\[
z'=\beta(T-x-y)-6z-z^2
-\frac{2\beta C z^3}{[\beta(T-x-y)-2z]^2}.
\]

The physical region has positive concentrations and mean polymer length \(L=\beta(T-x-y)/z>2\); the long-polymer approximation needs more than this boundary inequality. The source reference is \((\alpha,\beta,T,C)=(0.2,1000,10,4.6)\).

An independently implemented four-equation source integration and the reduced three-equation integration agree to a maximum scaled difference of \(1.01\times10^{-8}\). The reproduced period is 0.746056, consistent with the source's approximately 0.75; observed mean lengths range from 14.606 to 39.379. The complex-step Jacobian error is below \(2.4\times10^{-15}\). These checks validate the implementation at the reference, not every later stiff regime.

## Declared survey and numerical correction

The initial Sobol design contains 8,192 proposals with

\[
\alpha\in[0.005,20],\quad\beta\in[10,10^4],\quad T\in[1.05,100],\quad C=4.6.
\]

Selection requires an unstable positive equilibrium with largest real eigenvalue above 0.001 and equilibrium mean length above 10. Including the reference gives 1,697 parameter tuples and two fresh mass-exchange perturbations each: 3,394 lanes. This filter does not search hidden attractors around a stable equilibrium. The observation follows 100 transient units and lasts 800 units, with four separate rate windows.

| Calculation | Complete | Unresolved interval | Signals passing all four windows |
| --- | ---: | ---: | ---: |
| Ordinary concentration coordinates, threshold 0.005 | 1,231 | 2,163 domain exits | 106 preliminary |
| Exact positive/logarithmic coordinates, threshold 0.005 | 3,394 | 0 | 114 preliminary |
| Same corrected run, exploratory threshold 0.001 | 3,394 | 0 | 260 preliminary |
| Continuation of all 260 weak signals: warm 500, observe 8,000, threshold 0.001 | 256 | 4 step-budget exits | 0 retained |

The first implementation loses positivity when polymer number is extremely small. An initial short solver agreement does not protect such a later regime. The correction changes coordinates exactly, without a population floor, clipping or a changed vector field:

\[
u=\log z,\quad v=\log(L-2),\quad E=e^v,\quad L=2+E,
\]
\[
F=E-4-2\beta C/E^2-e^u,
\]
\[
u'=F,\qquad v'=\frac{2\beta(x-1)-LF}{E},\qquad
x'=\alpha\left(T-x-\frac{e^uL}{\beta}\right)-2xe^u.
\]

Finite-time growth rates depend on coordinates; these preliminary counts must not be interpreted as invariant Lyapunov classifications. The same physical starts and parameters were used for the full corrected recheck. Initial field and tangent comparisons against an independent implicit solver pass, and every corrected pilot lane completes.

## Independent return checks

The source and the three strongest selected raw-coordinate signals (IDs 1958, 3170 and 3316) close numerically as period-one orbits under an independent Radau integration of the exact log equations. The return criterion uses all three coordinates at solver-root upward crossings of \(x=1\), not an apparent closed curve in a plot.

Four selected slow cases (415, 548, 1832 and 2750) were also observed for 3,000 units from fresh starts. Their last 100 returns do not close within the declared \(10^{-7}\) scaled tolerance. They show long amplitude/period drift, including log polymer numbers below −2,700. This is not accepted chaos evidence, and the physical long-polymer approximation becomes questionable at such tiny number concentrations. Floating-point underflow of \(e^u\) is not a physical extinction event.

All four step-budget exits from the 8,000-unit continuation (152, 1531, 3163 and 3285) were subsequently continued from their saved log endpoints with an independent implicit solver. All four close numerically with lag one. Their original explicit-integration intervals remain incomplete; the endpoint checks do not retroactively complete them.

The largest mean rate among the 256 completed long continuations is 0.002372, but its four windows decrease through approximately 0.00530, 0.00336, 0.00061 and 0.00022. No completed long trajectory passes the weak threshold in every window.

## Interpretation and reproducibility

The outcome is **no retained candidate in this finite search**. It does not exclude other parameter regions, stable-equilibrium hidden attractors, the full polymer-size network, or the spatial model. It also does not establish exact periodicity from floating-point returns.

The packet retains all section events and the final 50 time units at original sampling for each independent return check. Full histories remain locally hash-indexed. A separate selection clarification corrects the overly broad boilerplate sentence in the copied slow-follow-up protocol without rewriting that run's original record.

Reproduction entry points are [the source control](../research/check_supramolecular_source.py), [declared pilot](../research/search_supramolecular.py), [log recheck](../research/recheck_supramolecular_log.py), [long continuation](../research/extend_supramolecular_log.py) and [independent returns](../research/check_supramolecular_returns.py). Their executed versions are also frozen with each run.
