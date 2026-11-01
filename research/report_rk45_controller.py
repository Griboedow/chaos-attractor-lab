"""Write the English research report from the frozen RK45 arrays."""
import argparse,json
from pathlib import Path
import numpy as np
from verify_rk23_controller import read


def main(args):
    root=args.root;js=lambda p:json.loads(p.read_text(encoding='utf-8'))
    scout=read(root/'oscillator/scout/screen.npz');follow=read(root/'oscillator/follow/follow.npz');long=read(root/'oscillator/five-million/follow.npz')
    cohort=read(root/'oscillator/cohort/cohort.npz');norm=js(root/'oscillator/error-norm-comparison/report.json')
    out=['# Twin rosettes and folded fans in the standard RK45 solver\n',
    '**Status:** strong finite numerical candidates in an existing algorithm. Historical uniqueness, asymptotic chaotic attraction and a new invariant-set equivalence class remain unresolved. The physical oscillator is stable.\n',
    '![Fixed actual-state clouds viewed by rotating rigid cameras](assets/rk45-oscillator-rotation.gif)\n',
    '[Large portrait](assets/rk45-oscillator-geometry.png) · [Initial-volume approach](assets/rk45-oscillator-approach.png) · [Independent evidence](assets/rk45-oscillator-evidence.png) · [Actual-library geometry](assets/rk45-library-clouds.png) · [Exact map and solution bounds](rk45-controller-structure.md) · [Prior-art and equivalence audit](rk45-priority.md)\n',
    '[Research manuscript (PDF)](../papers/rk45-controller/manuscript.pdf) · [Editable text, figures and review](../papers/rk45-controller/README.md)\n',
    '## Why these two cases are interesting\n',
    'D335 forms two intricate rosettes and retains **two positive finite growth estimates** in every long, independent high-precision and actual-library window. D448 forms strongly folded fans and retains one positive estimate in every corresponding check. These descriptive names are visual labels, not claims of a new attractor class.\n',
    'The unexpected-looking object lives inside the numerical solver. We integrate only `q\'=p`, `p\'=-q-2*zeta*p`; every physical solution decays to zero. The extra state is the proposed numerical step. The method, its coefficients, safety factors, componentwise scales, RMS norm and rejection policy are the existing SciPy RK45 implementation. Earlier numerical-analysis work already studies adaptive feedback, unstable controller equilibria and tolerance-scale oscillations.\n',
    '## Exact settings and five-million-step evidence\n',
    '| Case | Damping | Relative tolerance | Mean full finite spectrum | Minimum leading window |\n| --- | ---: | ---: | --- | ---: |']
    for ident,rates in zip(long['ids'],long['rates']):
        p=scout['params'][ident];out.append(f'| D{ident} | {p[0]:.17g} | {p[1]:.17g} | '+', '.join(f'{v:.9f}' for v in rates.mean(0))+f' | {rates[:,0].min():.9f} |')
    out+=['\nEach run starts afresh at normalized `(q,p,h)=(1,0,0.1)`, discards 100,000 accepted steps and measures five million in twenty equal windows. A full three-column QR calculation uses the analytic branch Jacobian in `(q,p,log h)`. The sum of the window rates matches the independently accumulated log determinant within 9.1e-14. No selected tail closes within the declared 1–512 lag test at tolerance 1e-9. This finite test is not proof of aperiodicity.\n',
    'D508 and D726 are retained secondary cases. Their five-million-step checks do not substitute for independent precision and library checks, which were performed for D1890, D335, D448 and D2006.\n',
    '## Independent arithmetic and the actual library\n',
    'The 70/90-digit implementation evaluates all seven Dormand–Prince stages with rational method coefficients. It imports neither the production polynomial map nor its Jacobian. After 100,000 warm-up and 5,000 alignment steps, it measures 100,000 steps in ten windows using two finite-separation directions in actual `(q,p,h)` coordinates. Separations are 1e-30 and 1e-40; parameters and controller constants are exact lifts of the declared binary64 values.\n',
    'The unmodified SciPy 1.18.1 runs each take 400,000 accepted steps from physical `(q,p)=(1,0)`, `atol=1e-6`, first step 0.1. They retain the floating-point clock. After 100,000 warm-up and 10,000 tangent-alignment steps, 29 windows remain. Tangents describe the mathematical map evaluated on library states, not the derivative of floating-point rounding. Maximum one-step residuals are below 9.4e-9 in normalized solution/log-step coordinates.\n',
    'An additional unaligned comparison uses 200,000 late states per method. Median actual-library-to-reference distances are 0.00250248 (D335) and 0.00132266 (D448), versus reference half-cloud medians 0.00325826 and 0.00188043. Maximum distances remain nonzero (0.35591 and 0.18432); no identical-invariant-set or shadowing claim follows.\n',
    '![Same cameras and scales for transcription and actual SciPy](assets/rk45-library-clouds.png)\n',
    '| Case | Minimum leading rate, 70 digits | 90 digits | Actual SciPy |\n| --- | ---: | ---: | ---: |']
    for ident in [335,448,1890,2006]:
        vals=[np.array(js(root/f'oscillator/precision/case-{ident}-{d}.json')['rates'])[:,0].min() for d in [70,90]]
        vals.append(read(root/f'oscillator/upstream/case-{ident}-analysis.npz')['rates'][:,0].min())
        out.append(f'| D{ident} | '+' | '.join(f'{v:.9f}' for v in vals)+' |')
    ident=335;r=long['rates'][list(long['ids']).index(ident)];mins=[r[:,1].min()]
    mins += [np.array(js(root/f'oscillator/precision/case-{ident}-{d}.json')['rates'])[:,1].min() for d in [70,90]]
    mins.append(read(root/f'oscillator/upstream/case-{ident}-analysis.npz')['rates'][:,1].min())
    out+=['\nFor D335, the minimum **second** rate is '+', '.join(f'{v:.9f}' for v in mins)+' in the five-million-step, 70-digit, 90-digit and actual-library checks respectively. All are positive. This is numerical support for hyperchaos, not a theorem about infinite-time Lyapunov exponents.\n',
    '![Long-run and independent checks](assets/rk45-oscillator-evidence.png)\n',
    '## Two rosettes, one orbit\n',
    'The original map has exact reflection symmetry `(q,p,h)->(-q,-p,h)`. In all 7,999,999 consecutive pairs of the dense D335 portrait, velocity changes sign; its smallest absolute value is about 0.17918. The two visible pieces are therefore exchanged by one sampled trajectory. They are not two coexisting attractors or a new object created by reflection. The finite check does not certify a global band decomposition.\n',
    'In actual D335 library steps, the median angular distance from a half-turn is 0.031778 radians. This matches the known interpretation discussed by [Humphries and Christodoulou, Section 5.2](https://math.mcgill.ca/humphries/research/papers/psthetarep.pdf), not a newly discovered mechanism.\n',
    'All portraits use eight million actual consecutive states after warm-up, with no nonlinear deformation or smoothing. The animation rotates a fixed consecutive million-point tail through 72 rigid orthographic cameras. Every frame retains every point and uses equal units; color represents visit density. The view rotates, rather than the dynamical set changing shape.\n',
    '## Initial states and nearby parameters\n',
    '| Case | Positive / 24 solution starts | Positive / 6 initial steps | Positive / 32 neighbors | Neighbors with two positive rates |\n| --- | ---: | ---: | ---: | ---: |']
    for ident in [335,448,2006,1890]:
        counts=[]
        for g in ['start','step','neighbor']:
            m=(cohort['ids']==ident)&(cohort['groups']==g);counts.append(int(cohort['flags'][m].sum()))
        m=(cohort['ids']==ident)&(cohort['groups']=='neighbor');two=np.sum(cohort['flags'][m]&(cohort['rates'][m,:,1].min(1)>0))
        out.append(f'| D{ident} | {counts[0]} / 24 | {counts[1]} / 6 | {counts[2]} / 32 | {two} / 32 |')
    out+=['\nThe 252 runs include four radii (0.01, 1, 100, 1e6) at six phases, six alternative first steps (0.001, 0.01, 1, 3, 10, 100), 32 scrambled Sobol parameter neighbors within relative log width ±0.02, and exact zero for each setting. They use 100,000 warm-up and 300,000 measured steps in ten windows. Positive means every leading window exceeds 0.001 and no tested numerical period is detected; the two-positive column additionally requires every second window to exceed zero. Full regular, weak and failed results remain in the archive.\n',
    'At exact zero, proposed steps grow tenfold until the declared state cutoff is reached. This expected control excludes a compact global attractor for the unrestricted full-state map.\n',
    '![Every initial point is retained in the approach view](assets/rk45-oscillator-approach.png)\n',
    'For each main portrait, 32,768 scrambled Sobol starts fill actual `q,p in [-4,4]`, `h in [0.25,8]`. All 65,536 complete 1,024 steps. Median distance to an independent 200,000-point reference falls from 2.65124 to 0.00242786 for D335, and from 2.42906 to 0.00131081 for D448. Reference half-cloud medians are 0.00325826 and 0.00188043. These are distances to finite samples, not certified distances to an invariant set or basin-volume estimates.\n',
    '## Changing the norm: a deliberately different algorithm\n',
    'A conjugate complex equation with a complex absolute-value error rule has phase symmetry. We compare its amplitude/step quotient, starting at amplitude `1/sqrt(1-zeta^2)`, with the same polynomial pair and controller factors. This replaces the error rule; it is not a change of coordinates applied to the complete original algorithm. Norm dependence is established prior work.\n',
    '| Case | Mean quotient spectrum | Tested quotient period (0 = none found) |\n| --- | --- | ---: |']
    for r in norm['runs']:out.append(f'| D{r["id"]} | '+', '.join(f'{v:.9f}' for v in r['means'])+f' | {r["quotient_period"]} |')
    out+=['\nThese fresh comparisons measure a million steps in twenty windows after 100,000 warm-up. D335 has negative finite quotient rates but no period within 512: it is not labeled a proven regular attractor. D448 retains a weak positive quotient estimate. D2006 also has a small positive mean, with a minimum leading window below the original 0.001 screening threshold. The quotient has an additional exact neutral phase direction when lifted to the full complex state. A quotient period is not necessarily a closed full orbit.\n',
    '## Complete search, exclusions and reproducibility\n']
    damping=js(root/'oscillator/damping/report.json')
    out += ['### Declared damping slices through D335\n', '| Slice | Settings | Leading criterion | Two positive | Periodic |\n| --- | ---: | ---: | ---: | ---: |']
    for r in damping['groups']:out.append(f'| {r["group"]} | {r["count"]} | {r["positive"]} | {r["two_positive"]} | {r["periodic"]} |')
    out += ['\nEach slice contains 257 fresh starts at the same relative tolerance, with damping widths ±0.018 and ±0.0005 around D335. The protocol uses 100,000 warm-up and 400,000 measured steps in twenty windows. Both grids were fixed before execution, and the exact base setting is repeated as a deterministic control. Two-positive means every second window exceeds zero. The retained points do not certify an open hyperchaotic interval or its boundaries.\n']
    good=follow['counters'][:,0]==0;periods={str(k):int(np.sum(follow['period']==k)) for k in np.unique(follow['period']) if k}
    out += [f'The scout contains 2,048 Sobol settings plus three declared controls, seed 2026110181, damping in [0.001,0.999] and relative tolerance from 1e-12 to 1e-1. It measures 50,000 steps after 20,000 warm-up. Of 2,051 settings, 2,050 complete and 1,061 pass the leading-window criterion. D802 leaves the declared solution magnitude cutoff 1e10 during warm-up; this is an observed cutoff exit, not proof of divergence.\n',
    f'Every one of the 1,061 flags receives a fresh million-step follow-up: {int(good.sum())} complete, {int(follow["flags"].sum())} retain the criterion, and {int(np.sum(follow["flags"]&(follow["rates"][:,:,1].min(1)>.001)))} retain two rates above 0.001 in every window. Detected period counts: `{json.dumps(periods,sort_keys=True)}`. No outcome is omitted because it is inconvenient.\n',
    'The production Jacobian matches 256 central-difference checks within 1.88e-9 and actual one-step SciPy evaluations within 9.57e-13. The norm-control Jacobian and polynomial stage formulas have their own 512 checks. Exact rational certificates give accepted-solution bounds for D335, D448 and D2006; the sufficient condition fails for D1890 and no bound is claimed for it. These do not establish a full-state trapping region.\n',
    'The [frozen packet](../results/rk45-controller/) contains all screen/follow outcomes, all 24 initial gallery clouds and the 13 complete two-positive follow-up portraits, six long runs, eight precision runs, four actual-library trajectories, all controls, volume histories, exact portrait histograms, animation tails, source snapshots and SHA-256 manifests. Raw eight-million-state portrait arrays remain local with hashes and regeneration information. Duplicate blocks were compared with their public aggregate before omission. A pre-execution reporting syntax error and the initial conservative symbolic draft are preserved.\n',
    '```sh\npython research/verify_rk45_controller.py --root results/rk45-controller --assets docs/assets\npython research/audit_rk45_structure.py --cases results/rk45-controller/source/structure/input-cases.json --out rk45-audit.json\n```\n',
    'Production runs use Python 3.12, NumPy, SciPy 1.18.1 and Numba. The independent replay uses mpmath; the exact audit uses SymPy 1.14. Installed source bytes are archived because an online manual may describe another release.\n',
    'The [priority review](rk45-priority.md) states what was read and what remains inaccessible. The strongest defensible contribution is the documented geometry and validation of these specific numerical regimes. A claim that nobody has seen an equivalent attractor would exceed the evidence.\n']
    args.out.write_text('\n'.join(out),encoding='utf-8')
    print('Wrote',args.out,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
