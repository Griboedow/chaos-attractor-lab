"""Locate all ten attracting phases next to the same optical saddle."""
import importlib.util
from pathlib import Path

import numpy as np
from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parents[3]
OLD = ROOT / "research/low_hanging_20261122/optics/verify.py"
spec = importlib.util.spec_from_file_location("cycle_certificate", OLD)
cycle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cycle)

# Columns approximate the expanding eigenvector and the real and imaginary
# parts of a stable eigenvector. Interval differentiation validates the frame.
EIGENFRAME = [
    ["-0.6342749658241229950287403482176795619600125169032560481318752017", "1", "0"],
    ["-0.5992321638468208834045237392137839474855161025076918549849558808",
     "0.199779632457877883809007953627244532152486389150760300054057484",
     "0.4774318212534440822348117450763334888402136839446826864431533788"],
    ["0.4884834506308933468252232703471826766189768524250065340515455",
     "0.03076749766783404148215354191526300524586216721512185214800194738",
     "0.05089509771793398501631454859689360662794138988158674198370264081"],
]


def column(items):
    return arb_mat([[v] for v in items])


def centered(radius):
    return arb(0, radius.upper())


def axes():
    return arb_mat([[arb(t) for t in row] for row in EIGENFRAME])


def complex_step(state, p, tangent=False):
    """Independent complex-field formula and directional derivatives."""
    x, y, detuning = state
    pump, a, phase, fast, gain, relax = p
    E = acb(x, y)
    intensity = x**2+y**2
    response = intensity/(1+intensity)
    rotation = a*acb(0, phase+fast*response+detuning).exp()
    propagated = rotation*E
    shifted = pump+propagated
    result = [shifted.real, shifted.imag, detuning+relax*(gain*response-detuning)]
    if not tangent:
        return result
    columns = []
    for e, t in [(acb(1),0), (acb(0,1),0), (acb(0),1)]:
        dr = 2*(E.conjugate()*e).real/(1+intensity)**2
        dE = rotation*e+acb(0,1)*propagated*(fast*dr+t)
        columns.append([dE.real,dE.imag,t+relax*(gain*dr-t)])
    return result, arb_mat(list(zip(*columns)))


def fixed_points(p):
    points = {}
    for key in ("gap", "sink"):
        period, seed, *_ = cycle.SEEDS[key]
        center = cycle.center(seed, period, p)
        box, _ = cycle.certify(key, center, p)
        points[key] = (center, column(box))
    return points


def local_sections(p, saddle, model):
    # Coordinates are scaled so the local block is exactly [-1,1]^3.
    r, s = arb("1e-25"), arb("1e-44")
    basis = axes() * arb_mat([[r,0,0], [0,s,0], [0,0,s]])
    box = saddle + basis * column([centered(arb(1))]*3)
    image, first = model(box.entries(), p, tangent=True)
    _, second = model(image, p, tangent=True)
    derivative = second*first
    jac = basis.inv() * derivative * basis
    unstable = jac[0,0] - abs(jac[0,1]) - abs(jac[0,2])
    stable = max(sum(abs(jac[i,j]) for j in range(3)).upper() for i in (1,2))
    assert unstable > arb("1.13") and stable < arb("0.83"), (unstable,stable)
    # One backward-bounded orbit crosses the entire unit block and intersects
    # this transverse slice. The section and subsequent enclosures contain it.
    return [saddle + basis*column([arb(u), centered(arb(1)), centered(arb(1))])
            for u in ("-0.50", "-0.51")]


def transport(p, state, sink_box, metric, radius, model):
    """Move whole enclosures in QR frames chosen numerically and verified by Arb."""
    ref = column([a.mid() for a in state.entries()])
    basis = arb_mat([[1,0,0],[0,1,0],[0,0,1]])
    offsets = state-ref
    peak = arb(0)
    for step in range(1, 5001):
        state = ref + basis*offsets
        _, derivative = model(state.entries(),p,tangent=True)
        mid_image = column(model(ref.entries(),p))
        new_ref = column([a.mid() for a in mid_image.entries()])
        linear = derivative*basis
        frame, _ = np.linalg.qr(np.array(linear.tolist(),dtype=float))
        new_basis = arb_mat([[arb(float(x)) for x in row] for row in frame])
        inv = new_basis.inv()
        delta = inv*(mid_image-new_ref)+(inv*linear)*offsets
        offsets = column([centered(abs(t)) for t in delta.entries()])
        ref, basis = new_ref, new_basis
        true_box = ref+basis*offsets
        peak = max(peak, max(a.rad().upper() for a in true_box.entries()))
        if peak > arb("0.1"):
            raise ValueError(f"Enclosure too wide after {step} steps: {peak}")
        difference = true_box - sink_box
        if cycle.norm(difference)<radius and (difference.transpose()*metric*difference)[0,0]<radius**2:
            assert peak < arb("1e-12")
            return step, float(peak)
    raise ValueError("No entry")


def verify(bits, model):
    ctx.prec=bits
    p=cycle.parameters()
    points=fixed_points(p)
    slices=local_sections(p, points["gap"][1], model)
    H, radius=cycle.capture_metric(points["sink"][0],p)
    returns = [transport(p, section, points["sink"][1], H, radius, model) for section in slices]
    assert [n for n,_ in returns] == [2274,2345]
    phases = {(-n+2*j)%10 for n,_ in returns for j in range(5)}
    assert phases == set(range(10))
    print(bits, model.__name__, returns)


if __name__ == "__main__":
    verify(384, cycle.step)
    verify(512, complex_step)
