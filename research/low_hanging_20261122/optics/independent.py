"""Cross-check the optical witnesses in complex arithmetic and optionally plot their detuning."""
import argparse
import json
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/low_hanging_20261122/optics"


def iterate(state, params):
    E, T = state
    P, a, d, k, g, h = params
    f = 1 - 1/(1+abs(E)**2)
    return P + a*E*mp.exp(mp.j*(d+k*f+T)), T+h*(g*f-T)


def main(figure):
    mp.mp.dps = 170
    rec = json.loads((ROOT / "results/wow_search_20261122/orbit/case-56729.json").read_text())
    names = ("pump", "field_retention", "base_detuning", "fast_shift", "memory_gain", "relaxation")
    p = [mp.mpf(rec["parameters"][name]) for name in names]
    proof = json.loads((OUT / "certificate.json").read_text())
    paths = {}
    for label, item in proof["witnesses"].items():
        x, y, T = [mp.mpf(int(v))*mp.power(2,k) for v,k in item["exact_start_as_mantissa_and_binary_exponent"]]
        E = mp.mpc(x, y)
        trace = [T]
        for n in range(1,4248):
            E, T = iterate((E,T),p)
            trace.append(T)
            if label == "minus" and n == item["sink_ellipsoid_entry_step"]:
                center = proof["cycles"]["sink"]["center_approx"]
                assert max(abs(E.real-center[0]), abs(E.imag-center[1]), abs(T-center[2])) < mp.mpf("1e-8")
            if label == "plus" and n >= 2200:
                assert mp.mpf("-8.71") < T < mp.mpf("-8.45")
                assert ((mp.mpf("-.46") < E.real < mp.mpf(".33")) if n%2 == 0
                        else (mp.mpf("1.88") < E.real < mp.mpf("2.82")))
        paths[label] = trace
    print("Complex replay passed for both trajectories through step 4247.")
    if not figure:
        return
    import matplotlib.pyplot as plt
    import numpy as np

    fig, (overview, detail) = plt.subplots(2,1,figsize=(9.4,5.4),height_ratios=[1.65,1],
                                         constrained_layout=True)
    for ax in (overview,detail):
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y',alpha=.18)
        ax.set_xlabel("Round trips")
    colors = {"minus":"#16697a", "plus":"#b1462b"}
    for label, desc in [("minus","launch reaching the 10-cycle"),("plus","launch reaching the alternating lobes")]:
        overview.plot(np.arange(4248),np.array(paths[label],dtype=float),lw=1.15,color=colors[label],label=desc)
    overview.set(xlim=(0,4247),ylabel="Detuning T",title="Two launches near the same period-two saddle")
    overview.legend(loc="center right",frameon=False,fontsize=9)
    end = 1150
    separation = [float(abs(a-b)) for a,b in zip(paths['plus'][:end],paths['minus'][:end])]
    detail.semilogy(np.arange(end),separation,color="#473a5f",lw=1.2)
    detail.set(xlim=(0,end-1),ylim=(1e-29,20),ylabel="Difference in T",yticks=[1e-28,1e-21,1e-14,1e-7,1])
    fig.savefig(OUT / "detuning-split.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure",action="store_true")
    main(parser.parse_args().figure)
