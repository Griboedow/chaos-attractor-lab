"""Compare raw recurrent orbits in four fixed native and principal-axis views."""

import argparse
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
from morphology import tile
from model import trace

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("ids", type=int, nargs="+")
    a = ap.parse_args()
    q = np.load(a.cache / "screen.npz")
    params = q["params"].astype("f8")
    for j in range(0, len(a.ids), 6):
        cases = a.ids[j : j + 6]
        canvas = Image.new("RGB", (230 * 4, 167 * len(cases)), (16, 24, 34))
        d = ImageDraw.Draw(canvas)
        for row, k in enumerate(cases):
            out = a.cache / f"candidate-{k}.npz"
            p = params[k]
            init = np.array([0.173, -0.341, 0.271])
            if not out.exists():
                pts, rates = trace(p, init, 100000, 65536, 4, 131072)
                np.savez_compressed(
                    out, points=pts.astype("f4"), params=p, initial=init, rates=rates
                )
            z = np.load(out)
            pts = z["points"]
            rr = z["rates"].mean(0)
            print(k, np.round(rr, 4), flush=True)
            for m in range(4):
                canvas.paste(tile(pts, (226, 143), m), (230 * m + 2, 167 * row + 1))
            d.text(
                (4, row * 167 + 145),
                f"Case {k}    full spectrum {rr[0]:.4f}, {rr[1]:.4f}, {rr[2]:.4f}",
                fill="white",
            )
        canvas.save(a.cache / f"long-visual-{j // 6}.jpg", quality=95)
