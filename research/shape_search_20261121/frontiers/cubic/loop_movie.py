"""Animate the two-gain material loop and the three phase-matched attracting orbits."""

import argparse
import json
from pathlib import Path

import numpy as np
from cubic_material_lines import (
    CASE,
    COLORS,
    CUTS,
    GAINS,
    SOURCE,
    advance,
    ellipse,
    evolve_block,
    sample,
)
from matplotlib.font_manager import FontProperties, findfont
from model import forcing
from PIL import Image, ImageDraw, ImageFont

OUT = Path(__file__).resolve().parents[4] / "results/cubic-material-lines"
W, H, AA = 580, 798, 2
LEFT, RIGHT = 30, W - 30
SCALE = (RIGHT - LEFT) / 4.9
XMIN, YMAX = -2.45, 1.35
PLANE_HEIGHT = round(2 * YMAX * SCALE)
INK, PAPER, MUTED = "#1c303d", "#fbfaf5", "#4f626a"
ZONES = ["#e7edf4", "#eef1ee", "#f9eade"]
SINK_COLORS = ["#884fa7", "#107f6e", "#bd493e"]
FRAMES_PER_TURN = 128
STEPS, BLOCKS = 2048, 256


def prepare_fonts():
    normal = findfont(FontProperties(family="DejaVu Sans"))
    bold = findfont(FontProperties(family="DejaVu Sans", weight="bold"))
    return {
        name: ImageFont.truetype(font, size * AA)
        for name, font, size in [
            ("title", bold, 27),
            ("gain", bold, 21),
            ("metric", normal, 18),
            ("caption", normal, 16),
            ("counter", bold, 17),
        ]
    }


class Canvas:
    def __init__(self, fonts):
        self.fonts = fonts
        self.im = Image.new("RGB", (AA * W, AA * H), PAPER)
        self.draw = ImageDraw.Draw(self.im)

    def text(self, x, y, value, style, color=INK, anchor="la"):
        self.draw.text(
            (AA * x, AA * y), value, font=self.fonts[style], fill=color, anchor=anchor
        )

    def rect(self, x0, y0, x1, y1, fill):
        self.draw.rectangle(
            (round(x0 * AA), round(y0 * AA), round(x1 * AA), round(y1 * AA)), fill=fill
        )

    def dot(self, x, y, radius, fill, edge=PAPER):
        self.draw.ellipse(
            (
                (x - radius) * AA,
                (y - radius) * AA,
                (x + radius) * AA,
                (y + radius) * AA,
            ),
            fill=fill,
            outline=edge,
            width=round(1.4 * AA),
        )


def pixels(z, top):
    return LEFT + (z.real - XMIN) * SCALE, top + (YMAX - z.imag) * SCALE


def plane(canvas, top, polygon, sinks):
    xcuts, _ = pixels(np.array(CUTS), top)
    for x0, x1, color in zip([LEFT, *xcuts], [*xcuts, RIGHT], ZONES):
        canvas.rect(x0, top, x1, top + PLANE_HEIGHT, color)
    for x in xcuts:
        for y in range(top + 2, top + PLANE_HEIGHT - 4, 9):
            canvas.rect(x - 0.4, y, x + 0.4, y + 3, "#96a5a8")
    info, arc, closed = sample(polygon)
    x, y = pixels(closed, top)
    sx, sy = pixels(sinks, top)
    if (
        min(x.min(), sx.min()) < LEFT + 2
        or max(x.max(), sx.max()) > RIGHT - 2
        or min(y.min(), sy.min()) < top + 2
        or max(y.max(), sy.max()) > top + PLANE_HEIGHT - 2
    ):
        raise ValueError("A displayed path would leave the common physical scale")
    points = np.rint(AA * np.column_stack((x, y))).astype(int)
    # Dropping consecutive copies of a screen pixel preserves the rasterized path.
    points = points[np.r_[True, np.any(points[1:] != points[:-1], axis=1)]]
    canvas.draw.line(points.ravel().tolist(), fill=INK, width=4 * AA)
    for cx, cy, color in zip(sx, sy, SINK_COLORS):
        canvas.dot(cx, cy, 6.1, color)

    # Equal arc length per displayed pixel. Average eight subpixel zone samples so
    # thin late folds survive downsampling, even when several land on one pixel.
    width = (RIGHT - LEFT) * AA
    distances = (np.arange(width * 4) + 0.5) * arc[-1] / (width * 4)
    zones = np.digitize(np.interp(distances, arc, closed.real), CUTS)
    palette = np.array([tuple(bytes.fromhex(color[1:])) for color in COLORS])
    strip = np.rint(palette[zones].reshape(width, 4, 3).mean(axis=1)).astype(np.uint8)
    strip = Image.fromarray(strip[None, :, :]).resize(
        (width, 11 * AA), Image.Resampling.NEAREST
    )
    canvas.im.paste(strip, (LEFT * AA, (top + PLANE_HEIGHT + 6) * AA))
    return info


def render_frame(polygons, sinks, frame, fonts, reference):
    canvas = Canvas(fonts)
    canvas.text(LEFT, 13, "One loop, two gains", "title")
    t = frame / FRAMES_PER_TURN
    canvas.text(RIGHT, 22, f"{t:.2f} periods", "counter", anchor="ra")
    for i, color in enumerate(SINK_COLORS):
        canvas.dot(LEFT + 5 + 12 * i, 61, 5.1, color)
    canvas.text(LEFT + 46, 49, "Attracting orbits", "caption", MUTED)
    canvas.rect(281, 59, RIGHT, 61, "#d5dcda")
    canvas.rect(281, 59, 281 + (RIGHT - 281) * t / 2, 61, "#496775")
    for i in range(3):
        x = 281 + i * (RIGHT - 281) / 2
        canvas.rect(x - 0.65, 56, x + 0.65, 64, "#667b82")
    for j, (gain, top, polygon, marked) in enumerate(
        zip(GAINS, (111, 454), polygons, sinks)
    ):
        info = plane(canvas, top, polygon, marked)
        canvas.text(LEFT, top - 34, f"κ = {gain:.2f}", "gain")
        canvas.text(
            RIGHT,
            top - 31,
            f"Loop length  {info['length']:,.1f}",
            "metric",
            anchor="ra",
        )
        if frame % FRAMES_PER_TURN == 0:
            known = reference[j][frame // FRAMES_PER_TURN]
            if (
                abs(info["length"] - known["length"]) > 1e-6 * known["length"]
                or info["crossings"] != known["crossings"]
            ):
                raise ArithmeticError(
                    "Stroboscopic loop differs from the material-line packet"
                )
    canvas.text(
        LEFT,
        765,
        "Each strip follows the loop through the three shaded zones.",
        "caption",
        MUTED,
    )
    return canvas.im.resize((W, H), Image.Resampling.LANCZOS)


def main(out):
    row = next(
        r
        for r in json.loads((SOURCE / "certificates.json").read_text())["rows"]
        if r["case"] == CASE
    )
    survey = next(
        r
        for r in json.loads((SOURCE / "gradient.json").read_text())
        if r["case"] == CASE
    )
    reference = json.loads((OUT / "measurements.json").read_text())["resolutions"][-1][
        "runs"
    ]
    drive = forcing(np.linspace(0, 2 * np.pi, 2 * STEPS + 1), row)
    dt = 2 * np.pi / STEPS
    polygons = [ellipse(512) for _ in GAINS]
    sinks = []
    for gain in GAINS:
        run = next(r for r in survey["runs"] if r["gain"] == gain)
        points = np.array([complex(*z) for z in run["starts_last_period"]])
        sinks.append(points[np.argsort(points.real)])
    initial_sinks = [s.copy() for s in sinks]
    frames, fonts = [], prepare_fonts()
    for f in range(2 * FRAMES_PER_TURN + 1):
        if f:
            first = ((f - 1) % FRAMES_PER_TURN) * (STEPS // FRAMES_PER_TURN)
            stop = first + STEPS // FRAMES_PER_TURN
            for i, gain in enumerate(GAINS):
                for k in range(first, stop, STEPS // BLOCKS):
                    polygons[i] = evolve_block(
                        polygons[i], drive, gain, k, k + STEPS // BLOCKS, dt, 0.0175
                    )
                sinks[i] = advance(sinks[i], drive, gain, first, stop, dt)
        frames.append(render_frame(polygons, sinks, f, fonts, reference))
        if f % FRAMES_PER_TURN == 0:
            if (
                max(np.max(abs(s - start)) for s, start in zip(sinks, initial_sinks))
                > 1e-7
            ):
                raise ArithmeticError("A marked sink missed its one-period return")
            print(f"Verified turn {f // FRAMES_PER_TURN}", flush=True)
    durations = [50] * len(frames)
    durations[0], durations[-1] = 450, 1000
    frames[0].save(
        out,
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=0,
        format="WEBP",
        lossless=False,
        quality=84,
        method=4,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT / "loop.webp")
    main(parser.parse_args().out)
