"""Diagnostic of certified orientation signs, without drawing invented orbits."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from audit_product_orientation import interval


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    reports = [json.loads(p.read_text(encoding='utf-8')) for p in [args.primary, args.independent]]
    assert all(r['status'] == 'passed' for r in reports)
    ink, muted, orange, blue = '#233f4b', '#71848b', '#c66c31', '#2b7c9a'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'text.color': ink,
                         'axes.labelcolor': ink, 'axes.edgecolor': '#cbd4d7', 'xtick.color': muted, 'ytick.color': muted})
    fig = plt.figure(figsize=(14, 6.9), facecolor='white')
    grid = fig.add_gridspec(1, 2, left=.045, right=.965, bottom=.24, top=.79, width_ratios=[.82, 1.18], wspace=.23)
    ax, bx = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    ax.set(xlim=(-1.38, 1.38), ylim=(-1.3, 1.3), aspect='equal')
    ax.axis('off')
    angles = np.pi/2-2*np.pi*np.arange(23)/23
    pos = np.column_stack([np.cos(angles), np.sin(angles)])
    for i in range(23):
        color = orange if i in [4, 5, 6] else blue
        a, b = pos[i], pos[(i+1) % 23]
        direction = (b-a)/np.linalg.norm(b-a)
        ax.add_patch(FancyArrowPatch(a+.065*direction, b-.065*direction, arrowstyle='-|>', mutation_scale=9, linewidth=1.4, color=color))
        ax.add_patch(Circle(a, .056, facecolor=color, edgecolor='white', linewidth=.7, zorder=4))
        ax.text(*(a*1.16), str(i), ha='center', va='center', color=color, fontsize=9)
    ax.text(0, .34, '23 windows', ha='center', fontsize=18)
    ax.text(0, .10, '92 original iterations', ha='center', fontsize=12, color=muted)
    ax.text(0, -.20, '3 orientation reversals', ha='center', fontsize=13, color=orange)
    ax.text(0, -.44, 'One reversed return', ha='center', fontsize=14, fontweight='bold', color=orange)
    ax.set_title('A   Certified return route', loc='left', fontsize=13, pad=18)
    for j, report in enumerate(reports):
        for row in report['windows']:
            lo, hi = map(float, interval(row['determinant']))
            color = orange if row['sign'] == -1 else blue
            x = row['window']+(-.11 if j == 0 else .11)
            bx.plot([x, x], [lo, hi], color=color, alpha=1 if j == 0 else .4, linewidth=2.7)
            bx.plot([x], [(lo+hi)/2], 'o' if j == 0 else 's', color=color, alpha=1 if j == 0 else .4, markersize=3)
    bx.axhline(0, color=muted, linestyle='--', linewidth=1)
    bx.set(xlim=(-.7, 22.7), ylim=(-3.25, 4.9), xlabel='Source window index', ylabel=r'Certified bounds for $\det D(F^4)$')
    bx.set_xticks(np.arange(0, 23, 2))
    bx.spines[['top', 'right']].set_visible(False)
    bx.grid(axis='y', color='#e8edef', linewidth=.7)
    bx.set_axisbelow(True)
    bx.set_title('B   Every interval stays away from zero', loc='left', fontsize=13, pad=18)
    bx.text(.02, .97, 'Solid: factorized posterior\nPale: full eight-site check', transform=bx.transAxes, va='top', color=muted, fontsize=10)
    fig.text(.045, .944, 'A reversed return inside the certified chaotic subsystem', fontsize=21)
    fig.text(.045, .886, 'Product alphabet · Case 988 · least period 92 · two interval formulas and an exact audit', fontsize=12, color=muted)
    fig.text(.045, .155, r'$\det D(F^{92})(p)<0$  excludes neighborhood conjugacy to orientation-constant local diffeomorphisms.', fontsize=12)
    fig.text(.045, .100, 'This includes standard three-dimensional Henon maps. Flow time maps and smooth Euclidean flow return maps are also excluded.', fontsize=10, color=muted)
    fig.text(.045, .048, 'The circle is a schematic itinerary, not a trajectory. Attraction, cloud identification and intrinsic equivalence remain open.', fontsize=10, color=muted)
    fig.savefig(args.image, dpi=180, facecolor='white')
    plt.close(fig)
    root = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, str(root/'tools/render_dates.py'), 'record', args.image.resolve().relative_to(root).as_posix(), '--role', 'diagnostic'], check=True)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    (args.out/'figure.png').write_bytes(args.image.read_bytes())
    result = dict(status='rendered', inputs={str(p): sha(p) for p in [args.primary, args.independent]}, image_sha256=sha(args.image),
                  scope='Schematic itinerary and interval determinant bounds. Neither a continuous orbit nor the dimensions of covering windows are depicted. Floating conversion is for plotting; exact endpoints are retained in the certificates.')
    (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['primary', 'independent', 'image', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    main(p.parse_args())
