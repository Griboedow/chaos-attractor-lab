"""Locate covering-window centers and show the two certified return routes."""
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
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.patches import Circle, FancyArrowPatch


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    raw = json.loads(args.input.read_text(encoding='utf-8'))
    with np.load(args.states) as data:
        states = data['states']
        assert data['params'].tolist() == raw['params'] and data['points'].tolist() == raw['points']
    meta = json.loads(args.metadata.read_text(encoding='utf-8'))
    view = np.array(meta['basis'])@np.array(meta['projection']).T
    center = np.array(meta['center'])
    projected = (states-center)@view
    windows = (np.array(raw['centers'])-center)@view
    half = meta['halfwidth']
    hist, xe, ye = np.histogram2d(projected[:, 0], projected[:, 1], bins=1100, range=[[-half, half]]*2)
    ink, muted, orange, blue = '#233f4b', '#71848b', '#c66c31', '#2b7c9a'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'text.color': ink, 'axes.labelcolor': ink, 'axes.edgecolor': '#cbd4d7', 'xtick.color': muted, 'ytick.color': muted})
    fig = plt.figure(figsize=(14, 6.4), facecolor='white')
    grid = fig.add_gridspec(1, 2, width_ratios=[1.04, 1], left=.055, right=.97, top=.81, bottom=.21, wspace=.17)
    ax = fig.add_subplot(grid[0])
    bx = fig.add_subplot(grid[1])
    cmap = LinearSegmentedColormap.from_list('cloud', ['#ffffff', '#bfd0d3', '#678b95', '#355965'])
    cap = np.quantile(hist[hist > 0], .995)
    ax.imshow(hist.T, origin='lower', extent=[-half, half, -half, half], cmap=cmap, norm=PowerNorm(.30, vmin=0, vmax=cap), interpolation='nearest')
    ax.scatter(windows[1:, 0], windows[1:, 1], color=orange, s=35, edgecolors='white', linewidths=.7, zorder=4)
    ax.scatter(*windows[0], color=blue, s=100, marker='D', edgecolors='white', linewidths=1, zorder=5)
    ax.annotate(r'$N_0$', windows[0], xytext=(12, 8), textcoords='offset points', color=blue, fontsize=12, fontweight='bold')
    ax.set(xlim=(-half, half), ylim=(-half, half), xlabel='View coordinate 1', ylabel='View coordinate 2')
    ax.set_aspect('equal')
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_title('A   Location in the numerical cloud', loc='left', fontsize=13, pad=15)
    bx.set(xlim=(0, 1), ylim=(0, 1))
    bx.axis('off')
    bx.set_title('B   Two repeatable routes', loc='left', fontsize=13, pad=15)
    positions = {'0': (.20, .64), '1': (.66, .64), '22': (.66, .26)}
    for name, pos in positions.items():
        color = blue if name == '0' else orange
        bx.add_patch(Circle(pos, .065, facecolor='white', edgecolor=color, linewidth=2))
        bx.text(*pos, r'$N_{'+name+'}$', ha='center', va='center', color=color, fontsize=13)
    arrow = lambda a, b, color, rad=0: bx.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=14, color=color, linewidth=1.8, connectionstyle=f'arc3,rad={rad}'))
    arrow((.25, .69), (.61, .69), orange, -.18)
    arrow((.66, .565), (.66, .335), orange)
    arrow((.60, .26), (.17, .57), orange, -.25)
    arrow((.15, .69), (.25, .69), blue, -2.4)
    bx.text(.20, .96, 'Stay at the same window', ha='center', va='center', color=blue, fontsize=11)
    bx.text(.20, .88, r'23 applications of $G$', ha='center', color=muted, fontsize=10)
    bx.text(.78, .46, '20 intermediate\nwindows', ha='center', va='center', color=muted, fontsize=10)
    bx.text(.36, .23, 'Return excursion', color=orange, ha='center', fontsize=11)
    bx.text(.44, .09, r'$G=F^4$   ·   Each route takes 92 original steps', ha='center', color=ink, fontsize=11)
    fig.text(.055, .945, 'Certified branching in the four-wing candidate', fontsize=22)
    fig.text(.055, .893, 'Product alphabet · Case 988 · a compact invariant subsystem with positive entropy', fontsize=12, color=muted)
    fig.text(.055, .105, 'Window centers are enlarged; most overlap near N0. Background: two million computed states.', fontsize=10, color=muted)
    fig.text(.055, .063, 'The right panel is a schematic of covering relations. It does not draw continuous trajectories or window dimensions.', fontsize=10, color=muted)
    fig.text(.055, .021, 'The certificate proves branching recurrence. Attraction and equality with the background cloud remain unproved.', fontsize=10, color=muted)
    args.image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.image, dpi=180, facecolor='white')
    plt.close(fig)
    subprocess.run([sys.executable, str(root/'tools/render_dates.py'), 'record', args.image.resolve().relative_to(root).as_posix(), '--role', 'diagnostic'], check=True)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    (args.out/'figure.png').write_bytes(args.image.read_bytes())
    report = dict(status='rendered', inputs={str(p): sha(p) for p in [args.input, args.states, args.metadata]}, image_sha256=sha(args.image), states=len(states), windows=len(windows), camera=view.tolist(), center=center.tolist(), histogram_cap=float(cap), scope='Diagnostic location of exact window centers over a numerical cloud, with an orthogonal equal-unit projection; symbols enlarged. Schematic return graph, not a continuous trajectory or attraction proof.')
    (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['input', 'states', 'metadata', 'image', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    main(p.parse_args())
