"""Choose diverse finite portraits against archived reference shape descriptors."""
import argparse
import csv
import json
import os
import sys
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from PIL import Image, ImageDraw
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from compare_attractor_shapes import fingerprint, distance


def draw(paths, ids, output, columns=6):
    tile = 210
    canvas = Image.new('RGB', (columns*tile, (len(ids)+columns-1)//columns*tile), '#f6f5f0')
    pen = ImageDraw.Draw(canvas)
    for i, case in enumerate(ids):
        x = i % columns * tile
        y = i // columns * tile
        points = np.load(paths[case])['points']
        points = points[::max(1,len(points)//40000)]
        center = points.mean(0)
        value, axes = np.linalg.eigh(np.cov(points.T))
        screen = (points-center) @ axes[:, ::-1][:, :2]
        scale = np.max(np.abs(screen))*1.02 + 1e-10
        pixels = np.clip(((screen / scale + 1) * (tile-2)/2).astype(int), 0, tile-2)
        density = np.zeros((tile, tile))
        np.add.at(density,(pixels[:,1], pixels[:,0]),1)
        shade = np.log1p(density) / np.log1p(max(2,np.percentile(density[density>0],99.7)))
        shade = np.minimum(1,shade)**.55
        rgb = np.asarray([246,245,240])- shade[:,:,None]*np.asarray([222,177,159])
        canvas.paste(Image.fromarray(rgb.astype('uint8')), (x,y))
        pen.text((x+3,y+2), str(case), fill='black')
    canvas.save(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--in', dest='sources', type=Path, nargs='+', required=True)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--retain', type=int, default=54)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    with np.load(args.references) as archive:
        refs = {k:archive[k] for k in archive.files}
    paths = {}
    for source in args.sources:
        for f in source.glob('case-*.npz'):
            paths[int(f.stem.split('-')[-1])] = f
    kept = {}
    fp = {}
    for case,path in sorted(paths.items()):
        with np.load(path) as f:
            z, rates = f['points'],f['rates']
        if min(rates) <= .008:
            continue
        fp[case] = fingerprint(z)
        nearest = min((float(distance(fp[case], r)), name) for name,r in refs.items())
        within = float(distance(fingerprint(z[:len(z)//2]), fingerprint(z[len(z)//2:],1)))
        kept[case] = dict(nearest=nearest, within=within, min_rate=float(min(rates)))
    candidates = sorted(kept, key=lambda i: kept[i]['nearest'][0] - .45*kept[i]['within'], reverse=True)
    ids = []
    n = min(len(candidates),args.retain)
    diversity = {i: .7 for i in candidates}
    while len(ids) < n:
        def score(i):
            return min(diversity[i], .65*kept[i]['nearest'][0]) - .22*kept[i]['within']
        winner = max(candidates, key=score)
        candidates.remove(winner)
        ids.append(winner)
        for i in candidates:
            diversity[i] = min(diversity[i],distance(fp[winner],fp[i]))
    (args.out/'selection.json').write_text(json.dumps(dict(ids=ids, reference_count=len(refs), rows=kept),indent=2)+'\n')
    draw(paths, ids, args.out/'contact.png')
    print(ids)


if __name__ == '__main__':
    main()
