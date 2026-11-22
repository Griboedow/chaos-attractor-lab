"""Index the finite search, five raw numerical arrays, published image and independent controls."""
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image

from search import RAW, OUTPUT
from selected import CASE


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def summarize_array(name, key, recipe):
    path = RAW / name
    with np.load(path) as src:
        description = {k: dict(shape=list(src[k].shape), dtype=str(src[k].dtype)) for k in src.files}
    return dict(cloud_relative='boundary/' + name, key=key, shape=description[key]['shape'],
                arrays=description, bytes=path.stat().st_size, sha256=digest(path), recipe=recipe)


def main():
    search_path = RAW / 'candidate-search.json'
    screen = load_json(search_path)
    fields = ('seeds', 'converged', 'parameter_clusters_including_lower_period', 'checked', 'passed')
    counts = {field: sum(m[field] for m in screen['counts']) for field in fields}
    assert len(screen['rows']) == counts['passed']
    primary = summarize_array('case-6902-stratified-inverse-sphere.npz', 'points', 'inverse.py then stratified.py')
    sphere = summarize_array('case-6902-inverse-sphere.npz', 'points', 'inverse.py --depth 8')
    plane = summarize_array('case-6902-inverse-plane.npz', 'points', 'inverse.py --depth 8')
    forward = summarize_array('case-6902-forward.npz', 'labels', 'inverse.py --depth 8')
    microscopic = summarize_array('case-6902-microscopic.npz', 'labels', 'microscope.py --width 1000 --height 800')
    inverse_check = load_json(OUTPUT / 'inverse-check.json')
    stratified = load_json(OUTPUT / 'stratified-inverse.json')
    for data, record in [(primary, stratified), (sphere, inverse_check['sphere']),
                         (plane, inverse_check['plane']), (forward, inverse_check['forward'])]:
        assert (data['sha256'], data['shape']) == (record['sha256'], record['shape'])
    png = OUTPUT / 'fifteen-newton-basins.png'
    with Image.open(png) as im:
        dimensions = list(im.size)
        im.verify()
    records = ['certificate.json', 'independent-decimal.json', 'inverse-check.json',
               'stratified-inverse.json', 'microscopic-check.json', 'scale-control.json']
    result = dict(case=CASE, system='Newton map of p_a(z) = z^6 - az + a',
                  parameter='Unique zero of N_a^9(0) in the parameter disk in certificate.json',
                  search=dict(**counts, degree=[4,5,6], period=list(range(3,10)), settings_per_cell=8192,
                              cluster_accounting='Totals per degree/period cell. Counts include lower-period returns and repeat appearances across cells.',
                              search_file='boundary/' + search_path.name, sha256=digest(search_path)),
                  primary_inverse_cloud=primary, full_depth_eight_clouds=[sphere, plane],
                  forward_basin_grids=[forward, microscopic],
                  image=dict(path=png.name, shape_width_height=dimensions, sha256=digest(png),
                             drawn='2026-11-22'),
                  controls={name:digest(OUTPUT / name) for name in records})
    (OUTPUT / 'manifest.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
