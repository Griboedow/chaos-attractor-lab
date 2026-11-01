"""Check finite geometric rankings, controls, source bytes and every plotted state."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from compare_raar_catalog import sha
from compare_attractor_shapes import fingerprint,distance,normalize
from compare_affine_shapes import whiten

def main():
    repo=Path(__file__).resolve().parents[1];root=repo/'results/raar-catalog-comparison';load=lambda p:json.loads(p.read_text())
    manifest=load(root/'manifest.json')['files'];assert {r['path'] for r in manifest}=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
    for row in manifest:p=root/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes']
    report=load(root/'comparison/report.json');source=repo/'results/raar-circles/render/case-1051-animation.npz';assert sha(source)==report['source_sha256']
    with np.load(source) as f:z=f['states']
    refs={}
    for name,digest in report['reference_files'].items():
        path=repo/'results/attractor-novelty/reference/trajectories'/name;assert sha(path)==digest;refs[path.stem]=np.load(path)['points']
    assert len(refs)==108
    for label,transform in [('similarity',normalize),('affine',whiten)]:
        fp=fingerprint(transform(z));ranks=sorted((float(distance(fp,fingerprint(transform(v)))),n) for n,v in refs.items());expected=report['representations'][label]['all_rankings'];assert [n for d,n in ranks]==[r['reference'] for r in expected];np.testing.assert_allclose([d for d,n in ranks],[r['distance'] for r in expected],rtol=0,atol=1e-13)
        within=distance(fingerprint(transform(z[:500000])),fingerprint(transform(z[500000:]),1));assert abs(within-report['representations'][label]['within_time_halves_fingerprint'])<1e-13
        for expected in report['representations'][label]['subsample_rankings']:
            seed=expected['seed'];rng=np.random.default_rng(seed);fp=fingerprint(transform(z[rng.choice(len(z),6000,False)]),seed);ranks=sorted((float(distance(fp,fingerprint(transform(v[rng.choice(len(v),6000,False)]),seed))),n) for n,v in refs.items());assert ranks[0][1]==expected['nearest'] and abs(ranks[0][0]-expected['distance'])<1e-13
        with np.load(root/'comparison/aligned-samples.npz') as f:
            for item in report['representations'][label]['nearest']:
                key=label+'--'+item['reference'];a=f[key+'--candidate'];b=f[key+'--reference'];rms=np.sqrt((np.mean(cKDTree(b).query(a)[0]**2)+np.mean(cKDTree(a).query(b)[0]**2))/2);assert abs(rms-item['best_local_icp_rms'])<1e-13
    linear=np.array([[2.,.7,-.2],[0.,.4,.1],[.3,0.,1.1]]);control=distance(fingerprint(whiten(z)),fingerprint(whiten(z@linear+[3,-7,11])));assert control<1e-12 and abs(control-report['exact_affine_control'])<1e-13
    render=load(root/'render/render.json')
    with np.load(root/'render/plotted-clouds.npz') as f:
        for item in render['columns']:
            i=item['column'];transform=normalize if item['representation']=='similarity' else whiten;a=transform(z);b=transform(refs[item['reference']]);aa=(a[f[f'candidate_indices_{i}']]-item['candidate_center'])/item['candidate_scale']@np.array(item['rotation'])+item['translation'];bb=(b[f[f'reference_indices_{i}']]-item['reference_center'])/item['reference_scale'];assert np.array_equal(aa,f[f'candidate_{i}']) and np.array_equal(bb,f[f'reference_{i}']);assert len(aa)==len(bb)==8000
            for q in [aa,bb]:xy=q@f['projection'];assert np.all(xy>=f['low']) and np.all(xy<=f['high'])
    assert sha(root/'render/raar-catalog-comparison.png')==sha(repo/'docs/assets/raar-catalog-comparison.png')
    print('PASS: 108 reference hashes, complete and resampled rankings, temporal/affine controls, six alignment residuals and all 48000 plotted states')

if __name__=='__main__':main()
