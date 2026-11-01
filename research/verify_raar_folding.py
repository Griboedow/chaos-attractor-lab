"""Audit frozen bytes, complete controls, image counts and interval-certificate inequalities."""
import hashlib,json
from decimal import Decimal
from pathlib import Path
import numpy as np
from PIL import Image

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def interval(s):return tuple(Decimal(v.strip()) for v in s.strip('[]').split(','))
def sign(s):
    a,b=interval(s);assert a<=b
    if a>0:return 1
    if b<0:return -1
    raise AssertionError('Undetermined sign')
def variations(s):assert 0 not in s;return sum(a!=b for a,b in zip(s,s[1:]))
def main():
    repo=Path(__file__).resolve().parents[1];root=repo/'results/raar-folding';load=lambda p:json.loads((root/p).read_text());manifest=load('manifest.json')['files']
    assert {r['path'] for r in manifest}=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
    for r in manifest:p=root/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256']
    r=load('seeds/report.json');assert r['radii']==[1e-5,1e-7,1e-9] and r['seed_counts']==[4096,16384] and len(r['runs'])==6
    assert sum(v['completed'] for v in r['runs'])==61440
    for row in r['runs']:
        with np.load(root/'seeds'/(row['key']+'.npz')) as d:
            assert not d['status'].any() and len(d['status'])==row['completed']==row['seed_count'];assert np.array_equal(d['params'],r['params']);assert d['snapshots'].shape==(10,7,row['seed_count'],3);assert np.isfinite(d['snapshots']).all()
            assert np.array_equal(np.quantile(d['reference_distances'],[.5,.9,.99,1],axis=1).T,row['distance_quantiles']);assert np.array_equal(np.quantile(d['reference_split_distances'],[.5,.9,.99,1]),row['reference_split_quantiles'])
    with np.load(root/'seeds/evolution-histograms.npz') as d:old=d['histograms']
    with np.load(root/'hires/evolution-histograms.npz') as d:hi=d['histograms'];assert not d['clipped'].any()
    assert hi.shape==(65,1400,1800) and np.all(hi.sum((1,2))==114688);assert np.array_equal(hi.reshape(65,700,2,900,2).sum((2,4)),old)
    meta=r['histogram']
    with np.load(root/'seeds/radius-1e-07-n16384.npz') as d:
        for k,points in zip(d['returns'],d['snapshots']):
            xy=(points.reshape(-1,3)-meta['camera_center'])@np.array(meta['projection']);h,_,_=np.histogram2d(xy[:,0],xy[:,1],bins=(1800,1400),range=[meta['x_bounds'],meta['y_bounds']]);assert np.array_equal(h.T[::-1],hi[k])
    render=load('render/analysis.json');assert render['uniform_image_scale']==1 and render['density_area_factor']==4
    l,t,rr,b=render['fixed_pixel_crop'];assert np.array_equal(hi.sum((1,2)),hi[:,t:b,l:rr].sum((1,2)))
    with Image.open(root/'render/raar-cycle-unfolding.gif') as im:assert im.n_frames==65 and im.size==(1100,900)
    sturm=load('inverse-proof/report.json');ind=load('inverse-independent/report.json');assert sturm['status']==ind['status']=='passed';assert ind['total_preimages']==64
    with np.load(repo/'results/raar-circles/render/case-1051-animation.npz') as d:cloud=d['states'];params=d['params']
    assert [int(c['label'].split('-')[1]) for c in sturm['cases']]==np.linspace(10000,989999,25,dtype=int).tolist()
    counts=[]
    for case,other in zip(sturm['cases'],ind['cases']):
        assert case['label']==other['label'];index=int(case['label'].split('-')[1]);assert np.array_equal(case['target'],cloud[index]);assert np.array_equal(case['parameters'],params)
        total=variations(case['sturm_negative_infinity_signs'])-variations(case['sturm_positive_infinity_signs']);assert total==case['total_real_roots']==len(case['roots']);sign(case['projected_circle_target_exclusion']);[sign(v) for v in case['pi_chart_exclusions']]
        valid={}
        for rootcase in case['roots']:
            assert variations(rootcase['sturm_left_signs'])-variations(rootcase['sturm_right_signs'])==1
            good=sign(rootcase['radial_A'])>0 and sign(rootcase['radial_B'])>0;assert good==rootcase['valid_preimage']
            if good:valid[rootcase['root']]=rootcase
        counts.append(len(valid));assert len(valid)==case['number_of_preimages']==other['number_of_preimages'];assert sorted(valid)==sorted(v['sturm_root'] for v in other['preimages'])
        for x in other['preimages']:
            assert interval(x['contraction'])[1]<1;assert all(sign(v)>0 for v in x['axis_distances'])
            for outer,box,k in zip(valid[x['sturm_root']]['preimage_box'],x['box'],x['krawczyk']):
                a,b=interval(outer);c,d=interval(box);e,f=interval(k);assert a<c<e<=f<d<b
    assert counts.count(2)==18 and counts.count(4)==7
    for folder,names in [('render',['raar-cycle-unfolding.gif','raar-cycle-unfolding.png','raar-cycle-seed-controls.png']),('inverse-render',['raar-four-preimages.png'])]:
        for name in names:assert sha(root/folder/name)==sha(repo/'docs/assets'/name)
    print(json.dumps(dict(status='passed',files=len(manifest),seed_trajectories=61440,complete_targets=25,certified_preimages=64,all_render_histograms_checked=65)))
if __name__=='__main__':main()
