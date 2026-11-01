"""Audit the fold certificate's strict interval inequalities and actual-state figure."""
import json
from decimal import Decimal
from pathlib import Path
import numpy as np
from verify_raar_folding import sha,interval,sign
def main():
    repo=Path(__file__).resolve().parents[1];root=repo/'results/raar-fold-certificate';load=lambda p:json.loads((root/p).read_text());files=load('manifest.json')['files'];assert {r['path'] for r in files}=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
    for r in files:p=root/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==r['bytes']
    p=load('proposal/input.json');r=load('proof/report.json');assert r['status']=='passed' and r['input_sha256']==sha(root/'proof/input.json') and sha(root/'proof/input.json')==sha(root/'proposal/input.json')
    a,b=interval(r['t_box']);c,d=interval(r['interval_newton']);assert 0<a<c<=d<b<1;sign(r['line_derivative']);sign(r['nonzero_minor']);sign(r['kernel_transversality']);assert all(sign(v)>0 for v in r['axis_distances']);assert [sign(v) for v in r['endpoint_determinants']]==[-1,1];assert interval(r['distance_to_reference_index_993053'])[1]<Decimal('.000042');assert Decimal(r['independent_hessian_max_residual'])<Decimal('1e-150')
    with np.load(repo/'results/raar-circles/render/case-1051-animation.npz') as f:z=f['states'];assert np.array_equal(f['params'],p['params'])
    assert np.array_equal(z[p['indices']],p['segment_endpoints']);assert sha(repo/'results/raar-circles/render/case-1051-animation.npz')==p['source_sha256'];assert p['indices']==[857629,993053]
    with np.load(root/'proposal/determinants.npz') as f:
        determinants=f['determinants'];assert len(determinants)==1000000 and not f['flags'].any();assert (determinants<0).sum()==p['negative_count']==748512;assert (determinants>0).sum()==p['positive_count']==251488;k=f['nearest_opposite_distances'].argmin();assert [int(f['negative_indices'][k]),int(f['nearest_positive_indices'][k])]==p['indices'];assert f['nearest_opposite_distances'][k]==p['pair_distance']
    with np.load(root/'render/sign-histograms.npz') as f:
        xy=(z-f['center'])@f['projection']
        for panel,points in enumerate([xy[:-1],xy[1:]]):
            for j,s in enumerate([-1,1]):
                q=points[determinants[:-1]*s>0];h,_,_=np.histogram2d(q[:,0],q[:,1],bins=(1000,700),range=[f['x_bounds'],f['y_bounds']]);assert np.array_equal(h.T[::-1],f['histograms'][panel,j])
        assert np.all(f['histograms'].sum((1,2,3))==999999)
    assert sha(root/'render/raar-certified-fold.png')==sha(repo/'docs/assets/raar-certified-fold.png');print('PASS: strict fold/interval-Newton inequalities, exact source endpoints, sample proximity and all 999999 plotted input/output pairs')
if __name__=='__main__':main()
