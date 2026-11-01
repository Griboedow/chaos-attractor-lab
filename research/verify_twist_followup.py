"""Check follow-up bytes and reconstruct finite numerical summaries from data."""
import hashlib,json
from pathlib import Path
import numpy as np
from twist_map import step,step_tangent


def main():
    root=Path(__file__).resolve().parents[1];packet=root/'results/twist-followup'
    manifest=json.loads((packet/'manifest.json').read_text());assets=json.loads((packet/'asset-manifest.json').read_text())
    for base,entries in [(packet,manifest['files']),(root,assets['files'])]:
        for entry in entries:
            path=base/entry['path'];assert path.stat().st_size==entry['bytes']
            assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256'],str(path)
    for ident,expected in [(1709,12314),(3254,0)]:
        d=np.load(packet/f'approach/case-{ident}.npz');assert len(d['initial'])==32768
        returns=np.array([np.max(np.linalg.norm(d['final24'][:,lag:]-d['final24'][:,:-lag],axis=2),axis=1) for lag in [1,2,3]]).T
        assert np.max(abs(returns-d['return_residual']))<1e-14
        assert np.count_nonzero(returns[:,2]<1e-8)==expected and not (returns[:,:2]<1e-8).any()
    d=np.load(packet/'approach/case-3254.npz');outliers=np.flatnonzero(d['nearest_reference_distance'][-1]>.02)
    follow=np.load(packet/'approach-outliers/histories.npz');assert np.array_equal(follow['ids'],outliers) and len(outliers)==52
    assert np.min(follow['windows'][:,:,0])>.2526
    for name in ['approach-outliers/histories.npz','unstable-seeds/case-1709.npz','unstable-seeds/case-3254.npz']:
        data=np.load(packet/name);p=data['parameters'];r=p[0] if p.ndim==1 else p[:,0,None]
        assert np.max(abs(data['windows'].sum(2)-3*np.log(r)))<1e-10
        assert (data['windows'][:,:,0]>0).all()
    roots=json.loads((packet/'periodic-skeleton/analysis.json').read_text())['runs']
    for row in roots:
        p=np.array(row['parameters'])
        for orbit in row['periodic_orbits']:
            z=np.array(orbit['states'][0]);q=z.copy()
            for _ in range(orbit['least_period']):q=step(q,p)
            assert np.linalg.norm(q-z)<2e-10
    splitting=json.loads((packet/'splitting/analysis.json').read_text())['runs']
    for row in splitting:
        data=np.load(packet/'splitting'/f"case-{row['id']}.compact.npz");idx=data['original_sample_indices'];angle=np.arcsin(np.minimum(1.,abs(np.sum(data['normal']*data['stable'],axis=1))))
        assert np.max(abs(angle-data['angle_radians'][idx]))<1e-12
        for length in [100,1000]:
            cumulative=np.r_[0.,np.cumsum(data['log_plane_area'])];minimum=np.min(cumulative[length:]-cumulative[:-length])/length
            expected=next(b['minimum_mean_log_area'] for b in row['blocks'] if b['length']==length)
            assert abs(minimum-expected)<1e-12
        for k in np.linspace(0,len(idx)-1,128).astype(int):
            _,j=step_tangent(data['states'][k],np.eye(3),data['parameters']);normal=data['normal'][k]
            actual=np.log(abs(np.linalg.det(j))*np.linalg.norm(np.linalg.solve(j.T,normal)))
            assert abs(actual-data['log_plane_area'][idx[k]])<1e-10
    prior=np.load(packet/'prior-replay/histories.npz');assert len(prior['starts'])==33
    assert (prior['return_residual'][:,0]<1e-8).all() and (prior['windows'][:,:,0]<0).all()
    for i in range(6):assert np.load(packet/f'prior-comparison/panel-{i}.npz')['histogram'].sum()==1000000
    print(json.dumps(dict(files=len(manifest['files']),assets=len(assets['files']),ensemble_starts=65536,outlier_continuations=52,unstable_seeds=512,prior_fixed_point_replays=33,checks='passed')))


if __name__=='__main__':main()
