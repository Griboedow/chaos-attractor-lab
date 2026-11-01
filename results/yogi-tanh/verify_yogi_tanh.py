"""Verify smooth-Yogi hashes, classifications, negative controls and render point accounting."""
import argparse,json
from pathlib import Path
import numpy as np
from PIL import Image
from verify_optimizer_followups import check_manifest,arrays,sha


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))


def check_spectra(root,folder,n):
    d=arrays(root/folder/'histories.npz');a=read(root/folder/'analysis.json')
    assert a['declared']==a['completed']==n and len(a['runs'])==n and np.all(d['status']==0)
    assert np.isfinite(d['windows']).all() and np.isfinite(d['tails']).all()
    for i,r in enumerate(a['runs']):
        w=d['windows'][i];z=d['tails'][i];tol=1e-9*(1+np.max(abs(z)))
        period=next((lag for lag in range(1,129) if np.max(abs(z[lag:]-z[:-lag]))<tol),None)
        lead=bool(w[:,0].min()>.005 and np.ptp(z[:,0])>1e-4 and period is None)
        assert r['period_up_to128']==period and r['leading']==lead and r['two']==(lead and w[:,1].min()>.001)
        assert np.allclose(w.mean(0),r['spectrum'],rtol=0,atol=1e-14) and np.allclose(w.min(0),r['window_minima'],rtol=0,atol=1e-14)
        assert np.max(abs(w.sum(1)-d['determinant_windows'][i]))<1e-7
    return d,a


def main(args):
    root=args.root;check_manifest(root);d=arrays(root/'scout/scout.npz');a=read(root/'scout/analysis.json')
    assert a['declared']==a['completed']==8256 and np.all(d['status']==0)
    assert np.isfinite(d['windows']).all() and np.isfinite(d['tails']).all()
    periods=[]
    for z in d['tails']:
        tol=1e-9*(1+np.max(abs(z)))
        periods.append(next((lag for lag in range(1,33) if np.max(abs(z[lag:]-z[:-lag]))<tol),0))
    assert np.array_equal(periods,d['periods'])
    good=(d['windows'].min(1)>.005)&(np.ptp(d['tails'][:,:,0],axis=1)>1e-4)&(d['periods']==0)
    assert np.array_equal(np.flatnonzero(good),a['primary_ids']) and good.sum()==488
    followed,fa=check_spectra(root,'follow',488);assert np.array_equal(followed['ids'],a['primary_ids']) and fa['leading']==471 and fa['two']==110
    long,la=check_spectra(root,'five-million',8);assert la['leading']==7 and la['two']==5
    co,ca=check_spectra(root,'cohorts',420)
    for r in ca['summary']:
        rows=[v for v in ca['runs'] if v['id']==r['id'] and v['kind']==r['kind']]
        assert len(rows)==r['total'] and sum(v['leading'] for v in rows)==r['leading'] and sum(v['two'] for v in rows)==r['two']
    precision=read(root/'precision/analysis.json');assert precision['declared']==precision['completed']==6
    for r in precision['runs']:
        w=arrays(root/'precision'/f"case-{r['id']}-{r['digits']}digits.npz")['windows']
        assert np.array_equal(w,np.array(r['windows'])) and np.allclose(w.min(0),r['minimum_windows'],rtol=0,atol=1e-14)
        assert np.all(w[:,0]>0)
        if r['id'] in [39,7239]:assert np.all(w[:,1]>0)
    for folder,analysis,n in [('upstream','upstream-analysis',6),('upstream-strong','upstream-strong-analysis',4)]:
        a=read(root/folder/'analysis.json');b=read(root/analysis/'analysis.json');w=arrays(root/analysis/'windows.npz')
        assert a['declared']==a['completed']==b['declared']==b['completed']==n
        for i,r in enumerate(a['runs']):
            ident=r['id'];f=root/folder/f'case-{ident}.npz';z=arrays(f)
            assert z['states'].shape==(300000,3) and np.isfinite(z['states']).all()
            assert np.max(abs(z['first20_physical_reference']-z['first20_physical_actual'])/(1+abs(z['first20_physical_reference'])))<1e-9
            row=next(t for t in b['runs'] if t['id']==ident);j=list(w['ids']).index(ident)
            assert row['source_sha256']==sha(f) and row['step_max_scaled_residual']<1e-10
            assert np.allclose(w['windows'][j].mean(0),row['spectrum'],rtol=0,atol=1e-14)
            assert np.max(abs(w['windows'][j].sum(1)-w['determinant_windows'][j]))<1e-7
            if ident in [39,7239]:assert np.all(w['windows'][j,:,:2]>0)
            if ident==2850:assert np.max(abs(z['states'][-1,:2]))<1e-100 and row['exact_fixed_line_linear_rates'] is not None
    provenance=read(root/'copy-provenance.json');split=provenance['lossless_gallery_split'];assert len(split)==36
    for r in split:
        f=root/r['destination'];g=arrays(f);assert sha(f)==r['sha256'] and g['states'].shape==(200000,3) and np.isfinite(g['states']).all() and int(g['ids'])==r['id']
    render=read(root/'render/analysis.json');hist=arrays(root/'render/portrait-histograms.npz')['histograms'];assert hist.shape==(2,850,1300)
    assert np.array_equal(hist.sum((1,2)),[8000000,8000000]) and np.array(render['frame_point_counts']).shape==(72,2)
    assert np.all(np.array(render['frame_point_counts'])==1000000)
    for r in render['runs']:
        basis=np.array(r['basis']);assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12 and np.linalg.det(basis)>0 and r['clipped']==0
        z=arrays(root/'render'/f"case-{r['id']}-animation.npz")['states'];assert z.shape==(1000000,3) and np.isfinite(z).all()
    with Image.open(root/'render/yogi-tanh-rotation.gif') as gif:assert gif.n_frames==72 and gif.size==(1400,660)
    a=read(root/'structure/analysis.json');assert a['declared']==a['completed']==1024 and a['adam_epsilon_outside_local_derivative_identity']
    print('PASS: frozen packet; 8256 scouts; 488 followups; 8 long runs; 420 cohorts; 6 precision runs; 10 actual library histories; 36 gallery clouds; exact render accounting; structure controls',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);main(p.parse_args())
