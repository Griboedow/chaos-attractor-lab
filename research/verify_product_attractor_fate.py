"""Check finite fate records, replay a declared subset, and check full-site tangents."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from numba import set_num_threads
from search_product_attractor_fate import fast_step,follow


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name in [Path(__file__).name,'search_product_attractor_fate.py','verify_product_coverings.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    read=lambda p:json.loads(p.read_text());protocol=read(args.run/'protocol.json')
    assert sha(args.source)==protocol['source_sha256']
    with np.load(args.source) as f:beta,tau=f['params'];levels=f['levels']
    with np.load(args.run/'fate.npz') as f:data={k:f[k] for k in f.files}
    assert np.array_equal(data['starts'],protocol['starts']) and not np.any(data['flags'])
    keys=['tail','rates','minima','minimum_states','minimum_times','ends','flags','bounds']
    for offset in range(0,256,16):
        with np.load(args.run/f'block-{offset:03}.npz') as f:
            for key in ['starts']+keys:assert np.array_equal(data[key][offset:offset+16],f[key])
    selected=np.array([64*g+i for g in range(4) for i in [0,16,32,48]])
    replay=follow(data['starts'][selected],beta,tau,levels,protocol['warm'],protocol['count'],protocol['block'])
    for key,values in zip(keys,replay):assert np.array_equal(values,data[key][selected]),key
    periods=np.zeros(256,int)
    for k in range(256):
        for lag in list(range(1,129))+[256,512]:
            if np.max(np.linalg.norm(data['tail'][k,-512:]-data['tail'][k,-512-lag:-lag],axis=1))<1e-9:periods[k]=lag;break
    assert np.array_equal(periods,data['periods']) and not periods.any()
    assert data['rates'].shape==(256,20) and data['rates'].min()>0
    for x,minimum in zip(data['minimum_states'],data['minima'].min(1)):
        assert fast_step(x,np.array([1.,0.,0.]),beta,tau,levels)[2]==minimum
    sys.path.append(str(args.mpmath_path));import mpmath as mp
    import verify_product_coverings as full
    mp.mp.dps=mp.iv.dps=90;full.mp=mp;independent=full.Model(read(args.input))
    rng=np.random.default_rng(2026110203)
    probes=np.vstack([data['starts'],data['minimum_states'],rng.uniform(-2.2,2.2,(128,3))])
    errors=np.zeros(4)
    for x in probes:
        w=rng.normal(size=3);w/=np.linalg.norm(w)
        y,z,r,g,flag=fast_step(x,w,beta,tau,levels);assert flag==0
        f,J,radius=independent.value_jacobian([mp.iv.mpf(float(v)) for v in x])
        f=np.array([float(v.mid) for v in f]);J=np.array([[float(v.mid) for v in row] for row in J]);jw=J@w;length=np.linalg.norm(jw)
        errors=np.maximum(errors,[np.max(abs(y-f)),np.max(abs(z-jw/length)),abs(r-float(radius.mid)),abs(g-np.log(length))])
    assert errors.max()<1e-10
    result=dict(status='passed',source_sha256=sha(args.source),fate_sha256=sha(args.run/'fate.npz'),
        all_16_batches_match_combined_archive=True,all_period_tests_recomputed=True,all_minimum_states_checked=True,
        replay_indices=selected.tolist(),observed_steps_replayed=int(len(selected)*protocol['count']),all_8_replayed_arrays_exact=True,
        full_site_probes=len(probes),full_site_digits=90,maximum_errors=dict(zip(['value','normalized_tangent','normalization','log_growth'],errors.tolist())),
        minimum_leading_block=float(data['rates'].min()),minimum_normalization=float(data['minima'].min()),
        scope='Finite numerical validation only. Sixteen of 256 trajectories replayed in full; all batch concatenations, period tests and closest-normalization states checked. Independent full-site expression and AD derivatives checked at 640 states; no asymptotic or basin-volume conclusion.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','source','run','mpmath-path','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
