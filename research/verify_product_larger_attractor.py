"""Replay and independent-expression checks for the larger-attractor probes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from probe_product_core_recurrence import run,value
from search_product_alphabet import step
from search_soft_sphere import independent


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name in [Path(__file__).name,'probe_product_core_recurrence.py','search_product_alphabet.py',
                 'search_soft_sphere.py','follow_spherical_constellations.py','search_soft_constellations.py',
                 'search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    read=lambda p:json.loads(p.read_text());rows=[]
    for source in sorted(args.probe.glob('*.npz')):
        replay=args.replay/source.name
        with np.load(source) as a,np.load(replay) as b:
            assert a.files==b.files
            for key in a.files:assert a[key].dtype==b[key].dtype and a[key].shape==b[key].shape and a[key].tobytes()==b[key].tobytes()
            rows.append(dict(file=source.name,arrays=len(a.files),sha256=sha(source),same_environment_replay=True))
    a=read(args.probe/'report.json');b=read(args.replay/'report.json');a.pop('seconds');b.pop('seconds');assert a==b
    source=np.load(args.source);beta,tau=source['params'];levels=source['levels'];points=source['points']
    r=read(args.recurrence/'report.json');data=np.load(args.recurrence/'recurrence.npz');targets=data['targets'];checks=[]
    rng=np.random.default_rng(2026110161);max_error=0.
    probes=np.vstack([rng.uniform(-.5,.5,(128,3)),data['closest_states_0'],data['closest_states_1']])
    for x in probes:
        y,_=value(x,beta,tau,levels);z=independent(x,beta,points,tau);s=step(x,beta,tau,levels,2)[0]
        max_error=max(max_error,float(np.max(abs(y-z))),float(np.max(abs(s-z))))
    assert max_error<1e-12
    from fractions import Fraction
    raw=read(args.input);inverse=np.array([[float(Fraction(v)) for v in row] for row in raw['exact_windows'][0]['inverse']])
    for k,start in enumerate(data['starts']):
        h=data[f'history_{k}'];counts=data[f'counts_{k}'];closest=data[f'closest_states_{k}']
        assert np.all(np.diff(h,axis=0)<=0) and np.all(np.diff(counts,axis=1)<=0)
        assert np.allclose(np.linalg.norm(closest-targets,axis=1),h[-1],rtol=2e-15,atol=1e-18)
        replay=run(start,beta,tau,levels,targets,inverse,r['warm'],r['block'],r['block'],data['radii'])
        assert np.array_equal(replay[1][0],h[0]) and np.array_equal(replay[4][0],counts[0])
        assert h[-1].tolist()==r['starts'][k]['closest_distances'] and counts.sum(0).tolist()==r['starts'][k]['counts_by_radius']
        checks.append(dict(start=k,replayed_steps=r['block'],first_block_distances_and_counts_exact=True,closest_state_distances_checked=True))
    result=dict(status='passed',probe_replays=rows,report_equal_except_runtime=True,recurrence_checks=checks,
        independent_formula_probes=len(probes),maximum_formula_error=max_error,
        scope='All branch/perturbation arrays replay bit for bit. Recurrence first blocks replay; all saved closest-state distances and cumulative counts are checked. Full-site nested reflections independently check floating formulas. These tests are not an infinite-time theorem.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','source','probe','replay','recurrence','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
