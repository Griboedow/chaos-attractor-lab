"""Verify the frozen membership packet, dependency audits and numerical figure."""
import argparse
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import numpy as np
from audit_raar_coverings import bounds,sha
from audit_raar_critical_membership import main as audit
from propose_raar_critical_crossing import follow
from certify_raar_fold import jet as analytic_jet
from check_raar_critical_fold_ad import jet as automatic_jet
from certify_raar_cycles import exactfloat


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);repo=Path(__file__).resolve().parents[1];packet=repo/'results/raar-critical-membership';base=packet/'runs'
    manifest=json.loads((packet/'manifest.json').read_text())
    paths={r['path'] for r in manifest['files']}
    assert paths=={p.relative_to(packet).as_posix() for p in packet.rglob('*') if p.is_file() and p.name!='manifest.json'}
    for row in manifest['files']:
        p=packet/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes']
    for name,old in [('raar-critical-thin-fold-region01','raar-critical-membership-audit03'),('raar-critical-thin-fold-ad01','raar-critical-membership-audit04')]:
        out=args.out/name;audit(SimpleNamespace(base=base,out=out,fold=name))
        assert (out/'report.json').read_bytes()==(base/old/'report.json').read_bytes()
    figure=base/'raar-critical-membership-render02';analysis=json.loads((figure/'analysis.json').read_text());data=np.load(figure/'figure-data.npz')
    raw=json.loads((base/'raar-critical-thin-window01/input.json').read_text());S=np.array(raw['bases'][0])*np.array(raw['radii'][0]);center=np.array(raw['centers'][0]);p=np.array(raw['params'])
    assert np.array_equal(p,data['params']) and np.array_equal(data['normalized_u'],np.linspace(-1,1,8193))
    starts=center+data['normalized_u'][:,None]*S[:,0];assert np.array_equal(starts,data['starts'])
    saved,dets,axes=follow(starts,p,42)
    assert np.array_equal(saved[-1],data['states42']) and np.array_equal(dets[-1],data['determinants']) and np.array_equal(axes.min(0),data['axis_minima'])
    assert float(axes.min())==analysis['minimum_sampled_axis']
    with gzip.open(base/'raar-critical-thin-fold-region01/volume.jsonl.gz','rt') as stream:critical=[r for r in map(json.loads,stream) if r['classification']=='fold_if_critical']
    low=[min(bounds(v['image'][i])[0] for v in critical) for i in range(3)];high=[max(bounds(v['image'][i])[1] for v in critical) for i in range(3)]
    assert analysis['critical_box_exact']==[[str(v) for v in low],[str(v) for v in high]]
    box=np.array([[np.nextafter(float(v),-np.inf) for v in low],[np.nextafter(float(v),np.inf) for v in high]])
    assert np.array_equal(box,data['critical_box'])
    intervals=[]
    for record in critical:
        m=np.zeros(3);r=np.ones(3)
        for k in range(0,len(record['path']),2):
            j=int(record['path'][k]);r[j]*=.5;m[j]+=r[j]*(1 if record['path'][k+1]=='1' else -1)
        intervals.append([m[0]-r[0],m[0]+r[0]])
    assert np.array_equal(np.array(intervals),data['critical_u_intervals'])
    assert np.max(abs(data['camera_rotation'].T@data['camera_rotation']-np.eye(2)))<1e-12
    assert sha(repo/'results/raar-circles/render/case-1051-animation.npz')==analysis['reference_sha256']
    assert sha(repo/'docs/assets/raar-critical-membership.png')==sha(figure/'raar-critical-membership.png')
    # Supplemental formula control, not a replacement for the full interval runs.
    mp.mp.dps=100;mp.iv.dps=100;iv=mp.iv;pi=[exactfloat(v,iv) for v in p];components=0
    for state in data['states42'][::128]:
        x=[exactfloat(v,iv) for v in state];a=analytic_jet(x,pi,iv);b=automatic_jet(x,pi,iv)
        av=[a[0]]+a[1]+sum(a[2],[])+sum([sum(v,[]) for v in a[3]],[])+a[4]
        bv=[b[0]]+b[1]+sum(b[2],[])+sum([sum(v,[]) for v in b[3]],[])+b[4]
        for left,right in zip(av,bv):assert max(left.a,right.a)<=min(left.b,right.b);components+=1
    result={'status':'passed','manifest_sha256':sha(packet/'manifest.json'),'files':len(paths),'primary_and_secondary_accounting_reports':'byte-identical replay','guide_states':len(starts),'steps':42,'figure_arrays':'bitwise replay','analytic_automatic_jet_component_intersections':components,'scope':'Frozen-byte verification, complete accounting replay, figure-state replay and supplemental derivative controls. Full interval enclosure regeneration remains the separate published commands; the mathematical continuum argument is stated explicitly.'}
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args())
