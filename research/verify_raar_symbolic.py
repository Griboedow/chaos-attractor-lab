"""Audit both complete covering outputs, the located cycle and figure provenance."""
import argparse,json,tempfile
from decimal import Decimal
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.spatial import cKDTree
import audit_raar_coverings
from archive_raar_symbolic import sha
from verify_raar_cycles import inclusion,interval

def main(args):
    repo=Path(__file__).resolve().parents[1];root=repo/'results/raar-symbolic-dynamics';runs=root/'runs';js=lambda p:json.loads(p.read_text(encoding='utf-8'))
    manifest=js(root/'manifest.json')['files'];actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts};assert actual=={r['path'] for r in manifest}
    for r in manifest:p=root/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],p
    for r in js(root/'provenance.json'):assert sha(root/r['path'])==r['sha256']
    raw=js(runs/'raar-covering-refined01/input.json');assert sha(runs/'raar-covering-refined01/input.json')=='d6330277f1f10f5c7bc2118acb69dba5a10fbd025761c8f5b83beb5344b13b60'
    with tempfile.TemporaryDirectory(prefix='raar-symbolic-audit-',dir=args.scratch) as td:
        for mode,total in [('primary',16467),('ad',9449)]:
            out=Path(td)/mode;audit_raar_coverings.main(SimpleNamespace(base=runs,mode=mode,edges=list(range(8)),out=out));current=js(out/'report.json');assert current==js(runs/f'raar-covering-audit-{mode}01/report.json');assert current['total_accepted_cells']==total and current['complete_graph'] and Decimal(current['global_minimum_axis'])>Decimal('.10367')
    cyclefile=runs/'raar-symbolic-cycle-proof01/symbolic-period-49.json';cycle=js(cyclefile);inclusion(cycle);assert cycle['period']==49 and cycle['parameters']==raw['params'];assert [r['divisor'] for r in cycle['proper_divisor_exclusions']]==[1,7]
    for r in cycle['proper_divisor_exclusions']:assert any(a>0 or b<0 for a,b in map(interval,r['coordinate_difference']))
    assert all(interval(v)[0]>0 for v in cycle['axis_distances']);assert interval(cycle['characteristic_at_plus_one'])[1]<0
    proof=js(runs/'raar-symbolic-cycle-proof01/report.json');assert proof['status']=='passed' and sha(runs/'raar-symbolic-cycle-proof01/certify_raar_cycles.py')==proof['source_sha256'] and sha(runs/'raar-symbolic-cycle-proof01/input.json')==proof['input_sha256'] and proof['cases']==[cycle]
    audit=js(runs/'raar-symbolic-cycle-audit01/report.json');assert audit['status']=='passed' and audit['cycle_certificate_sha256']==sha(cyclefile) and audit['windows_sha256']==sha(runs/'raar-covering-refined01/input.json');assert audit['itinerary']==list(range(7)) and audit['unstable_dimensions']==1 and audit['stable_dimensions']==2
    assert all(-1<a<=b<1 for row in audit['coordinate_boxes'] for a,b in map(interval,row));assert all(interval(v)[0]>0 for v in audit['remaining_quadratic_schur_margins']);a,b=interval(audit['unstable_multiplier_enclosure']);assert Decimal('579.0978')<a<b<Decimal('579.0982')
    numerical=js(runs/'raar-symbolic-orbits01/report.json');assert len(numerical['cases'])==7
    for row in numerical['cases']:
        assert row['solver_success'] and row['accepted_numerical'] and row['maximum_physical_residual']<1.8e-14 and row['maximum_window_coordinate']<.938
        data=np.load(runs/'raar-symbolic-orbits01'/f"word-{row['word']}.npz");assert len(data['orbit'])==row['steps']==49*len(row['word']) and np.array_equal(data['params'],raw['params']);assert data['residuals'].max()==row['maximum_physical_residual'] and abs(data['coordinates']).max()==row['maximum_window_coordinate']
    assert js(runs/'raar-covering-ad-edge0-01/interrupted.json')['status']=='interrupted_incomplete'
    assert len(js(runs/'raar-covering-proposal01/report.json')['cases'])==36 and len(js(runs/'raar-return-search01/report.json')['attempts'])==2474
    final=runs/'raar-symbolic-render02';fig=np.load(final/'figure-data.npz');view=js(final/'analysis.json');reference=repo/'results/raar-circles/render/case-1051-animation.npz';assert sha(reference)==view['reference_sha256'];z=np.load(reference)['states'];xy=(z-fig['center'])@fig['rotation'];h,_,_=np.histogram2d(xy[:,0],xy[:,1],bins=[fig['xedges'],fig['yedges']]);assert np.array_equal(h,fig['histogram']) and int(h.sum())==view['histogram_count']==view['reference_points']==1000000
    assert fig['rotation'].shape==(3,2) and np.max(abs(fig['rotation'].T@fig['rotation']-np.eye(2)))<1e-12;assert np.array_equal(fig['orbit49'],np.array(audit['orbit_decimal'],float)) and len(fig['orbit49'])==49 and len(fig['orbit7'])==7
    for orbit in [fig['orbit49'],fig['orbit7']]:
        q=(orbit-fig['center'])@fig['rotation'];assert np.all(q[:,:2]>[fig['xedges'][0],fig['yedges'][0]]) and np.all(q[:,:2]<[fig['xedges'][-1],fig['yedges'][-1]])
    tree=cKDTree(fig['orbit7'])
    for word in view['words']:
        orbit=np.load(runs/'raar-symbolic-orbits01'/f'word-{word}.npz')['orbit'];distance=tree.query(orbit)[0];assert np.array_equal(distance,fig['word_'+word]);assert distance.min()>1e-4 and distance.max()<.7
    initial=np.load(runs/'raar-symbolic-render01/figure-data.npz');assert initial.files==fig.files
    for key in fig.files:assert np.array_equal(initial[key],fig[key])
    for name in ['geometry','graph','itineraries']:assert sha(final/f'raar-symbolic-{name}.png')==sha(repo/'docs/assets'/f'raar-symbolic-{name}.png')
    result=dict(status='passed',files=len(manifest),manifest_sha256=sha(root/'manifest.json'),primary_cells=16467,second_cells=9449,located_least_period=49,reference_states=1000000,numerical_words=7,scope='Complete output/partition/hash audit and exact figure-data replay. Interval enclosures are supplied by the two separately executed engines; this audit is not a third enclosure implementation.')
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path);p.add_argument('--output',type=Path);main(p.parse_args())
