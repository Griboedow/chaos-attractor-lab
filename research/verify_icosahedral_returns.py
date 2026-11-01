"""Verify finite return evidence and repeat independent numerical controls."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/icosahedral-returns'
p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path,required=True);args=p.parse_args()
manifest=json.loads((BASE/'manifest.json').read_text())
for row in manifest['files']+manifest['inputs']:
    file=ROOT/row['path'];assert file.stat().st_size==row['bytes'],row['path']
    assert hashlib.sha256(file.read_bytes()).hexdigest()==row['sha256'],row['path']
for original in (BASE/'source-original').glob('*.py'):
    expected=original.read_text(encoding='utf-8-sig').replace('from pathlib import Path','from pathlib import Path\nfrom return_paths import INPUT,STRUCTURE,CYCLES,OUTPUT',1).replace('base=Path(__file__).parent','base=OUTPUT')
    for stage,target in [('icosa-v1','INPUT'),('icosa-itinerary-final','STRUCTURE'),('icosa-cycle-refinement','CYCLES'),('icosa-directed-cycles','CYCLES')]:
        expected=expected.replace(f"base/'{stage}/",f"{target}/'{stage}/")
    expected=expected.replace("np.load(base/'icosa-long/177-dense.npy',mmap_mode='r')[::128]","np.load(INPUT/'icosa-long/sample.npz')['tail'][list(np.load(INPUT/'icosa-long/sample.npz')['ids']).index(177)][::8]")
    assert expected==(ROOT/'research/icosahedral_returns'/original.name).read_text(),original.name

scout=np.load(BASE/'icosa-return-scout/data.npz');assert scout['starts'].shape==(513,3)
assert int(scout['event_counts'].sum())==1972
B=scout['basis'];normal=scout['normal'];inverse=scout['inverse_basis']
assert np.max(abs(normal@B))<1e-14 and np.max(abs(inverse@B-np.eye(2)))<1e-14
grid=np.load(BASE/'icosa-return-grid/data.npz');x=grid['inputs'];y=grid['outputs'];J=grid['derivatives'];strips=grid['strips'];V=float(grid['stable_half_width'])
assert x.shape==y.shape==(2,17,17,2) and np.isfinite(y).all() and np.isfinite(J).all()
assert np.array_equal(strips,np.array([[-8e-6,2e-6],[.002105,.00212]])) and V==.0002
report=json.loads((BASE/'icosa-return-grid/report.json').read_text())
for k,row in enumerate(report['rows']):
    assert np.array_equal(x[k,:,:,0],np.broadcast_to(np.linspace(*strips[k],17)[:,None],(17,17)))
    assert np.array_equal(x[k,:,:,1],np.broadcast_to(np.linspace(-V,V,17)[None,:],(17,17)))
    expansion=abs(J[k,:,:,0,0])-.1*abs(J[k,:,:,0,1])
    slope=(abs(J[k,:,:,1,0])+.1*abs(J[k,:,:,1,1]))/expansion
    recomputed={'left_margin_beyond_both_targets':y[k,0,:,0].min()-strips[:,1].max(),'right_margin_beyond_both_targets':strips[:,0].min()-y[k,-1,:,0].max(),'stable_margin':V-abs(y[k,:,:,1]).max(),'sampled_minimum_unstable_expansion':expansion.min(),'sampled_maximum_image_cone_slope':slope.max(),'minimum_crossing_speed':grid['crossing_speeds'][k].min()}
    for key,value in recomputed.items():assert np.isclose(value,row[key],rtol=1e-13,atol=0),key
    assert min(recomputed['left_margin_beyond_both_targets'],recomputed['right_margin_beyond_both_targets'],recomputed['stable_margin'])>0
    assert expansion.min()>447 and slope.max()<.025
controls=json.loads((BASE/'icosa-return-controls/report.json').read_text())
assert controls['implementation_max_error']<2e-15 and len(controls['rows'])==10
for row in controls['rows']:
    e=[v['output_difference'] for v in row['fixed_steps']]
    assert all(a/b>16.29 for a,b in zip(e[:-1],e[1:])) and e[-1]<8.5e-11
    assert row['adaptive_output_difference']<1.5e-12 and row['adaptive_derivative_difference']<2.3e-8
assert max(r['max_absolute_error'] for r in controls['derivative_controls'])<1.91e-5
with Image.open(ROOT/'docs/assets/icosahedral-returns.png') as im:assert im.size==(2520,1620)

args.scratch.mkdir(parents=True,exist_ok=True)
with tempfile.TemporaryDirectory(prefix='icosa-returns-',dir=args.scratch) as tmp:
    output=Path(tmp);shutil.copytree(BASE/'icosa-return-scout',output/'icosa-return-scout')
    # Pin all inputs to the manifest-checked packet, ignoring external overrides.
    os.environ.update(CHAOS_ICOSA_INPUT=str(ROOT/'results/icosahedral-flow'),CHAOS_ICOSA_STRUCTURE_INPUT=str(ROOT/'results/icosahedral-structure'),CHAOS_ICOSA_CYCLES_INPUT=str(ROOT/'results/icosahedral-cycles'),CHAOS_ICOSA_RETURNS_OUTPUT=tmp,PYTHONDONTWRITEBYTECODE='1')
    sys.path.insert(0,str(ROOT/'research/icosahedral_returns'))
    import icosa_poincare as flow
    from icosa_fast_return import rhs,rk_return
    # Arrays are materialized; release ZIP handles before Windows temp cleanup.
    flow.d.close();flow.frame.close()
    random=np.random.default_rng(2026110157);implementation_error=0.
    for _ in range(30):
        state=np.r_[random.normal(size=3),random.normal(size=6)]
        implementation_error=max(implementation_error,float(np.max(abs(rhs(0,state)-flow.field(0,state)))))
    assert implementation_error<1e-12
    fresh=[]
    for k,i,j in [(0,0,0),(0,16,16),(0,8,8),(1,0,0),(1,16,16),(1,8,8)]:
        uv=x[k,i,j];q,t,derivative,speed=flow.evaluate(uv,rtol=3e-14,atol=3e-16,max_step=.015)
        discrepancy=float(np.linalg.norm(q-y[k,i,j]));assert discrepancy<1e-9
        assert abs(t-grid['times'][k,i,j])<1e-8 and np.max(abs(derivative-J[k,i,j]))<1e-5
        fixed,when,crossing=rk_return(uv,.00125)
        assert np.linalg.norm(fixed-q)<1e-9 and abs(when-t)<1e-8
        fresh.append(discrepancy)
        print(json.dumps({'checked_input':uv.tolist(),'adaptive_discrepancy':discrepancy}),flush=True)
    run=subprocess.run([sys.executable,str(ROOT/'research/icosahedral_returns/icosa_spherical_prior.py')],env=os.environ,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    actual=json.loads((output/'icosa-spherical-prior/report.json').read_text());saved=json.loads((BASE/'icosa-spherical-prior/report.json').read_text())
    assert actual==saved
    assert actual['edge_count']==60 and abs(actual['maximum_direction_dot']+1)<1e-14
    assert actual['our_sampled_radius_range'][1]-actual['our_sampled_radius_range'][0]>.4
    assert actual['our_normalized_sextic_range'][1]-actual['our_normalized_sextic_range'][0]>.5
print(json.dumps({'files':len(manifest['files']),'inputs':len(manifest['inputs']),'grid_starts':578,'all_sampled_margins_positive':True,'fresh_adaptive_max_discrepancy':max(fresh),'known_graph_edges':60,'covering_chaos_priority':'not certified'}))
