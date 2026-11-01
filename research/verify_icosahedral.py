"""Verify the frozen checkpoint and its declared finite evidence, not priority."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,tempfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/icosahedral-flow'
parser=argparse.ArgumentParser();parser.add_argument('--scratch',type=Path,required=True);args=parser.parse_args()
manifest=json.loads((BASE/'manifest.json').read_text())
for row in manifest['files']+manifest['inputs']:
    f=ROOT/row['path'];assert hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256'],row['path']
d=np.load(BASE/'icosa-v1/input.npz');scout=np.load(BASE/'icosa-v1/data.npz')
positive=np.where((scout['windows'].min(1)>.01)&(scout['speed']>1e-4))[0].tolist()
assert positive==json.loads((BASE/'icosa-v1/summary.json').read_text())['positive_ids']
assert len(positive)==51 and scout['tail'].shape==(256,4000,3)
refine=np.load(BASE/'icosa-spectra/data.npz');long=np.load(BASE/'icosa-long/sample.npz')
for ident in [177,183,43]:
    mask=refine['ids']==ident;assert np.all(refine['windows'][mask,:,0]>.01)
    assert np.max(np.abs(refine['windows'][mask].sum(2)+3*d['gamma'][ident]))<5e-6
    j=int(np.where(long['ids']==ident)[0][0]);assert np.all(long['windows'][j,:,0]>.01)
    assert np.max(np.abs(long['windows'][j].sum(1)+3*d['gamma'][ident]))<2e-8
failed=(refine['ids']==85)&(refine['h']==.00625)
assert np.max(np.abs(refine['windows'][failed]+d['gamma'][85]))<1e-8
for folder in ['icosa-dop853','icosa-dop853-177']:
    report=json.loads((BASE/folder/'summary.json').read_text());w=np.array(report['windows'])
    assert np.all(w[:,0]>.01) and report['trace_error']<1e-10
    assert np.max(np.abs(w.sum(1)+3*d['gamma'][report['id']]))<1e-10
basin=np.load(BASE/'icosa-basins/data.npz');positive=(basin['windows'][:,:,0].min(1)>.01)&(basin['speed']>1e-5);near=np.linalg.norm(basin['end'],axis=1)<1e-6
for row in json.loads((BASE/'icosa-basins/summary.json').read_text())['rows']:
    mask=(basin['ids']==row['id'])&(basin['labels']==row['kind'])
    assert int(mask.sum())==row['count'] and int((mask&positive).sum())==row['positive'] and int((mask&near).sum())==row['origin']
for row in json.loads((BASE/'icosa-final/render.json').read_text()):
    raw=np.load(BASE/f"icosa-final/portrait-{row['id']}.npz");cam=raw['camera']
    assert np.max(np.abs(cam@cam.T-np.eye(3)))<1e-12 and np.linalg.det(cam)>0
    assert int(raw['histogram'].sum())==row['histogram_sum']==row['states']==1280000
    cams=np.array(row['animation']['cameras']);assert np.max(np.abs(cams@cams.transpose(0,2,1)-np.eye(3)))<1e-12
    assert row['animation']['states']==320000 and not row['group_replication']
    from PIL import Image
    with Image.open(BASE/f"icosa-final/icosahedral-{row['id']}.gif") as image:assert image.n_frames==72
reiter=json.loads((BASE/'reiter-replay/checks.json').read_text());assert reiter['checkpoint_error']<6e-7
reference=np.load(BASE/'reiter-replay/trajectory.npz');assert reference['points'].shape==(320000,3) and np.isfinite(reference['points']).all()
for original in (BASE/'source-original').glob('*.py'):
    expected=original.read_text(encoding='utf-8-sig').replace('from pathlib import Path','from pathlib import Path\nimport sys\nfrom icosa_paths import ROOT,RESEARCH,OUTPUT\nsys.path.insert(0,str(RESEARCH/\'unexpected_search\'))',1)
    expected=expected.replace('root=Path(__file__).parent','root=OUTPUT').replace('base=Path(__file__).parent','base=OUTPUT').replace("repo=Path('D:/Codex/chaos-attractor-lab')",'repo=ROOT').replace("sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')",'sys.path.insert(0,str(RESEARCH))')
    assert expected==(ROOT/'research/icosahedral_search'/original.name).read_text(),original.name
cycle=json.loads((BASE/'icosa-cycle/report.json').read_text())
assert cycle['shooting_residual']<1e-10 and max(cycle['nontrivial_multiplier_moduli'])<1
assert abs(cycle['log_abs_determinant_over_period']-cycle['trace_target'])<1e-10
# Re-evaluate the interval inequalities rather than trusting saved success flags.
args.scratch.mkdir(parents=True,exist_ok=True)
with tempfile.TemporaryDirectory(prefix='icosa-verify-',dir=args.scratch) as tmp:
    temp=Path(tmp);(temp/'icosa-v1').mkdir();import shutil,os
    shutil.copy2(BASE/'icosa-v1/input.npz',temp/'icosa-v1/input.npz')
    env=dict(os.environ,CHAOS_ICOSA_OUTPUT=str(temp),PYTHONDONTWRITEBYTECODE='1')
    result=subprocess.run([sys.executable,str(ROOT/'research/icosahedral_search/icosa_absorbing.py')],env=env,capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    actual=json.loads((temp/'icosa-absorbing/certificate.json').read_text());saved=json.loads((BASE/'icosa-absorbing/certificate.json').read_text())
    assert actual==saved
print(json.dumps({'verified_files':len(manifest['files']),'verified_inputs':len(manifest['inputs']),'scout_positives':51,'long_survivors':[177,183,43],'rejected_long_candidate':85,'interval_absorbing_bounds':'re-evaluated','historical_novelty':'not certified'}))
