"""Recheck finite graph evidence and independent numerical periodic controls."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree,distance
from PIL import Image

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/icosahedral-cycles'
p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path,required=True);args=p.parse_args()
manifest=json.loads((BASE/'manifest.json').read_text())
for row in manifest['files']+manifest['inputs']:
    assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
for original in (BASE/'source-original').glob('*.py'):
    expected=original.read_text(encoding='utf-8-sig').replace('from pathlib import Path','from pathlib import Path\nimport sys\nfrom cycle_paths import ROOT,RESEARCH,INPUT,STRUCTURE,OUTPUT\nsys.path.insert(0,str(RESEARCH/\'unexpected_search\'))',1).replace('base=Path(__file__).parent','base=OUTPUT')
    for key in ['icosa-v1','icosa-long']:
        expected=expected.replace(f"base/'{key}/",f"INPUT/'{key}/").replace(f"base/f'{key}/",f"INPUT/f'{key}/")
    for key in ['icosa-itinerary-final','icosa-graph-geometry']:
        expected=expected.replace(f"base/'{key}/",f"STRUCTURE/'{key}/")
    expected=expected.replace("sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')",'sys.path.insert(0,str(RESEARCH))').replace("runpy.run_path(str(base/'icosa_taylor_precision.py'))","runpy.run_path(str(Path(__file__).parent/'icosa_taylor_precision.py'))")
    assert expected==(ROOT/'research/icosahedral_cycles'/original.name).read_text(),original.name
rec=json.loads((BASE/'icosa-recurrences/candidates.json').read_text());assert len(rec['selected'])==40 and rec['raw_qualifying_pairs']==4769
attempts=json.loads((BASE/'icosa-relative-cycles/report.json').read_text())['rows'];assert len(attempts)==24
assert [r['index'] for r in attempts if r['accepted_relative_root']]==[0,3,10,18,19,20]
structure=ROOT/'results/icosahedral-structure';foci=np.load(structure/'icosa-itinerary-final/data.npz')['foci'];geo=np.load(structure/'icosa-graph-geometry/geometry.npz');data=np.load(BASE/'icosa-directed-cycles/data.npz')
assert distance.pdist(foci).min()>.6
def count(ev):
    C=np.zeros((30,30),int);np.add.at(C,(ev[:-1,1],ev[1:,1]),1);return C
wide=data['events_0.2'];C=count(wide);A=C>0
assert np.array_equal(C,data['counts']) and len(wide)==1992 and A.sum()==60
assert np.array_equal(A.astype(int)+A.T.astype(int),geo['edge_mask'].astype(int))
assert np.all(A.sum(0)==2) and np.all(A.sum(1)==2)
assert np.array_equal(data['events_0.16'][:,1],wide[:,1])
assert np.array_equal(data['events_0.25'],data['events_0.3']) is False  # Entry times differ.
assert np.array_equal(data['events_0.25'][:,1],data['events_0.3'][:,1])
assert np.array_equal(data['events_0.25'][:-1,1],wide[:,1])
for name in ['triangles','pentagons']:
    for face in geo[name]:
        sub=A[np.ix_(face,face)];assert np.all(sub.sum(0)==1) and np.all(sub.sum(1)==1)
narrow=data['events_0.12'];skipped=0
for (t0,a),(t1,b) in zip(narrow[:-1],narrow[1:]):
    if A[a,b]:continue
    left=np.searchsorted(wide[:,0],t0,side='right')-1;right=np.searchsorted(wide[:,0],t1,side='right')-1;path=wide[left:right+1,1]
    assert len(path)==3 and path[0]==a and path[-1]==b and np.all(A[path[:-1],path[1:]])
    skipped+=1
assert skipped==194 and (count(narrow)>0).sum()==118
co=np.load(ROOT/'results/icosahedral-flow/icosa-v1/input.npz');c=co['centers'][177];u=co['axes'][177];s=co['sigma'][177];g=co['gamma'][177]
def field(t,x):
    r=x-c;return -g*x+(np.exp(-.5*np.sum(r*r,axis=1)/s**2)[:,None]*np.cross(u,r)).sum(0)
rows=json.loads((BASE/'icosa-cycle-refinement/report.json').read_text())['rows'];closures=[]
for row in rows:
    orbit=np.load(BASE/f"icosa-cycle-refinement/cycle-{row['index']}.npz");points=orbit['points'];assert points.shape==(40000,3)
    dist,idx=cKDTree(foci).query(points);ids=idx[dist<.2];seq=ids[np.r_[True,ids[1:]!=ids[:-1]]]
    if seq[0]==seq[-1]:seq=seq[:-1]
    assert len(seq)==row['segments'] and np.all(A[seq,np.roll(seq,-1)])
    assert any(set(seq)==set(f) for f in geo['triangles' if len(seq)==3 else 'pentagons'])
    errors=[r['closure_error'] for r in row['rk4']];assert all(a>10*b for a,b in zip(errors[:-1],errors[1:]))
    assert abs(row['unstable_multiplier'])>400 and row['shooting_residual']<1e-10
    assert np.isclose(row['stable_modulus_inferred_from_divergence'],np.exp(-3*g*row['period'])/abs(row['unstable_multiplier']),rtol=1e-12,atol=0)
    times=np.r_[np.arange(0,40000,400)*row['period']/40000,row['period']]
    sol=solve_ivp(field,[0,row['period']],row['point'],method='DOP853',rtol=2e-13,atol=2e-15,max_step=.005,t_eval=times)
    assert sol.success and np.max(np.linalg.norm(sol.y[:,:-1].T-points[::400],axis=1))<1e-7
    error=float(np.linalg.norm(sol.y[:,-1]-row['point']));assert error<1e-7;closures.append(error)
    print(json.dumps({'checked_cycle':row['index'],'adaptive_closure_error':error}),flush=True)
render=json.loads((BASE/'icosa-cycle-visuals-final/render.json').read_text())
for panel in render['panels']:
    cam=np.array(panel['camera']);assert np.max(abs(cam@cam.T-np.eye(3)))<1e-12 and np.linalg.det(cam)>.999999
with Image.open(ROOT/'docs/assets/icosahedral-cycles.png') as image:assert image.width>2000 and image.height>1000
args.scratch.mkdir(parents=True,exist_ok=True)
with tempfile.TemporaryDirectory(prefix='icosa-cycles-',dir=args.scratch) as tmp:
    out=Path(tmp);shutil.copytree(BASE/'icosa-cycle-refinement',out/'icosa-cycle-refinement');shutil.copytree(BASE/'icosa-taylor-precision',out/'icosa-taylor-precision')
    env=dict(os.environ,CHAOS_ICOSA_INPUT=str(ROOT/'results/icosahedral-flow'),CHAOS_ICOSA_CYCLES_OUTPUT=str(out),PYTHONDONTWRITEBYTECODE='1')
    run=subprocess.run([sys.executable,str(ROOT/'research/icosahedral_cycles/icosa_taylor_checks.py')],env=env,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    actual=json.loads((out/'icosa-taylor-checks/report.json').read_text());saved=json.loads((BASE/'icosa-taylor-checks/report.json').read_text());assert actual==saved
print(json.dumps({'files':len(manifest['files']),'inputs':len(manifest['inputs']),'missed_intermediate_visits':skipped,'directed_edges':60,'reintegrated_closure_errors':closures,'periodic_existence_chaos_priority':'not certified'}))
