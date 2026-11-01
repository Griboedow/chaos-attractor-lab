"""Verify the follow-up packet and repeat all local interval root checks."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile
import numpy as np
from scipy.spatial import ConvexHull
from PIL import Image

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/icosahedral-structure'
p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path,required=True);args=p.parse_args()
manifest=json.loads((BASE/'manifest.json').read_text())
for row in manifest['files']+manifest['inputs']:
    assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
for original in (BASE/'source-original').glob('*.py'):
    expected=original.read_text(encoding='utf-8-sig').replace('from pathlib import Path','from pathlib import Path\nimport sys\nfrom structure_paths import ROOT,RESEARCH,INPUT,OUTPUT\nsys.path.insert(0,str(RESEARCH/\'unexpected_search\'))',1).replace('base=Path(__file__).parent','base=OUTPUT')
    for key in ['icosa-v1','icosa-long','icosa-final']:
        expected=expected.replace(f"base/'{key}/",f"INPUT/'{key}/").replace(f"base/f'{key}/",f"INPUT/f'{key}/")
    expected=expected.replace("sys.path.insert(0,'D:/Codex/chaos-attractor-lab/research')",'sys.path.insert(0,str(RESEARCH))')
    assert expected==(ROOT/'research/icosahedral_structure'/original.name).read_text(),original.name
for precision in [60,90]:
    c=json.loads((BASE/f'icosa-equilibria/certificate-{precision}.json').read_text())
    assert [v['certified_distinct_roots'] for v in c['cases']]==[161,101]
    for case,counts in zip(c['cases'],[(1,80,80),(1,50,50)]):
        assert tuple(sum(r['positive_real_eigenvalue_count']==i for r in case['roots']) for i in range(3))==counts
it=np.load(BASE/'icosa-itinerary-final/data.npz');geo=np.load(BASE/'icosa-graph-geometry/geometry.npz');events=it['events'];counts=np.zeros((30,30),int)
for a,b in zip(events[:-1,1],events[1:,1]):counts[a,b]+=1
assert np.array_equal(counts,it['transition_counts']) and counts.sum()==1796 and (counts>0).sum()==118
assert len(np.unique(events[:,1]))==30 and np.all(counts[~geo['edge_mask']]==0)
v=geo['directions'];edge=np.abs(v@v.T-np.cos(np.pi/5))<1e-8;np.fill_diagonal(edge,False)
assert np.array_equal(edge,geo['edge_mask']) and edge.sum()==120 and np.all(edge.sum(1)==4)
faces={tuple(np.flatnonzero(np.abs(v@e[:3]+e[3])<1e-8)) for e in ConvexHull(v).equations}
assert sum(len(f)==3 for f in faces)==20 and sum(len(f)==5 for f in faces)==12
runs=np.load(BASE/'icosa-saddle-starts/data.npz');assert runs['windows'].shape==(12,4,3) and runs['windows'][:,:,0].min()>.29
gamma=np.load(ROOT/'results/icosahedral-flow/icosa-v1/input.npz')['gamma'][177]
assert np.max(np.abs(runs['windows'].sum(2)+3*gamma))<2e-8
render=json.loads((BASE/'icosa-mechanism-visuals/render.json').read_text());cam=np.array(render['camera']);assert np.max(np.abs(cam@cam.T-np.eye(3)))<1e-12
with Image.open(BASE/'icosa-mechanism-visuals/icosahedral-motion.gif') as image:assert image.n_frames==160
args.scratch.mkdir(parents=True,exist_ok=True)
with tempfile.TemporaryDirectory(prefix='icosa-roots-',dir=args.scratch) as tmp:
    out=Path(tmp);(out/'icosa-equilibria').mkdir();shutil.copy2(BASE/'icosa-equilibria/inventory.json',out/'icosa-equilibria/inventory.json')
    env=dict(os.environ,CHAOS_ICOSA_INPUT=str(ROOT/'results/icosahedral-flow'),CHAOS_ICOSA_STRUCTURE_OUTPUT=str(out),PYTHONDONTWRITEBYTECODE='1')
    result=subprocess.run([sys.executable,str(ROOT/'research/icosahedral_structure/certify_icosa_equilibria.py'),'--precision','60'],env=env,capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    actual=json.loads((out/'icosa-equilibria/certificate-60.json').read_text());saved=json.loads((BASE/'icosa-equilibria/certificate-60.json').read_text());assert actual['cases']==saved['cases']
print(json.dumps({'files':len(manifest['files']),'inputs':len(manifest['inputs']),'recertified_roots':262,'observed_transitions':1796,'observed_directed_edges':118,'chaos_and_priority':'not certified'}))
