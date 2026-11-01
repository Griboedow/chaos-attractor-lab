from pathlib import Path
import sys,json,hashlib,ast
import numpy as np
from scipy.stats import qmc
from PIL import Image,ImageDraw
root=Path('D:/Codex/chaos-attractor-lab');scratch=Path(__file__).parent
sys.path.insert(0,str(root/'research'))
from build_attractor_references import simulate

base=root/'results/gaussian-vortex';checks={}
for name,replay in [('long-runs.npz','portable-replay'),('render-run.npz','portable-render'),('render-histogram.npz','portable-render')]:
    a=np.load(base/name);b=np.load(scratch/replay/name)
    assert a.files==b.files
    for key in a.files:np.testing.assert_array_equal(a[key],b[key])
    checks[name]={'all_array_keys_identical':a.files}
for name in ['gaussian-vortex-triskelion.png','gaussian-vortex-triskelion.gif']:
    a=hashlib.sha256((root/'docs/assets'/name).read_bytes()).hexdigest()
    b=hashlib.sha256((scratch/'portable-images'/name).read_bytes()).hexdigest()
    assert a==b
    checks[name]={'identical_sha256':a}
a=json.loads((base/'render.json').read_text());b=json.loads((scratch/'portable-render/render.json').read_text())
assert a==b;checks['render_metadata_identical']=True
u=qmc.Sobol(11,scramble=True,seed=20261108).random_base2(10)
p=np.empty_like(u);p[:,0]=.03+.6*u[:,0];p[:,1:4]=(u[:,1:4]*2-1)*12;p[:,4:7]=.45+u[:,4:7];p[:,7]=.6+1.4*u[:,7];p[:,8:11]=(u[:,8:11]*2-1)*.7
p[:256,2:4]=p[:256,1,None];p[:256,5:7]=p[:256,4,None];p[:256,9:11]=p[:256,8,None]
np.testing.assert_array_equal(p,np.load(base/'search-summary.npz')['parameters'])
checks['sobol_parameter_array_identical']=True
out=scratch/'catalog-replay';(out/'trajectories').mkdir(parents=True,exist_ok=True)
run=simulate('Lorenz',out);assert run['status']=='ok'
np.testing.assert_array_equal(np.load(out/'trajectories/Lorenz.npz')['points'],np.load(root/'results/attractor-novelty/reference/trajectories/Lorenz.npz')['points'])
checks['portable_catalog_lorenz_identical']=True
sources=[]
for pattern in ['*gaussian_vortex.py','check_vortex_regular.py','build_*references.py','compare_*shapes.py','render_novelty_comparison.py','import_kicked_top_references.py']:
    for path in (root/'research').glob(pattern):
        ast.parse(path.read_text(encoding='utf-8'));sources.append(str(path.relative_to(root)).replace('\\','/'))
checks['parsed_sources']=sorted(set(sources))
checks['scope']='Separate portable-wrapper replay: exact arrays and identical rendered PNG/GIF, plus source parsing and regenerated catalog control. Not a new chaos or priority proof.'
(base/'replay-checks.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8')
gif=Image.open(root/'docs/assets/gaussian-vortex-triskelion.gif')
sheet=Image.new('RGB',(900,640),'#070b12')
for j,index in enumerate([0,12,24,36,48,60]):
    gif.seek(index);frame=gif.convert('RGB').resize((300,320));sheet.paste(frame,((j%3)*300,(j//3)*320))
sheet.save(scratch/'animation-review.png')
print(json.dumps({'exact_replays':len(checks)-2,'parsed_sources':len(checks['parsed_sources']),'animation_frames':gif.n_frames}))
