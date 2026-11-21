"""Collect numerical returns, escape attempts and independent exact section orientation bounds."""
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
from native_section import CACHE,RAW,RESULT,P,CASES
from return_tube import graph
from ensemble import starts

def load(p):return json.loads(p.read_text(encoding='utf-8'))
def save(name,value):(RESULT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def numbers(s):return [float(v) for v in s]
def box(x):return (x,x)
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mult(x,y):
    a=[v*w for v in x for w in y];return min(a),max(a)
def orientation():
    # Exact-rational interval arithmetic; graph and field coefficients are the registered binary64 values.
    F0=lambda f:F.from_float(float(f));coeff=[F0(t) for t in graph().coef]
    z=(F('.98045'),F('.99945'));v=add(z,box(-F0(.99)));x=box(coeff[-1])
    for a in reversed(coeff[:-1]):x=add(mult(x,v),box(a))
    x=add(x,(-F('3e-5'),F('3e-5')))
    w=add(box(F0(P[3])),add(mult(box(F0(P[4])),z),mult(box(F0(P[5])),mult(x,x))))
    dy=mult(x,w)
    assert dy[0]>0
    points=np.load(RAW/'stable-orbit-section.npy')[:,1:]
    bounds=np.array([[float(x[0]),float(z[0])],[float(x[1]),float(z[1])]])
    sep=float(np.linalg.norm(np.maximum(0,np.maximum(bounds[0]-points,points-bounds[1])),axis=1).min())
    return dict(interpretation='Exact-rational interval evaluation of the continuous vector field at native y=0 across the entire proposed chart. Display intervals have been widened outward. Separation uses numerical cycle-section coordinates.',
        chart_x_outward=[np.nextafter(float(x[0]),-np.inf),np.nextafter(float(x[1]),np.inf)],
        chart_z_outward=[np.nextafter(float(z[0]),-np.inf),np.nextafter(float(z[1]),np.inf)],
        positive_normal_native_outward=[np.nextafter(float(dy[0]),-np.inf),np.nextafter(float(dy[1]),np.inf)],
        exact_lower_normal_as_rational=[str(dy[0].numerator),str(dy[0].denominator)],
        minimum_box_to_numerical_cycle_section_distance=sep)
def main():
    rs=load(RAW/'ensemble-summary.json');runs=rs['runs']
    groups=[]
    for key in sorted(set((r['kind'],r['step'],r['duration']) for r in runs)):
        group=[r for r in runs if (r['kind'],r['step'],r['duration'])==key];cap=[r for r in group if r['capture_time'] is not None]
        live=[r for r in group if r['capture_time'] is None]
        groups.append(dict(kind=key[0],h=key[1],duration=key[2],runs=len(group),captures=len(cap),
            final_capture_times=[r['capture_time'] for r in cap],
            minimum_noncapture_section_distance=min(r['cycle_distances_quantiles'][0] for r in live),
            events=sum(r['events'] for r in group)))
    stress=load(RAW/'transverse-stress.json')
    sg=[]
    for h,radius in sorted(set((v['step'],v['radius']) for v in stress)):
        group=[r for r in stress if r['step']==h and r['radius']==radius];counts=Counter(v['fate'] for v in group)
        sg.append(dict(h=h,radius=radius,total=len(group),cloud=counts['cloud'],cycle=counts['cycle'],other=counts['other']))
    save('fates.json',dict(capture_definition='Every final crossing in a consecutive suffix of at least 128 is within 2e-5 native Euclidean distance of the four positive-orientation section points of one rotated stable orbit. The time is the first event in this final suffix.',
        reference_positive_section_points=rs['section_points'],chaotic_phase_seed_method='64 different pseudo-random samples without replacement from one saved chaotic trajectory. The first 32 are used at two matched steps. The last 32 are used at a longer horizon.',seed=564720261121,
        starts=starts(),ensemble=groups,transverse_stress=sg,
        stress_definition='16 equally spaced z coordinates at each signed native x offset; zero counted once. 14,000 time units. Cloud outcome requires 30 of 150 final crossings in the chart, all inside its defined tube.'))
    p=load(RAW/'tube-n4097-l17-h0.00625.json');p['coarse_check']=load(RAW/'tube-n601-l9-h0.0125.json')
    p['adaptive']=load(RAW/'adaptive.json');p['short_equivariance_map']=load(RAW/'reduced-n2049-l17-h0.00625.json')
    p['exact_section_orientation']=orientation()
    save('return-map.json',p)
    proof={}
    for f in sorted(RAW.glob('*.json')):
        if f.name.startswith(('arb-','affine-')):
            proof[f.stem]=load(f)
    save('interval-returns.json',dict(python_flint_version='0.9.0',certificate_scope='Local declared input boxes only; root and flow enclosures include all real initial states in those boxes. No domain-wide return assertion.',records=proof))
    files=[]
    for f in sorted(RAW.iterdir()):
        if f.is_file() and f.suffix in ('.npz','.npy','.json'):
            files.append(dict(cache_relative=f.relative_to(CACHE).as_posix(),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),bytes=f.stat().st_size))
    save('manifest.json',dict(cache_root=str(CACHE),case=180844,parameters=P.tolist(),coordinate='Cartesian (x,y,z), native units',
          data_files=files,source_clouds={k:dict(cache_relative=CASES[k]['cloud_relative'],sha256=CASES[k]['sha256']) for k in ('stable','coexisting_chaotic_cloud')}))
    print('Saved',len(files),'raw records. Groups',groups,'stress',sg)

if __name__=='__main__':main()
