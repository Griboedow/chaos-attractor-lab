"""Complete accounting and final-window recurrence; no infinite-time classification."""
from pathlib import Path
import hashlib,json
import numpy as np

def read(path):return json.loads(Path(path).read_text())

rows=[]
for folder,shapes in [('rectangular-long01',[(2,3),(3,2),(2,4),(4,2)]),('dual-long01',[(2,2),(2,3),(3,2)]),('balanced-long01',[(3,2)])]:
    for m,n in shapes:
        filename=Path(folder)/f'{m}x{n}.json';data=read(filename);arr=np.load(Path(folder)/f'{m}x{n}-final-50000-states.npz')
        for r,z in zip(data['runs'],arr['states']):
            if not np.isfinite(z).all():
                rows.append(dict(stage=folder,shape=[m,n],case=r['case'],completed=False,accepted_steps=r['completed_steps']));continue
            tail=z[-6144:];scale=1+float(np.max(abs(tail)));tolerance=1e-9*scale;found=None;residual=None
            for lag in range(1,2049):
                e=float(np.max(abs(tail[-4096:]-tail[-4096-lag:-lag])))
                if e<tolerance:found=lag;residual=e;break
            rows.append(dict(stage=folder,shape=[m,n],case=r['case'],completed=True,last_growth=r['blocks'][-1]['rate'],recurrence_lag=found,recurrence_residual=residual,tolerance=tolerance,max_abs_state=float(np.max(abs(tail))),parent_sha256=hashlib.sha256(filename.read_bytes()).hexdigest()))
campaigns=[]
for folder in ('pilot01','dual01','balanced01'):
    datasets=[read(p) for p in sorted(Path(folder).glob('*.json'))];runs=[r for d in datasets for r in d['runs']]
    if folder=='pilot01':
        captures=sum(bool(r['capture_step']) for r in runs);completed=len(runs)-captures;escapes=0;flags=sum(not r['capture_step'] and r['finite_rate']>1e-5 for r in runs)
    else:
        captures=sum(r['status']=='numerical capture' for r in runs);completed=sum(r['status']=='completed' for r in runs);escapes=sum(r['status']=='escape threshold' for r in runs);flags=sum(bool(r['rates']) and r['rates']['mean']>.005 for r in runs)
    campaigns.append(dict(folder=folder,inputs=len(runs),numerical_captures=captures,completed=completed,escape_threshold=escapes,positive_flags=flags,noncontracting_volume=sum(abs(r.get('determinant',0))>=1 for r in runs)))
controls=read('controls-balanced01/summary.json')['runs'];groups=[]
for case in (150,987):
    rr=[r for r in controls if r['case']==case]
    for category in ('original','ordinary','endpoint','fresh','rate'):
        selected=[r for r in rr if r['label']==category or r['label'].startswith(category+'-')]
        groups.append(dict(case=case,category=category,count=len(selected),completed=sum(r['status']=='completed' for r in selected),escaped=sum(r['status']=='escape threshold' for r in selected),captured=sum(r['status']=='numerical capture' for r in selected),positive=sum(bool(r['rates']) and r['rates']['mean']>0 for r in selected)))
result=dict(campaigns=campaigns,final_windows=rows,controls=groups,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path('analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(campaigns=campaigns,recurrences=[r for r in rows if r.get('recurrence_lag')],controls=groups)))
