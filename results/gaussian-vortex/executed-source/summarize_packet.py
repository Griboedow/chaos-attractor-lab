from pathlib import Path
import json
import numpy as np
root=Path('D:/Codex/chaos-attractor-lab');base=root/'results/gaussian-vortex'
data=np.load(base/'long-runs.npz');render=np.load(base/'render-run.npz')
rows=[]
for i in range(6):
    rows.append({'candidate_id':int(data['ids'][i]),'h':float(data['stepsizes'][i]),
                 'spectrum':data['windows'][i].mean(0).tolist(),
                 'minimum_leading_window':float(data['windows'][i,:,0].min()),
                 'maximum_abs_divergence_error':float(np.max(np.abs(data['windows'][i].sum(1)+3*data['parameters'][i,0]))),
                 'max_observed_radius':float(data['max_radius'][i])})
cohorts={}
for name in ['fresh-starts','neighbors']:
    d=np.load(base/f'{name}.npz');passed=d['windows'][:,:,0].min(1)>.005
    cohorts[name]={'count':len(passed),'positive_all_windows':int(passed.sum()),
                  'below_gate_indices':np.where(~passed)[0].tolist(),
                  'gate':'Every 500-time-unit leading-exponent window > 0.005',
                  'maximum_abs_divergence_error':float(np.max(np.abs(d['windows'].sum(2)+3*d['parameters'][:,0,None])))}
summary={'long_runs':rows,'render_spectrum':render['windows'].mean(0).tolist(),
         'render_minimum_leading_window':float(render['windows'][:,0].min()),
         'render_max_observed_radius':float(render['max_radius']),
         'cohorts':cohorts}
(base/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
