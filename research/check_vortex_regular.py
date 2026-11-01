"""Continue every regular-looking fresh-start outcome with adaptive event timing."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from check_gaussian_vortex import independent_field


def main():
    base=Path(__file__).resolve().parents[1]/'results/gaussian-vortex'
    data=np.load(base/'fresh-starts.npz')
    ids=np.where(data['windows'][:,:,0].mean(1)<.001)[0]
    rows=[]
    for i in ids:
        p=data['parameters'][i]
        component=int(np.argmax(np.var(data['tail'][i],axis=0)))
        level=float(data['tail'][i,:,component].mean())
        def crossing(t,z):
            return z[component]-level
        crossing.direction=1
        sol=solve_ivp(lambda t,z:independent_field(z,p)[0],(0,1000.),data['end'][i],
                      method='DOP853',rtol=2e-11,atol=2e-13,events=crossing,max_step=.2)
        assert sol.success
        points=sol.y_events[0];times=sol.t_events[0]
        assert len(points)>24
        options=[]
        for lag in range(1,9):
            errors=np.linalg.norm(points[-20:]-points[-20-lag:-lag],axis=1)
            periods=times[-20:]-times[-20-lag:-lag]
            options.append({'section_return_lag':lag,'max_return_error':float(errors.max()),
                            'mean_period':float(periods.mean()),'period_range':float(np.ptp(periods))})
        first=next((r for r in options if r['max_return_error']<1e-7),None)
        rows.append({'fresh_start_index':int(i),'component':component,'section_level':level,
                     'continuation_time':1000.,'crossings':len(points),'shortest_tested_close_return':first,
                     'all_lags':options,'last20_section_points':points[-20:].tolist(),
                     'last20_section_times':times[-20:].tolist()})
    (base/'regular-controls.json').write_text(json.dumps({'method':'DOP853','rtol':2e-11,'atol':2e-13,
        'max_step':.2,'scope':'Finite numerical periodic-return evidence, without interval existence or stability certification.',
        'outcomes':rows},indent=2)+'\n',encoding='utf-8')
    print(json.dumps([{'index':r['fresh_start_index'],'return':r['shortest_tested_close_return']} for r in rows]),flush=True)


if __name__=='__main__':
    main()
