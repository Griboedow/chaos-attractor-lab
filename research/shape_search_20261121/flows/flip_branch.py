"""Locate the first flip of the seven-crossing cycle and search its attracting descendants."""
import argparse
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from band_continuation import continuation,dynamics
from cycle_scan import resample
from exact_knot import exact_diagram,prepare
from field import integrate,ring_tangent
from knot_diagram import alexander
from orbits import refine


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--packet',type=Path,default=Path('results/shape_search_20261121/flows'))
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((a.packet/'manifest.json').read_text())
    base=next(r for r in manifest['stable_cycles'] if r['case']==28756)
    target=next(r for r in manifest['cases'] if r['case']==180844)
    start=np.array(base['initial']);t=base['period'];p0=np.array(base['parameters'])
    p1=np.array(target['parameters']);scale=.89
    def parameters(x):
        p=p0.copy()
        p[8]=(1-x)*p0[8]+x*p1[8]
        return p
    def track(frac,guess,period):
        p=parameters(frac)
        axis=int(np.argmax(abs(ring_tangent(np.r_[guess,np.zeros(3)],p)[:3])))
        sol=refine(axis,guess,period,scale,p,1)
        if sol is None:raise RuntimeError('Shooting failed at '+str(frac))
        seed,per=sol;stats,_=dynamics(seed,per,p,True)
        return seed,per,stats
    _,left=continuation(p0,parameters(.63),start,t,scale)
    near,t_left=left
    frac_lo,frac_hi=.63,.725
    left,t_left,m_left=track(frac_lo,near,t_left)
    right,t_right,m_right=track(frac_hi,left,t_left)
    assert m_left['dominant'][0]>-1 and m_right['dominant'][0]<-1
    bracket=[]
    for i in range(23):
        center=(frac_lo+frac_hi)/2
        ratio=(center-frac_lo)/(frac_hi-frac_lo)
        trial,period,m=track(center,(1-ratio)*left+ratio*right,(1-ratio)*t_left+ratio*t_right)
        bracket.append(dict(fraction=center,e=float(parameters(center)[8]),**m))
        if m['dominant'][0]>-1:
            frac_lo,left,t_left,m_left=center,trial,period,m
        else:
            frac_hi,right,t_right,m_right=center,trial,period,m
    report=dict(base=28756,final=180844,parameter_index=8,parameter_name='e',
                bracket=[dict(fraction=x,e=float(parameters(x)[8]),**m) for x,m in [(frac_lo,m_left),(frac_hi,m_right)]],
                bisection_steps=bracket,trials=[])
    save=lambda:(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    previous,previous_t=left,t_left
    for frac in (.700,.720,.728,.736,.744,.752,.760,.768,.776,.784,.792,.80,.81,.82,.84,.86,.88,.89,.895,.900,.905,.910,.915,.920,.925,.93,.935,.94,.945,.95,.955,.96,.97,.98,.99,1):
        p=parameters(frac)
        # An offset in physical space excites the unstable transverse mode of the original branch.
        loc,per,main=track(frac,previous,previous_t)
        previous,previous_t=loc,per
        z0=loc+scale*np.array([.003,-.0015,.0027])
        step=.01
        rates,z=integrate(p,z0,step,round(4000/step),round(230/step),6,2,ring_tangent)
        row=dict(fraction=frac,e=float(p[8]),primary_period=per,primary_refined_start=loc.tolist(),primary_multiplier=main['dominant'][0],
                 rate_windows=rates.tolist(),initial=z0.tolist())
        axis=int(np.argmax(np.ptp(z,axis=0)))
        level=np.median(z[:,axis])
        hits=np.flatnonzero((z[:-1,axis]<level)&(z[1:,axis]>=level))
        alpha=(level-z[hits,axis])/(z[hits+1,axis]-z[hits,axis])
        crossings=z[hits]+alpha[:,None]*(z[hits+1]-z[hits]);times=(hits+alpha)*step*2
        candidates=[]
        for lag in range(1,min(120,len(times)-8)):
            periods=times[-8:]-times[-8-lag:-lag]
            if periods.mean()>per*18:break
            distance=np.linalg.norm(crossings[-8:]-crossings[-8-lag:-lag],axis=1)/scale
            if max(distance)<.0015 and np.ptp(periods)<.002*periods.mean():
                candidates.append((lag,float(periods.mean()),float(max(distance))))
        for lag,period,error in candidates[:3]:
            sol=refine(axis,crossings[-1-lag],period,scale,p,1)
            if sol is None:continue
            seed,period=sol;stats,dense=dynamics(seed,period,p,True)
            if stats['modulus']>=.998 or stats['closure']/scale>1e-6:continue
            # A near-return after half of a candidate period can signal convergence to the original branch.
            phase=np.linspace(0,period,9,endpoint=False)
            gaps=[np.max(np.linalg.norm(dense.sol(phase)[:3]-dense.sol((phase+period/k)%period)[:3],axis=0))/scale for k in range(2,33)]
            if min(gaps)<1e-3:continue
            row['attracting_cycle']=dict(**stats,refined_start=seed.tolist(),phase_shift_floor=float(min(gaps)),
                                          period_ratio=period/per,original_section_error=error)
            if 1.7<period/per<2.3 and not any('witness_file' in rr for rr in report['trials']):
                orbit=dense.sol(np.linspace(0,period,32768,endpoint=False))[:3].T
                arc=resample(orbit,640);views=[]
                for camera in (17,35,81,192,261,713,967,1217,2114,3117,5032,7101):
                    d=alexander(arc,camera)
                    if d and d['crossings']<85:views.append((d['crossings'],camera))
                cert_views=[]
                for _,camera in sorted(views)[:2]:
                    vertices,mat,shift,rms=prepare(arc,camera)
                    try: exact=exact_diagram(vertices,mat)
                    except ValueError:continue
                    exact.pop('minimum_depth',None);cert_views.append(dict(integer_matrix=mat.tolist(),**exact))
                if len(cert_views)==2 and cert_views[0]['coefficients']==cert_views[1]['coefficients']:
                    cert=dict(quantization_per_rms=10**9,center=shift.tolist(),rms=float(rms),
                              integer_vertices=vertices.tolist(),views=cert_views)
                    file='cycle-'+str(round(frac*1000))
                    np.savez_compressed(a.out/(file+'.npz'),points=orbit)
                    (a.out/(file+'.json')).write_text(json.dumps(cert,separators=(',',':'))+'\n')
                    row.update(witness_file=file+'.json',cloud_file=file+'.npz',coefficients=cert_views[0]['coefficients'])
            break
        report['trials'].append(row);save()
        print(frac,main['dominant'][0],row.get('attracting_cycle',{}).get('period_ratio'),flush=True)


if __name__=='__main__':
    main()
