"""Adversarial section perturbations along the transverse direction of the return tube."""
import json
import numpy as np
from scipy.spatial import cKDTree
from native_section import RAW,RESULT,trace,packed,orbit_section
from return_tube import graph,CENTER,ZLO,ZHI,WIDTH
from compact_tube import ZL,ZH
from summary import metrics

def main():
    g=graph();p=packed();nodes=orbit_section()
    radii=[0.,3e-5,.0003,.001,.003,.01,.02,.04,.08]
    zgrid=np.linspace(ZL+.00003,ZH-.00003,16);rows=[]
    for h in [.02,.01]:
        for radius in radii:
            sel=[]
            for j,z in enumerate(zgrid):
                for sign in ([1] if radius==0 else [-1,1]):
                    d=sign*radius;x=g(z-CENTER)+d
                    events=trace(np.array([x,0.,z]),h,round(14000/h),p)
                    m=metrics(events,nodes)
                    last=events[-150:,1:]
                    mask=(last[:,0]>.78)&(last[:,0]<.83)&(last[:,1]>.96)&(last[:,1]<1.03)
                    v=last[mask]
                    hits=len(v)
                    inside=(len(v)>0 and np.all((v[:,1]>ZLO)&(v[:,1]<ZHI)&(abs(v[:,0]-g(v[:,1]-CENTER))<WIDTH)))
                    fate='cycle' if m['capture_time'] is not None else ('cloud' if hits==30 and inside else 'other')
                    sel.append(dict(step=h,radius=radius,sign=sign,z=float(z),x=float(x),fate=fate,
                         copy=m['captured_copy'],time=m['capture_time'],five_band_hits=hits))
            rows+=sel
            print(h,radius,*[sum(r['fate']==m for r in sel) for m in ('cycle','cloud','other')],flush=True)
            (RAW/'transverse-stress.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':main()
