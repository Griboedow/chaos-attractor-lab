"""Check relative closures with state-only adaptive solvers and register final chaotic-band cores."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from adaptive import polar_ring
from periodic_catalog import shifts
from relative_cores import rotation


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True);ap.add_argument('--packet',type=Path,required=True)
    a=ap.parse_args();output=a.packet/'nested-skeletons.json'
    m=json.loads(output.read_text());by_order={r['period_order']:r for r in m['cores'] if r['period_order']==2}
    records=[]
    for file in sorted(a.source.glob('order-*.json')):
        d=json.loads(file.read_text());record=dict(d)
        if d['completed']:
            p=np.array(d['parameters']);v=np.array(d['initial']);t=d['sector_period']
            R=rotation(2*np.pi*d['relative_rotation_index']/7);F=lambda time,loc:polar_ring(loc,p)
            checks=[]
            for method,rtol in [('RK45',2e-12),('DOP853',4e-13)]:
                got=solve_ivp(F,(0,t),v,method=method,rtol=rtol,atol=rtol/100)
                if not got.success:raise RuntimeError(got.message)
                checks.append(dict(method=method,rtol=rtol,absolute_rotated_sector_closure=float(np.linalg.norm(got.y[:,-1]-R@v))))
            pts=np.load(d['path'])['points'];floor,k=shifts(pts)
            record.update(independent_sector_closures=checks,minimum_phase_shift_floor_rms=floor,worst_divisor=k)
            if max(z['absolute_rotated_sector_closure'] for z in checks)<2e-6 and floor>.001:
                item={k:v for k,v in record.items() if k!='history'};by_order[d['period_order']]=item
        records.append(record)
    output.write_text(json.dumps(dict(case=180844,cores=[by_order[k] for k in sorted(by_order)]),indent=2)+'\n')
    (a.packet/'relative-core-validation.json').write_text(json.dumps(records,indent=2)+'\n')
    print([(r['period_order'],r['dominant_multiplier']) for r in by_order.values()])


if __name__=='__main__':
    main()
