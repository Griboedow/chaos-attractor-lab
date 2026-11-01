"""Finite-trajectory separation of the leading plane and strongest stable line.

Forward inverse-transpose iteration finds the normal to the leading plane;
backward inverse iteration finds the strongest contracting tangent direction.
This diagnoses a necessary feature of dominated splitting, never certifies it.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit
from twist_map import step,step_tangent


@njit
def compute(p,start,warm,observed,pad):
    z=start.copy()
    for _ in range(warm):z=step(z,p)
    total=observed+2*pad;states=np.empty((total+1,3));jac=np.empty((total,3,3));normal=np.empty((total+1,3));stable=np.empty_like(normal)
    states[0]=z;normal[0]=np.array([1.,2.,3.])/np.sqrt(14.)
    for i in range(total):
        z,j=step_tangent(z,np.eye(3),p);states[i+1]=z;jac[i]=j
        v=np.linalg.solve(j.T,normal[i]);normal[i+1]=v/np.linalg.norm(v)
    stable[-1]=np.array([3.,1.,2.])/np.sqrt(14.)
    for i in range(total-1,-1,-1):
        v=np.linalg.solve(jac[i],stable[i+1]);stable[i]=v/np.linalg.norm(v)
    angles=np.empty(observed);area=np.empty(observed);line=np.empty(observed);flowdefect=np.empty(observed)
    for k in range(observed):
        i=k+pad;angles[k]=np.arcsin(min(1.,abs(np.dot(normal[i],stable[i]))))
        area[k]=np.log(p[0]**3*np.linalg.norm(np.linalg.solve(jac[i].T,normal[i])))
        v=jac[i]@stable[i];line[k]=np.log(np.linalg.norm(v));v/=np.linalg.norm(v)
        flowdefect[k]=min(np.linalg.norm(v-stable[i+1]),np.linalg.norm(v+stable[i+1]))
    return states[pad:pad+observed],normal[pad:pad+observed],stable[pad:pad+observed],angles,area,line,flowdefect


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    d=json.loads(args.source.read_text());selected=d['candidates'][:2]
    protocol=dict(source=str(args.source),sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        candidates=[r['id'] for r in selected],warm=100000,observed=200000,padding_each_end=10000,
        method='Forward inverse-transpose plane normal and backward inverse tangent stable line; binary64 analytic map derivative; no trajectory interpolation.',
        quantities='Angle between strongest stable line and leading2D plane; one-step log area in plane and log length in line; sums over every sliding block.',
        scope='Finite-trajectory necessary diagnostics only. No uniform hyperbolicity, pseudohyperbolicity, dominated splitting or basin proof.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['check_twist_splitting.py','twist_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    rows=[]
    for selected_row in selected:
        ident=selected_row['id'];p=np.array(selected_row['parameters']);x,n,s,angle,area,line,defect=compute(p,np.array([.1,.2,.3]),100000,200000,10000)
        path=out/f'case-{ident}.npz';np.savez_compressed(path,states=x,normal=n,stable=s,angle_radians=angle,log_plane_area=area,log_stable_length=line,parameters=p)
        blocks=[]
        ca=np.r_[0.,np.cumsum(area)];cl=np.r_[0.,np.cumsum(line)]
        for length in [1,5,20,100,1000]:
            a=ca[length:]-ca[:-length];l=cl[length:]-cl[:-length]
            blocks.append(dict(length=length,minimum_mean_log_area=float(a.min()/length),maximum_mean_log_line=float(l.max()/length),
                fraction_nonpositive_area=float(np.mean(a<=0)),fraction_noncontracting_line=float(np.mean(l>=0))))
        row=dict(id=ident,angle_degrees_quantiles=np.rad2deg(np.quantile(angle,[0,.0001,.001,.01,.1,.5,1])).tolist(),
            quantiles=[0,.0001,.001,.01,.1,.5,1],mean_log_area=float(area.mean()),mean_log_line=float(line.mean()),
            maximum_line_invariance_defect=float(defect.max()),blocks=blocks,file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        assert row['maximum_line_invariance_defect']<1e-8
        rows.append(row);(out/'analysis.json').write_text(json.dumps(dict(runs=rows),indent=2)+'\n')
        print(row,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
