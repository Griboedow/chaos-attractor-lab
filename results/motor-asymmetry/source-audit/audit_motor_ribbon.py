"""Source-form, solver and row-distinctness checks for a stored motor trace."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from motor_rows import field_jac,source_force_constraint_field


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    a=np.load(args.source);p=a['parameters'];x=a['states'];worst=0.
    ids=np.rint(np.linspace(0,len(x)-1,17)).astype(int)
    for i in ids:
        f,j=field_jac(x[i],p,3);independent=source_force_constraint_field(x[i],p,3)
        jc=np.column_stack([source_force_constraint_field(x[i]+1e-25j*np.eye(8)[k],p,3).imag/1e-25 for k in range(8)])
        worst=max(worst,float(np.max(abs(f-independent))/(1+np.max(abs(f)))),float(np.max(abs(j-jc))/(1+np.max(abs(j)))))
    assert worst<1e-11
    rows=[]
    for i in ids[[0,8,16]]:
        ends={}
        for method in ['DOP853','Radau']:
            field=(lambda t,y:field_jac(y,p,3)[0]) if method=='DOP853' else (lambda t,y:source_force_constraint_field(y,p,3))
            sol=solve_ivp(field,(0,2),x[i],method=method,rtol=2e-11,atol=2e-13,max_step=.005)
            assert sol.success;ends[method]=sol.y[:,-1]
        error=float(max(abs(ends['DOP853']-ends['Radau'])));assert error<1e-8
        rows.append(dict(sample=int(i),time=float(a['time'][i]),error=error))
    positions=np.column_stack([x[:,6:8],-x[:,6:8].sum(axis=1)])
    distances=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        squared=(x[:,i]-x[:,j])**2+(x[:,3+i]-x[:,3+j])**2+(positions[:,i]-positions[:,j])**2
        distances.append(dict(pair=[i,j],minimum=float(np.sqrt(squared.min())),rms=float(np.sqrt(squared.mean()))))
    record=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        parameters=p.tolist(),source_field_and_complex_step_jacobian_relative_error=worst,sample_ids=ids.tolist(),
        short_independent_source_Radau_vs_DOP853=rows,row_distances=distances,
        row_distance_scope='Finite observations in (c_j,s_j,X_j); not a proof of separation for all time.',
        maximum_disk_radius_squared=float(np.max(x[:,:3]**2+x[:,3:6]**2)),disk_bound=float(1+p[2]**2),
        source_rate_amplitude=float(-.1/np.sqrt(1+p[2]**2)),source_rate_harmonic_amplitude=.1,
        warning='No exact cluster reduction or novel topology is inferred from these numerical checks.')
    (out/'audit.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf8')
    for name in ['audit_motor_ribbon.py','motor_rows.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(json.dumps(record))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
