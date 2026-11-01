"""Replay numerical geometry and parameter slices; check independent formulas."""
import argparse,hashlib,itertools,json,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from search_product_alphabet import step
from search_soft_sphere import independent,step as full_step


def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);args.replays.mkdir(parents=True,exist_ok=False)
    for name in [Path(__file__).name,'search_product_alphabet.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py','follow_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:states=f['states'];beta,tau=f['params'];levels=f['levels'][:,:2]
    points=np.array(list(itertools.product(*levels)));rng=np.random.default_rng(2026110157)
    with np.load(args.critical/'critical-images.npz') as f:critical=f['points']
    probes=np.vstack([states[rng.choice(len(states),64,replace=False)],critical[rng.choice(len(critical),64,replace=False)]])
    errors=[];fd_errors=[];jac_errors=[]
    for x in probes:
        y,J,*_=step(x,beta,tau,levels,2);z,K,*_=full_step(x,beta,points,tau)
        errors.append(float(max(abs(y-z))));jac_errors.append(float(np.max(abs(J-K))/max(1,np.max(abs(J)))))
        h=1e-7;fd=np.column_stack([(independent(x+h*np.eye(3)[j],beta,points,tau)-independent(x-h*np.eye(3)[j],beta,points,tau))/(2*h) for j in range(3)])
        fd_errors.append(float(np.max(abs(J-fd))/max(1,np.max(abs(J)))))
    assert max(errors)<1e-11 and max(jac_errors)<1e-10 and max(fd_errors)<1e-5
    specs=[(args.scan,'scan_product_mechanism.py',['--source',str(args.source)],['scan.npz'],['cases','strictly_positive_cohorts','finite_periodic_cohorts','flagged_cohorts']),
           (args.critical,'trace_product_critical_surfaces.py',['--source',str(args.source),'--cycle',str(args.cycle),'--diagnostics',str(args.diagnostics/'diagnostics.npz')],['critical-images.npz','cluster-labels.npz'],['images','numerical_clusters','contraction_bound'])]
    replays=[]
    for folder,script,command,arrays,keys in specs:
        dest=args.replays/folder.name
        subprocess.run([sys.executable,str(folder/script),*command,'--out',str(dest)],check=True,capture_output=True,text=True)
        old=read(folder/'report.json');new=read(dest/'report.json')
        for key in keys:assert old[key]==new[key],key
        for name in arrays:
            with np.load(folder/name) as a,np.load(dest/name) as b:
                assert a.files==b.files
                for key in a.files:assert a[key].dtype==b[key].dtype and a[key].shape==b[key].shape and a[key].tobytes()==b[key].tobytes(),key
        replays.append(dict(run=folder.name,report_sha256=sha(folder/'report.json'),identical_arrays=arrays,identical_report_keys=keys))
    with np.load(args.diagnostics/'diagnostics.npz') as f:values=f['values'];signs=f['signs']
    # A separate NumPy tangent pass through the full-site Jacobian for 2048 initial states.
    vector=np.array([1.,.371,-.217]);vector/=np.linalg.norm(vector);log_errors=[]
    for k in range(2048):
        _,J,*_=full_step(states[k],beta,points,tau);vector=J@vector;length=np.linalg.norm(vector);vector/=length
        log_errors.append(abs(float(np.log(length))-values[k,2]))
    assert max(log_errors)<1e-9
    data=read(args.diagnostics/'report.json');active=(values[1000:,3:6]>1).sum(1)
    for row in data['conditional_stretch']:
        mask=active==row['active_coordinates'];assert mask.sum()==row['count']
        assert float(values[1000:,2][mask].mean())==row['mean_leading_log']
    # Exact rational upper bound on every site norm; hence on every posterior mean.
    M=Q(201,100);R=Q(40);b=Q(float(beta));squared=max(sum(Q(float(v))**2 for v in p) for p in points)
    assert squared<M*M and 0<b<1
    forcing=b+abs(1-2*b)*M;gap=(1-b)*R-forcing;assert gap>0
    widths=[]
    for lo,hi in levels:
        a=(hi-lo)/2; widths.append(float(tau/(2*a)*np.arccosh(np.sqrt(2*a*a/tau))))
    report=dict(status='passed',replays=replays,independent_probes=len(probes),seed=2026110157,
                maximum_formula_error=max(errors),maximum_jacobian_relative_error=max(jac_errors),maximum_finite_difference_relative_error=max(fd_errors),
                full_site_tangent_steps=2048,maximum_log_growth_error=float(max(log_errors)),
                absorbing_bound=dict(site_norm_upper=str(M),site_squared_norm_max=str(squared),beta=str(b),radius=str(R),forcing_upper=str(forcing),strict_boundary_gap=str(gap),
                                     scope='All indefinitely defined forward orbits eventually enter radius 40. The ball includes normalization singularities and is not a compact smooth-domain trapping certificate.'),
                numerical_switching_halfwidths=widths,scope=__doc__)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','cycle','diagnostics','critical','scan','out','replays']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
