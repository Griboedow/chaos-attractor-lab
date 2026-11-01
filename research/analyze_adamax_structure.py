"""Audit transient false positives and exact algebraic Adamax structure."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.stats import qmc
from adamax_map import step,step_tangent


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    scout=np.load(args.source/'adamax-scout/scout.npz')
    primary=(scout['status']==0)&(scout['windows'].min(1)>.005)&(np.diff(scout['extrema'][:,:,0],axis=1)[:,0]>1e-4)
    tail=scout['tails'];zero=np.all(tail[:,:,:2]==0,axis=(1,2));narrow=np.ptp(tail[:,:,0],axis=1)<1e-8
    audit=dict(primary_count=int(primary.sum()),exact_zero_x_m_tail_ids=np.flatnonzero(primary&zero).tolist(),
        narrow_x_tail_threshold=1e-8,narrow_x_tail_ids=np.flatnonzero(primary&narrow).tolist(),
        interpretation='A positive tangent exponent is insufficient: an exactly zero floating-point x,m trajectory stays on the invariant zero line. Narrow tails are flagged, not automatically proved regular.',
        examples=[dict(id=int(i),parameters=scout['parameters'][i].tolist(),windows=scout['windows'][i].tolist(),
            end=scout['ends'][i].tolist(),tail_extent=np.ptp(tail[i],axis=0).tolist(),observed_reset_count=int(scout['resets'][i])) for i in np.flatnonzero(primary&zero)])
    (out/'scout-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    rng=np.random.default_rng(2026110124);errors=[];cycle_rows=[]
    for _ in range(256):
        b,c=rng.uniform(.02,.98,2);k=(1-b)/(1+b);e=rng.uniform(.00001,.8)*k/2;a=k/2-e
        p=np.array([b,c,e]);positive=np.array([a,-k*a,k/2]);negative=np.array([-a,k*a,k/2])
        zp,Jp,rp=step_tangent(positive,np.eye(3),p);zn,Jn,rn=step_tangent(negative,np.eye(3),p)
        assert rp and rn
        B=np.array([[1-2*b-4*e/k,-2*b/k],[1-b,b]])
        error=max(np.max(abs(zp-negative)),np.max(abs(zn-positive)),np.max(abs((Jn@Jp)[:2,:2]-B@B)),
            abs(np.linalg.det(B)-b*(3-4*e/k)),abs(np.trace(B)-(1-b-4*e/k)))
        assert error<1e-11;errors.append(float(error))
    e=5e-6;critical=float(np.roots([3+4*e,-(4-4*e),1]).min())
    for b in [.1,critical,.3394757778476924,.36351296468637884,.9]:
        k=(1-b)/(1+b);B=np.array([[1-2*b-4*e/k,-2*b/k],[1-b,b]])
        cycle_rows.append(dict(beta1=b,amplitude=k/2-e,nonzero_two_step_multipliers=[[v.real,v.imag] for v in np.linalg.eigvals(B@B)],
            spectral_radius=float(max(abs(np.linalg.eigvals(B@B))))))
    structure=[]
    for ident in [847,239,4099,4096]:
        p=scout['parameters'][ident];b,c,e=p
        states=np.load(args.source/f'adamax-previews/case-{ident}.npy',mmap_mode='r')
        prev=states[:-1];current=states[1:];reset=abs(prev[:,0])+e>=c*prev[:,2]
        X,M,U=current[reset].T;sign=np.sign(prev[reset,0])
        residual=M-sign*U*(U-e)+X*U
        max_scaled=float(np.max(abs(residual)/(1+abs(M)+abs(X*U)+U*U)))
        assert max_scaled<1e-12
        age=-1;ages=[]
        for is_reset in reset:
            if is_reset:age=0
            elif age>=0:age+=1
            if age>=0:ages.append(age)
        counts=np.bincount(ages)
        Umax=float(states[:,2].max());L=int(np.floor(np.log(e/Umax)/np.log(c)))+1
        assert c**L*Umax<e
        row=dict(id=ident,reset_count=int(reset.sum()),transitions=len(reset),reset_sheet_max_scaled_residual=max_scaled,
            reset_age_counts=counts.tolist(),observed_maximum_reset_age=len(counts)-1,
            finite_envelope_u_max=Umax,conditional_memory_run_bound=L,
            conditional_bound_scope='If this u upper bound holds forever, no run of L consecutive memory steps is possible. Observed extrema alone do not prove a global bound.')
        if b<c:
            B=max(.73,c/(c-b));row['proved_cold_bias_corrected_bound']=dict(x_abs=B,m_abs=B,u_max=B+e)
        structure.append(row)
    # Close an observed four-cycle from a separately declared cold-start cohort.
    cohort=np.load(args.source/'adamax-cohorts/histories.npz');meta=json.loads((args.source/'adamax-cohorts/analysis.json').read_text())
    index=next(i for i,r in enumerate(meta['runs']) if r['id']==239 and r['kind']=='cold_start' and r.get('tail_period_up_to64')==4)
    p=cohort['parameters'][index]
    def return_map(z):
        for _ in range(4):z=step(z,p)
        return z
    solution=root(lambda z:return_map(z)-z,cohort['ends'][index],tol=1e-11)
    z=solution.x.copy();orbit=[];product=np.eye(3);margins=[]
    for _ in range(4):
        orbit.append(z.copy());margins.append(abs(abs(z[0])+p[2]-p[1]*z[2]));z,J,_=step_tangent(z,np.eye(3),p);product=J@product
    closure=float(np.max(abs(z-solution.x)));eig=np.linalg.eigvals(product);assert closure<1e-12 and max(abs(eig))<1
    cycle=dict(id=239,selected_cohort_index=index,initial_amplitude=meta['runs'][index]['amplitude'],parameters=p.tolist(),orbit=np.array(orbit).tolist(),
        closure=closure,solver_success=bool(solution.success),minimum_switch_margin=min(margins),
        nontrivial_lag_distances=[float(np.linalg.norm(np.array(orbit)[k]-orbit[0])) for k in [1,2,3]],
        multipliers=[[float(v.real),float(v.imag)] for v in eig],spectral_radius=float(max(abs(eig))),
        scope='Floating-point closed attracting four-cycle, coexisting numerically with the selected chaotic trajectory; not an interval existence proof.')
    result=dict(period_two_controls=dict(samples=len(errors),maximum_absolute_identity_error=max(errors),critical_beta1=critical,
        critical_scope='Linear unit-modulus crossing only. No nonlinear bifurcation coefficient or Neimark-Sacker theorem evaluated.',examples=cycle_rows),
        sheet_geometry=structure,coexisting_four_cycle=cycle)
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in ['analyze_adamax_structure.py','adamax_map.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    (out/'sources.json').write_text(json.dumps({str(p.relative_to(args.source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.source/'adamax-scout/scout.npz',args.source/'adamax-cohorts/histories.npz',args.source/'adamax-cohorts/analysis.json']},indent=2)+'\n')
    print('Primary',audit['primary_count'],'zero-line',audit['exact_zero_x_m_tail_ids'],'narrow',len(audit['narrow_x_tail_ids']),flush=True)
    print('Period-two identities',max(errors),'four-cycle closure',closure,'radius',cycle['spectral_radius'],flush=True)
    print('Reset ages',[(r['id'],r['observed_maximum_reset_age']) for r in structure],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
