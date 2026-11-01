"""Check elementary Yogi identities and an exact initial switching discrepancy."""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path
import numpy as np
from numba import njit
from yogi_map import step,step_tangent


@njit
def bounded_run(p,start,bias):
    z=start.copy();maxima=np.abs(z);maximum_step=0.;pb=1.;pc=1.
    for _ in range(10000):
        pb*=p[0];pc*=p[1];new=step(z,p,1-pb if bias else 1.,1-pc if bias else 1.)
        maximum_step=max(maximum_step,abs(new[0]-z[0]));z=new;maxima=np.maximum(maxima,np.abs(z))
    return maxima,maximum_step


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);rng=np.random.default_rng(2026110144)
    inverse_errors=[];symmetry_errors=[];line_errors=[];bound_checks=[]
    for _ in range(256):
        b,c=rng.uniform(.02,.98,2);e=10**rng.uniform(-6,-1);p=np.array([b,c,e]);x=rng.choice([-1.,1.])*rng.uniform(.1,2);M=rng.uniform(-1,1)
        V=rng.uniform(max(c,1-c)+.01,2-c-.01)*x*x;target=np.array([x-M/(np.sqrt(V)+e),M,V])
        old_m=(M-(1-b)*x)/b
        for sign in [-1,1]:
            old_v=V+(1-c)*sign*x*x;old=np.array([x,old_m,old_v])
            assert old_v>=0 and np.sign(old_v-x*x)==sign
            inverse_errors.append(float(np.max(abs(step(old,p)-target))))
        z=np.r_[rng.uniform(-2,2,2),rng.uniform(.01,2)];R=np.array([-1.,-1.,1.])
        symmetry_errors.append(float(np.max(abs(step(z*R,p)-step(z,p)*R))))
        v=rng.uniform(.001,2);fixed=np.array([0.,0.,v]);value,J=step_tangent(fixed,np.eye(3),p);den=np.sqrt(v)+e
        analytic=np.array([[1-(1-b)/den,-b/den,0.],[1-b,b,0.],[0.,0.,1.]])
        line_errors.append(float(np.max(abs(J-analytic))));assert np.array_equal(value,fixed)
    assert max(inverse_errors)<1e-12 and max(symmetry_errors)==0 and max(line_errors)<1e-12
    for i in range(256):
        c=rng.uniform(.01,.99);b=rng.uniform(.001,.999)*np.sqrt(c);e=10**rng.uniform(-8,-1);p=np.array([b,c,e])
        x0=rng.choice([-1.,1.])*10**rng.uniform(-4,2);v0=10**rng.uniform(-8,2);start=np.array([x0,0.,v0])
        delta=min(c,1-c);K=(1-b)/(np.sqrt(delta)*(1-b/np.sqrt(c)));B=max(abs(x0),K/(1-b));Vmax=max(v0,(2-c)*B*B)
        for bias in [False,True]:
            maxima,jump=bounded_run(p,start,bias);ratios=maxima/np.array([B,B,Vmax]);assert ratios.max()<=1+1e-12 and jump<=K*(1+1e-12)
            bound_checks.append(dict(parameters=p.tolist(),start=start.tolist(),bias_correction=bias,bound_B=B,bound_V=Vmax,bound_step=K,
                observed_maxima=maxima.tolist(),observed_maximum_step=float(jump)))
    source=np.load(args.source);z=source['starts'][7591];x,m,v=map(Fraction.from_float,map(float,z));gap=v-x*x
    boundary=dict(id=7591,stored_start=z.tolist(),binary64_margin=float(z[2]-z[0]*z[0]),
        exact_margin_numerator=gap.numerator,exact_margin_denominator=gap.denominator,exact_margin_decimal_approx=float(gap),
        binary64_first_sign=0,real_promoted_first_sign=-1,
        conclusion='The exact real point represented by the stored three binary64 coordinates lies strictly below the switching surface; binary64 multiplication rounds its margin to zero. The two algorithms therefore choose different first branches.')
    assert gap<0 and boundary['binary64_margin']==0
    result=dict(inverse_controls=len(inverse_errors),inverse_max_residual=max(inverse_errors),symmetry_controls=len(symmetry_errors),
        symmetry_max_residual=max(symmetry_errors),fixed_line_controls=len(line_errors),fixed_line_jacobian_max_residual=max(line_errors),
        cold_bound_controls=len(bound_checks),cold_bound_checks=bound_checks,boundary_initialization=boundary,
        scope='Numerical implementation controls accompany algebraic arguments. The rational initial margin is exact. No chaotic-attractor or novelty theorem is claimed.',
        prior_art='v_next>=c*v and bounded-gradient upper bounds already occur in Zaheer2018AppendixB. General related Adam bounds with beta2>beta1^2 appear in Dereich et al2025; the displayed scalar specialization is not claimed historically new.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'protocol.json').write_text(json.dumps(dict(seed=2026110144,source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),bound_steps_per_run=10000),indent=2)+'\n')
    for n in ['analyze_yogi_structure.py','yogi_map.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    print({k:v for k,v in result.items() if k not in ['cold_bound_checks','boundary_initialization']},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
