"""Read-only behavioral audit of the public 2025 pneumatic-limb helper.

The unmodified author file is supplied explicitly from its public replication
archive. This script does not execute the simulation notebooks or write there.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    spec=importlib.util.spec_from_file_location('pneumatic_author_helper',args.helper)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    n=17;theta=np.linspace(-np.pi/2,np.pi/2,15)
    positions=np.column_stack([3.75+3.75*np.sin(theta),3.75*np.cos(theta)])
    positions=np.vstack([[0.,-1.],positions,[7.5,-1.]])
    beta=np.array([1.2,.1,.2]);knots=np.array([-2*np.pi,-2*np.pi/3,-.3,0,.3,2*np.pi/3,2*np.pi])
    values=np.array([-beta[0],-beta[1],-beta[2],0,beta[2],beta[1],beta[0]])
    curve=PchipInterpolator(knots,values)
    state=np.r_[positions.ravel(),np.zeros(2*n)];memory=np.zeros((n,2))
    first_input=state.copy();first=helper.ode_system(601.,first_input,memory,n,1.,.4,beta,1.8,curve,.01,5.,1.5)
    first_memory=memory.copy()
    second_input=state.copy();second=helper.ode_system(601.,second_input,memory,n,1.,.4,beta,1.8,curve,.01,5.,1.5)
    active=np.zeros((n,2));active[:,0]=600.;active[:,1]=1.
    # t<600 suppresses metadata updates without changing the constitutive calls.
    # Instead use t=601 with separately copied memories to compare identical states.
    force0=helper.compute_torsional_spring_forces(positions.copy(),601.,active.copy(),beta,1.8,curve,n,5.,1.5)
    shifted=positions.copy();shifted[:,0]+=10.
    force1=helper.compute_torsional_spring_forces(shifted,601.,active.copy(),beta,1.8,curve,n,5.,1.5)
    # Change only the declared sigmoid-midpoint argument, keeping geometry/memory.
    force_midpoint=helper.compute_torsional_spring_forces(positions.copy(),601.,active.copy(),beta,1.8,curve,n,5.,100.)
    record=dict(source='https://zenodo.org/records/15399470',file=args.helper.name,
        source_sha256=hashlib.sha256(args.helper.read_bytes()).hexdigest(),
        probe='17-node source semicircle, zero velocities; t601, multiplier1.8,k5,declared midpoint1.5.',
        changed_input_state_max=float(np.max(abs(first_input-state))),
        first_call_changed_memory_rows=np.flatnonzero(np.any(first_memory!=0,axis=1)).tolist(),
        repeat_same_t_and_explicit_state_rhs_difference=float(np.max(abs(first-second))),
        translated_same_shape_torsional_force_difference=float(np.max(abs(force0-force1))),
        change_midpoint_argument_1_5_to_100_force_difference=float(np.max(abs(force0-force_midpoint))),
        paper_pchip_positive_values_at_0_03_2pi3_2pi=[0,1.2,.1,.2],
        notebook_pchip_positive_values_at_0_03_2pi3_2pi=values[3:].tolist(),
        conclusions=[
            'The callback mutates explicit solver state and an external discrete-memory array.',
            'Repeated RHS evaluation at identical t and explicit y depends on hidden mutable history; rejected solver stages are not rolled back by this helper.',
            'The sigmoid-midpoint argument is overwritten by a neighboring x coordinate; a rigid translation changes the torsional forces in the probe.',
            'The paper and notebook order the first and last nonzero positive PCHIP control values differently.',
            'These are implementation/specification gates, not evidence against the experimental oscillations and not a no-chaos theorem.'])
    assert record['changed_input_state_max']>0
    assert record['repeat_same_t_and_explicit_state_rhs_difference']>1e-12
    assert record['translated_same_shape_torsional_force_difference']>1e-5
    assert record['change_midpoint_argument_1_5_to_100_force_difference']==0
    (out/'audit.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--helper',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
