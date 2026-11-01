"""Compare actual-library phase increments with the known near-half-turn picture."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);rows=[]
    for ident in [335,448,1890,2006]:
        source=args.upstream/f'case-{ident}.npz'
        with np.load(source) as f:s=f['states'];params=f['params']
        zeta=params[0];omega=np.sqrt(1-zeta*zeta);h=s[100001:,4];z=h*(-zeta+1j*omega)
        R=1+z+z*z/2+z**3/6+z**4/24+z**5/120+z**6/600
        phase=np.angle(R);distance=np.pi-np.abs(phase)
        old=s[100000:-1,0]-1j*(s[100000:-1,1]+zeta*s[100000:-1,0])/omega
        new=s[100001:,0]-1j*(s[100001:,1]+zeta*s[100001:,0])/omega
        error=np.abs(new-R*old)/np.maximum(1,np.abs(new))
        row=dict(id=ident,points=len(phase),discarded_accepted_steps=100000,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            distance_from_half_turn_quantiles=np.quantile(distance,[0,.1,.5,.9,1]).tolist(),quantile_levels=[0,.1,.5,.9,1],
            fraction_negative_real_R=float(np.mean(R.real<0)),mean_log_modulus=float(np.mean(np.log(np.abs(R)))),max_complex_step_residual=float(error.max()))
        np.savez_compressed(args.out/f'case-{ident}.npz',phase=phase,modulus=np.abs(R),previous_accepted_step=h,params=params)
        rows.append(row)
    report=dict(runs=rows,scope='Finite phase statistics of actual-library accepted steps. Near-half-turn spurious oscillation is credited to earlier work, not proposed as a new mechanism.',prior_source='https://math.mcgill.ca/humphries/research/papers/psthetarep.pdf',prior_section='5.2, report page 20')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--upstream',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
