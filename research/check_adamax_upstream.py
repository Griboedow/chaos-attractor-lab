"""Compare normalized cold starts with the actual installed PyTorch Adamax."""
import argparse,hashlib,inspect,json
from pathlib import Path
import numpy as np
import torch
import torch.optim.adamax as upstream


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    alpha=.002;epsilon=1e-8;records=[];all_actual=[];all_reference=[]
    for b,c in [(.9,.999),(.9,.99),(.52,.15),(.15,.9)]:
        for initial in [.73,-.73,3.]:
            parameter=torch.tensor([alpha*initial],dtype=torch.float64,requires_grad=True)
            optimizer=torch.optim.Adamax([parameter],lr=alpha,betas=(b,c),eps=epsilon,weight_decay=0,foreach=False)
            x,m,u=initial,0.,0.;actual=[];reference=[]
            for t in range(1,51):
                parameter.grad=parameter.detach().clone();optimizer.step();state=optimizer.state[parameter]
                actual.append([float(parameter.detach()[0]/alpha),float(state['exp_avg'][0]/alpha),float(state['exp_inf'][0]/alpha)])
                m=b*m+(1-b)*x;u=max(c*u,abs(x)+epsilon/alpha);x-=m/((1-b**t)*u)
                reference.append([x,m,u])
            actual=np.array(actual);reference=np.array(reference)
            error=float(np.max(abs(actual-reference)/(1+abs(reference))));assert error<1e-8
            records.append(dict(beta1=b,beta2=c,normalized_initial=initial,maximum_scaled_state_error=error))
            all_actual.append(actual);all_reference.append(reference)
    source=Path(inspect.getsourcefile(upstream));np.savez_compressed(out/'histories.npz',actual=all_actual,reference=all_reference)
    (out/'analysis.json').write_text(json.dumps(dict(torch_version=torch.__version__,torch_git_version=torch.version.git_version,
        source_file=str(source),source_file_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dtype='float64',alpha=alpha,epsilon=epsilon,steps_per_run=50,runs=records,
        maximum_scaled_state_error=max(r['maximum_scaled_state_error'] for r in records),
        scope='Actual continuously bias-corrected upstream state updates versus an independent scalar normalization. Short controls, not a long trajectory or chaos claim.'),indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print('Passed',len(records),'cold starts; max scaled error',max(r['maximum_scaled_state_error'] for r in records),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
