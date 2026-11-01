"""Independent high-precision and actual-library controls for the AdaDelta audit."""
import argparse,hashlib,inspect,json,platform
from pathlib import Path
import mpmath as mp
import numpy as np
import torch
import torch.optim.adadelta as implementation


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);records=[]
    for digits in [70,90]:
        mp.mp.dps=digits;alpha=mp.mpf(1)/10;beta=mp.mpf(9)/10;E=mp.mpf(1)/10**8
        x=mp.mpf(1);d=mp.mpf(400);v=mp.mpf(1);samples=[]
        for t in range(1,5001):
            V=beta*v+(1-beta)*x*x
            change=x*mp.sqrt((d+E)/(V+E))
            D=beta*d+(1-beta)*change*change
            x=x-alpha*change;d=D;v=V
            if t in [1,2,20,5000]:samples.append(dict(step=t,state=[mp.nstr(q,digits) for q in [x,d,v]]))
        records.append(dict(digits=digits,parameters='Exact rationals alpha=1/10,beta=9/10,curvature=1,epsilon_squared=1/100000000.',samples=samples))
    theta=torch.tensor([1.],dtype=torch.float64,requires_grad=True)
    optimizer=torch.optim.Adadelta([theta],lr=.1,rho=.9,eps=1e-8,foreach=False)
    theta.grad=torch.zeros_like(theta);optimizer.step()
    state=optimizer.state[theta];state['step'].zero_();state['square_avg'].fill_(1.);state['acc_delta'].fill_(400.)
    history=[]
    for _ in range(5000):
        theta.grad=theta.detach().clone();optimizer.step()
        history.append([float(theta.detach()[0]),float(state['acc_delta'][0]),float(state['square_avg'][0])])
    history=np.array(history);np.savez_compressed(out/'actual-torch.npz',states=history)
    independent=np.array([float(q) for q in records[-1]['samples'][-1]['state']])
    assert np.max(abs(history[-1]-independent))<1e-8
    file=Path(inspect.getsourcefile(implementation));function=inspect.getsource(implementation._single_tensor_adadelta)
    (out/'inspected-function.py').write_text(function,encoding='utf-8')
    result=dict(high_precision=records,torch=dict(version=torch.__version__,dtype=str(theta.dtype),platform=platform.platform(),
        parameters=dict(lr=.1,rho=.9,eps=1e-8,foreach=False),initial_state=[1.,400.,1.],
        initialization='Actual optimizer initialized with a zero-gradient call, then its step counter and both memory buffers explicitly set to the stated source state before the observed 5000 gradient updates.',
        state_columns=['theta','acc_delta','square_avg'],state_after_two=history[1].tolist(),state_after5000=history[-1].tolist(),
        endpoint_max_difference_from90digit_exact_rationals=float(np.max(abs(history[-1]-independent))),
        installed_source_file=str(file),installed_source_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        inspected_function_sha256=hashlib.sha256(function.encode()).hexdigest()),
        scope='Finite implementation and arithmetic controls. Exact period-two exclusion is proved separately; finite traces do not prove global convergence.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(result['torch'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
