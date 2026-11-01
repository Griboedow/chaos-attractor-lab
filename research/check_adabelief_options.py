"""Upstream rectification/AMSGrad controls, plus exact spectral identities."""
import argparse,importlib.util,json
from pathlib import Path
import numpy as np


def main():
    import torch
    p=argparse.ArgumentParser();p.add_argument('--author-source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    spec=importlib.util.spec_from_file_location('author_options',a.author_source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    reports=[]
    for rectify,amsgrad in [(True,False),(False,True)]:
        theta=torch.tensor([.0073],dtype=torch.float64,requires_grad=True)
        opt=module.AdaBelief([theta],lr=.01,betas=(.52,.15),eps=1e-16,
           weight_decay=0,weight_decouple=False,rectify=rectify,amsgrad=amsgrad,print_change_log=False)
        tail=[]
        for t in range(20000):
            theta.grad=theta.detach().clone();opt.step()
            if t>=19000:tail.append(theta.item())
        reports.append(dict(rectify=rectify,amsgrad=amsgrad,iterations=20000,
                            final_parameter=theta.item(),maximum_last_1000=float(np.max(abs(np.array(tail))))))
    # An independent exact symbolic determinant check in the (x,m',s') factors.
    import sympy as sp
    b,c,u,v=sp.symbols('b c u v')
    factor1=sp.Matrix([[1,0,0],[1-b,b,0],[u,-u,c]])
    factor2=sp.Matrix([[1,-v,sp.Symbol('w')],[0,1,0],[0,0,1]])
    identity=sp.factor((factor2*factor1).det()-b*c)
    assert identity==0
    report=dict(runs=reports,symbolic_determinant_residual=str(identity),torch=torch.__version__,sympy=sp.__version__)
    (a.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
