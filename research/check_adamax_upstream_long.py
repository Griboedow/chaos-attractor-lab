"""Generate long states with actual PyTorch Adamax, not the normalized model."""
import argparse,hashlib,inspect,json,time
from pathlib import Path
import numpy as np
import torch
import torch.optim.adamax as upstream


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=np.load(args.source);ids=np.array([847,239,4099,4096]);parameters=source['parameters'][ids]
    alpha=.002;epsilon=1e-8;warm=100000;observed=200000
    tensors=[torch.tensor([alpha*.73],dtype=torch.float64,requires_grad=True) for _ in ids]
    groups=[dict(params=[tensor],betas=tuple(p[:2])) for tensor,p in zip(tensors,parameters)]
    optimizer=torch.optim.Adamax(groups,lr=alpha,eps=epsilon,weight_decay=0,foreach=False)
    upstream_file=Path(inspect.getsourcefile(upstream))
    protocol=dict(source=str(args.source.resolve()),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),
        torch_version=torch.__version__,torch_git_version=torch.version.git_version,source_file=str(upstream_file),
        source_file_sha256=hashlib.sha256(upstream_file.read_bytes()).hexdigest(),dtype='float64',alpha=alpha,epsilon=epsilon,
        initial_theta=alpha*.73,loss='Independent scalar quadratic theta**2/2 for each parameter group. Exact gradient theta assigned each step.',
        warmup=warm,observed=observed,bias_correction='Actual upstream optimizer at every step, including observation.',
        scope='Complete optimizer state trajectories. Statistical comparisons are necessary because chaotic trajectories separate under roundoff.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    data=np.empty((len(ids),observed+1,3));first20=np.empty((len(ids),20,3));before=time.perf_counter()
    for t in range(1,warm+observed+1):
        for tensor in tensors:tensor.grad=tensor.detach().clone()
        optimizer.step()
        if t<=20 or t>=warm:
            for i,tensor in enumerate(tensors):
                state=optimizer.state[tensor]
                z=[tensor.detach().item()/alpha,state['exp_avg'].item()/alpha,state['exp_inf'].item()/alpha]
                if t<=20:first20[i,t-1]=z
                if t>=warm:data[i,t-warm]=z
        if t%50000==0:print(t,'actual PyTorch steps',flush=True)
    assert np.isfinite(data).all()
    path=out/'histories.npz';np.savez_compressed(path,ids=ids,parameters=parameters,states=data,first20=first20)
    (out/'analysis.json').write_text(json.dumps(dict(completed=len(ids),declared=len(ids),steps_per_case=warm+observed,
        histories_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.perf_counter()-before),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
