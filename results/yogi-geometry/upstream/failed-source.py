"""Execute actual Optax Yogi on the scalar quadratic, in float64."""
import argparse,hashlib,inspect,json,platform,time
from pathlib import Path
import numpy as np
import jax
jax.config.update('jax_enable_x64',True)
import jax.numpy as jnp
import optax
from optax._src import transform


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    d=np.load(args.source);ids=np.array(args.ids);p=d['parameters'][ids];alpha=.01
    b=jnp.array(p[:,0]);c=jnp.array(p[:,1]);eps=jnp.array(p[:,2]*alpha)
    solver=optax.yogi(alpha,b1=b,b2=c,eps=eps)
    theta=jnp.full(len(ids),.73*alpha,dtype=jnp.float64);state=solver.init(theta)
    def update(carry,unused):
        theta,state=carry
        updates,state=solver.update(theta,state,theta)
        theta=optax.apply_updates(theta,updates)
        values=jnp.stack([theta/alpha,state[0].mu/alpha,state[0].nu/alpha**2],axis=-1)
        return (theta,state),values
    (short_theta,short_state),short=jax.jit(lambda th,st:jax.lax.scan(update,(th,st),None,length=20))(theta,state)
    short=np.asarray(short);reference=[];z=np.tile([.73,0.,.01],(len(ids),1))
    for t in range(1,21):
        x,m,v=z.T;mn=p[:,0]*m+(1-p[:,0])*x;vn=v-(1-p[:,1])*np.sign(v-x*x)*x*x
        xn=x-(mn/(1-p[:,0]**t))/(np.sqrt(vn/(1-p[:,1]**t))+p[:,2])
        z=np.column_stack([xn,mn,vn]);reference.append(z.copy())
    errors=np.max(abs(short-np.array(reference)),axis=(0,2));assert errors.max()<1e-7
    source=Path(inspect.getsourcefile(transform));source_text=inspect.getsource(transform.scale_by_yogi)
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),parameters=p.tolist(),
        optax_version=optax.__version__,jax_version=jax.__version__,numpy_version=np.__version__,platform=platform.platform(),
        upstream_file_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),function_sha256=hashlib.sha256(source_text.encode()).hexdigest(),
        alpha=alpha,theta0=.73*alpha,initial_second_moment=1e-6,normalized_initial_v=.01,
        epsilon=(alpha*p[:,2]).tolist(),warm=100000,observed=200000,
        implementation='Actual optax.yogi; jax_enable_x64=True; CPU; continuously bias-corrected throughout; gradient exactly theta for f(theta)=theta^2/2.',
        first20_max_scaled_state_error=errors.tolist(),
        scope='Library control uses its default initial accumulator, which differs from both autonomous scout starts. Basin outcomes may differ.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    (out/'check_yogi_upstream.py').write_bytes(Path(__file__).read_bytes())
    began=time.perf_counter()
    def run(th,st):
        def warm(k,carry):return update(carry,None)[0]
        carry=jax.lax.fori_loop(0,100000,warm,(th,st))
        endpoint=jnp.stack([carry[0]/alpha,carry[1][0].mu/alpha,carry[1][0].nu/alpha**2],axis=-1)
        last,history=jax.lax.scan(update,carry,None,length=200000)
        return endpoint,history
    endpoint,history=jax.jit(run)(theta,state);endpoint=np.asarray(endpoint);history=np.asarray(history)
    assert np.isfinite(history).all() and history.dtype==np.float64
    np.savez_compressed(out/'histories.npz',ids=ids,parameters=p,states=history,warm_endpoint=endpoint,first20=short)
    result=dict(completed=len(ids),first20_max_error=float(errors.max()),late_extents=np.ptp(history[-10000:],axis=0).tolist(),
        wall_seconds=time.perf_counter()-began,scope='Long history is actual library output; no pathwise equality of chaotic trajectories is expected.')
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',required=True);main(p.parse_args())
