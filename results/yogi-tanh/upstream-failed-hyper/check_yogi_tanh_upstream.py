"""Execute unmodified TensorFlow Addons Yogi(tanh) on a scalar quadratic."""
import os
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL','2')
os.environ.setdefault('TF_ENABLE_ONEDNN_OPTS','0')
import argparse,hashlib,inspect,json,platform,time
from pathlib import Path
import numpy as np
import tensorflow as tf
import tensorflow_addons as tfa


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    tf.config.threading.set_intra_op_parallelism_threads(2);tf.config.threading.set_inter_op_parallelism_threads(2)
    tf.keras.backend.set_floatx('float64')
    with np.load(args.source) as f:source={k:f[k] for k in f.files}
    sourcefile=Path(inspect.getfile(tfa.optimizers.Yogi));raw=sourcefile.read_bytes()
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=args.ids,
        tensorflow=tf.__version__,tensorflow_addons=tfa.__version__,numpy=np.__version__,python=platform.python_version(),
        installed_source=str(sourcefile),installed_source_sha256=hashlib.sha256(raw).hexdigest(),
        source_url='https://github.com/tensorflow/addons/blob/v0.22.0/tensorflow_addons/optimizers/yogi.py',
        dtype='float64, including hyperparameters; Keras default floatx set to float64. beta1 is a Python float so the upstream Python branch is legal inside a graph; other hyperparameters are float64 tensors. Unmodified optimizer implementation.',objective='theta^2/2',
        steps=args.steps,initialization='theta=.73*alpha,m=0,v=1e-6 from library slot creation; continuous library bias correction; no regularization.',
        scope='Actual library histories. Autonomous derivative analysis is performed separately after a long warmup; tiny parameter-rounding differences are recorded.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    rows=[]
    for ident in args.ids:
        before=time.perf_counter();b,c,e,k=source['parameters'][ident];alpha=float(np.sqrt(k/10));epsilon=float(alpha*e)
        variable=tf.Variable([.73*alpha],dtype=tf.float64)
        const=lambda x:tf.constant(x,dtype=tf.float64)
        optimizer=tfa.optimizers.Yogi(learning_rate=const(alpha),beta1=float(b),beta2=const(c),epsilon=const(epsilon),activation='tanh',initial_accumulator_value=1e-6)
        optimizer._create_all_weights([variable])
        initial=np.array([variable.numpy()[0]/alpha,optimizer.get_slot(variable,'m').numpy()[0]/alpha,optimizer.get_slot(variable,'v').numpy()[0]/alpha**2])
        hyper={name:float(optimizer._get_hyper(name,tf.float64).numpy()) for name in ['learning_rate','beta_1','beta_2','epsilon']}
        assert hyper==dict(learning_rate=alpha,beta_1=float(b),beta_2=float(c),epsilon=epsilon)
        # Independent physical-coordinate reference for the first20steps.
        physical=np.array([variable.numpy()[0],optimizer.get_slot(variable,'m').numpy()[0],optimizer.get_slot(variable,'v').numpy()[0]])
        control=[];actual=[]
        for t in range(1,21):
            theta,m,v=physical;M=b*m+(1-b)*theta;V=v+(1-c)*np.tanh(10*(theta*theta-v))*theta*theta
            physical=np.array([theta-alpha*np.sqrt(1-c**t)/(1-b**t)*M/(np.sqrt(V)+epsilon),M,V]);control.append(physical.copy())
            optimizer.apply_gradients([(variable.read_value(),variable)])
            actual.append([variable.numpy()[0],optimizer.get_slot(variable,'m').numpy()[0],optimizer.get_slot(variable,'v').numpy()[0]])
        control=np.array(control);actual=np.array(actual)
        first20_scaled=float(np.max(abs(control-actual)/(1+abs(control))))
        assert first20_scaled<1e-9
        @tf.function(autograph=False)
        def integrate(count):
            buffer=tf.TensorArray(tf.float64,size=count,clear_after_read=False)
            def body(t,history):
                optimizer.apply_gradients([(variable.read_value(),variable)])
                state=tf.stack([variable[0],optimizer.get_slot(variable,'m')[0],optimizer.get_slot(variable,'v')[0]])
                return t+1,history.write(t,state)
            _,buffer=tf.while_loop(lambda t,_:t<count,body,(tf.constant(0),buffer),parallel_iterations=1)
            return buffer.stack()
        remaining=integrate(tf.constant(args.steps-20)).numpy();states=np.r_[actual,remaining]/np.array([alpha,alpha,alpha*alpha])
        assert np.isfinite(states).all()
        effective=np.array([b,c,epsilon/alpha,10*alpha*alpha])
        np.savez_compressed(out/f'case-{ident}.npz',states=states,parameters=effective,scout_parameters=source['parameters'][ident],initial=initial,first20_physical_actual=actual,first20_physical_reference=control)
        row=dict(id=ident,parameters=effective.tolist(),scout_parameters=source['parameters'][ident].tolist(),alpha=alpha,epsilon=epsilon,hyperparameters=hyper,
            initial=initial.tolist(),first20_max_scaled_difference=first20_scaled,final=states[-1].tolist(),seconds=time.perf_counter()-before)
        rows.append(row);(out/'analysis.json').write_text(json.dumps(dict(runs=rows,completed=len(rows),declared=len(args.ids)),indent=2)+'\n');print(row,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--steps',type=int,default=300000);main(p.parse_args())
