"""A bounded comparison set, not an exhaustive census of strange attractors."""
from pathlib import Path
import os,json,time,concurrent.futures,argparse
os.environ.setdefault('NUMBA_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import dysts


def simulate(name,BASE):
    import dysts.flows as flows
    model=getattr(flows,name)()
    t=time.time()
    try:
        z=model.make_trajectory(10000,method='DOP853',rtol=1e-9,atol=1e-11,
                               postprocess=False,standardize=False,pts_per_period=100)
        if z is None or z.shape!=(10000,3) or not np.isfinite(z).all():
            raise ValueError('Incomplete or nonfinite trajectory')
        np.savez_compressed(BASE/'trajectories'/f'{name}.npz',points=z[2000:])
        return {'name':name,'status':'ok','seconds':time.time()-t,'period':float(model.period)}
    except Exception as e:
        return {'name':name,'status':'failed','error':str(e),'seconds':time.time()-t}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    BASE=parser.parse_args().out
    source=Path(dysts.__file__).parent/'data/chaotic_attractors.json'
    data=json.loads(source.read_text())
    selected=[name for name,v in data.items() if v['embedding_dimension']==3 and not v['delay'] and not v['nonautonomous'] and not v['unbounded_indices']]
    BASE.mkdir(parents=True,exist_ok=True);(BASE/'trajectories').mkdir(exist_ok=True)
    records=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        jobs={pool.submit(simulate,n,BASE):n for n in selected}
        for future in concurrent.futures.as_completed(jobs):
            row=future.result();records.append(row)
            print(json.dumps(row),flush=True)
            (BASE/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
    info={'source_commit':'2a03f1ae7b0680b0470458783dcb4664660e131a','selected':selected,
          'excluded':{n:{k:v[k] for k in ['embedding_dimension','delay','nonautonomous','unbounded_indices']} for n,v in data.items() if n not in selected},
          'n':10000,'discarded':2000,'points_per_period':100,'rtol':1e-9,'atol':1e-11,'method':'DOP853',
          'postprocess':False,'standardize':False,'initial_conditions':'catalog defaults',
          'metadata':{n:data[n] for n in selected}}
    (BASE/'protocol.json').write_text(json.dumps(info,indent=2)+'\n')

if __name__=='__main__':main()
