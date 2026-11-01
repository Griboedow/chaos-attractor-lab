"""Wide-domain grid closure with explicit monitoring of normalization zeros.

This is a sampled proposal. Reaching a pole cell rejects this proposed smooth
trap, but is not an orbit reaching the pole or a proof against other traps.
"""
import argparse,hashlib,json,time
from pathlib import Path
import itertools
import numpy as np
from scipy.optimize import brentq
from propose_product_nonuniform_grid import advance,occupied


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    for name in [Path(__file__).name,'propose_product_nonuniform_grid.py','search_product_alphabet.py','follow_spherical_constellations.py','search_soft_sphere.py','search_soft_constellations.py','search_raar_constellations.py','search_spherical_constellations.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:z=f['states'];beta,tau=f['params'];levels=f['levels'];order=int(f['order'])
    lists=[];roots=[]
    for i in range(3):
        c,a=levels[i,:2].mean(),np.diff(levels[i,:2])[0]/2
        def coordinate(x):return x+c+a*np.tanh(2*a*(x-c)/tau)
        low,high=-args.halfwidth,args.halfwidth
        ys=np.linspace(coordinate(low),coordinate(high),int(np.ceil((coordinate(high)-coordinate(low))/args.step))+1)
        lists.append(np.array([low]+[brentq(lambda x:coordinate(x)-y,low,high,xtol=1e-15) for y in ys[1:-1]]+[high]))
        def v(x):return 2*(c+a*np.tanh(2*a*(x-c)/tau))-x
        xx=np.linspace(low,high,10001);yy=v(xx);rr=[]
        for j in np.flatnonzero(yy[:-1]*yy[1:]<0):rr.append(brentq(v,xx[j],xx[j+1],xtol=1e-15))
        roots.append(rr)
    poles=np.array(list(itertools.product(*roots)));shape=np.array([len(v)-1 for v in lists],dtype=np.int64)
    assert np.prod(shape)<100000000
    edges=np.full((3,max(shape)+1),np.inf)
    for i in range(3):edges[i,:len(lists[i])]=lists[i]
    locate=lambda states:np.column_stack([np.searchsorted(lists[i],states[:,i],side='right')-1 for i in range(3)])
    ijk=locate(z);frontier=np.unique((ijk[:,0]*shape[1]+ijk[:,1])*shape[2]+ijk[:,2])
    pole_ijk=locate(poles);pole_keys=(pole_ijk[:,0]*shape[1]+pole_ijk[:,1])*shape[2]+pole_ijk[:,2]
    mask=np.zeros(np.prod(shape),bool);mask[frontier]=True;initial=len(frontier);history=[]
    save=lambda name,v:(args.out/name).write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(shape=shape.tolist(),poles=len(poles),initial=initial)),flush=True)
    for generation in range(args.generations):
        pending=frontier.copy()
        new,flag,pos,key,minr,maxrem=advance(frontier,mask,edges,shape,beta,tau,levels,order,2.,1e-8,args.budget)
        reached=np.flatnonzero(mask[pole_keys])
        row=dict(generation=generation,sources=len(frontier),processed=int(pos),added=len(new),status=int(flag),minimum_sampled_radius=float(minr),maximum_sampled_remainder=float(maxrem),pole_cells_reached=reached.tolist(),seconds=time.time()-begin)
        history.append(row);save('progress.json',history);print(json.dumps(row),flush=True)
        frontier=new
        if flag or len(reached) or not len(new):break
    keys=occupied(mask)
    np.savez_compressed(args.out/'proposal.npz',edges=edges,shape=shape,occupied=keys,pending=pending,new=frontier,poles=poles,pole_keys=pole_keys,params=[beta,tau],levels=levels,order=order)
    result=dict(status='pole_cell_reached' if len(reached) else 'numerically_closed' if not flag and not len(frontier) else 'not_closed',
        source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),halfwidth=args.halfwidth,step=args.step,shape=shape.tolist(),initial=initial,occupied=len(keys),poles=poles.tolist(),reached_pole_indices=reached.tolist(),history=history,seconds=time.time()-begin,scope=__doc__)
    save('report.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['source','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--step',type=float,default=.02);p.add_argument('--halfwidth',type=float,default=2.2)
    p.add_argument('--budget',type=int,default=20000000);p.add_argument('--generations',type=int,default=20)
    main(p.parse_args())
