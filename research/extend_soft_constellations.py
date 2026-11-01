"""Long fresh-start checks of the declared Gaussian-posterior candidates."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from follow_soft_constellations import follow


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    with np.load(args.source) as f:
        ids=np.array(args.ids);params=f['params'][ids];points=f['points'][ids];counts=f['counts'][ids];kinds=f['kinds'][ids];base_ids=f['base_ids'][ids];relative=f['relative_tau'][ids]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),selection=args.selection,warm=100000,observed=args.count,windows=20,initial=[.317,-.211,.419],scope='Exact binary64 input parameters define the real model. Finite-time numerical evidence, not asymptotic attraction or novelty.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
    for name in [Path(__file__).name,'follow_soft_constellations.py','search_soft_constellations.py','search_raar_constellations.py','follow_spherical_constellations.py','search_spherical_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    begin=time.time();rates,states,status,minima,stats=follow(params,points,counts,100000,args.count);rows=[]
    for j,ident in enumerate(ids):
        np.savez_compressed(args.out/f'case-{ident}.npz',params=params[j],points=points[j,:counts[j]],kind=kinds[j],base_id=base_ids[j],relative_tau=relative[j],states=states[j],rates=rates[j],status=status[j],minima=minima[j],soft_stats=stats[j])
        rows.append(dict(id=int(ident),base_id=int(base_ids[j]),beta=float(params[j,0]),tau=float(params[j,1]),status=status[j].tolist(),spectrum_qr_order=rates[j].mean(0).tolist(),minimum_leading_window=float(rates[j,:,0].min()),minimum_second_window=float(rates[j,:,1].min()),minimum_radius=float(minima[j,0]),soft_stats=stats[j].tolist(),parameter_hex=[float(x).hex() for x in params[j]],point_hex=[[float(x).hex() for x in v] for v in points[j,:counts[j]]]))
    report=dict(protocol,cases=rows,seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({**report,'cases':[{k:v for k,v in row.items() if k!='point_hex'} for row in rows]}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--count',type=int,default=5000000);p.add_argument('--selection',required=True);main(p.parse_args())
