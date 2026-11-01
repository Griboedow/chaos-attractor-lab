"""Declared long checks of visually selected spatial point-set RAAR candidates."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import set_num_threads
from follow_raar_constellations import follow


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);set_num_threads(6)
    with np.load(args.source) as f:
        ids=np.array(args.ids);p=f['params'][ids];points=f['points'][ids];counts=f['counts'][ids];kinds=f['kinds'][ids]
    protocol=dict(source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),ids=ids.tolist(),selection='Visual inspection of all 72 saved scout-follow-up portraits: fragmented columns, stepped fan, lattice sheets, parallel pillars, T-shaped folds and curved clusters.',warm=100000,observed=args.count,initial=[.317,-.211,.419],exact_model='Each beta and each supplied point coordinate is its exact binary64 value. Constructed point sets are numerical rotated/scaled/translated polyhedron vertices, not an asserted exact algebraic symmetry.')
    (args.out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    for name in [Path(__file__).name,'follow_raar_constellations.py','search_raar_constellations.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    begin=time.time();rates,states,status,minima=follow(p,points,counts,100000,args.count);rows=[]
    for j,i in enumerate(ids):
        np.savez_compressed(args.out/f'case-{i}.npz',states=states[j],params=p[j],points=points[j,:counts[j]],rates=rates[j],status=status[j],minima=minima[j],kind=kinds[j])
        rows.append(dict(id=int(i),status=status[j].tolist(),spectrum_xy_then_z=rates[j].mean(0).tolist(),minimum_leading_xy_window=float(rates[j,:,0].min()),minimum_gap=float(minima[j,0]),minimum_axis=float(minima[j,1]),beta_hex=float(p[j,0]).hex(),point_hex=[[float(x).hex() for x in v] for v in points[j,:counts[j]]]))
    report=dict(protocol,cases=rows,seconds=time.time()-begin)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps({**report,'cases':[{k:v for k,v in row.items() if k!='point_hex'} for row in rows]}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--ids',type=int,nargs='+',required=True);p.add_argument('--count',type=int,default=5000000);main(p.parse_args())
