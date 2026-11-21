"""Recheck frozen entropy packets with both interval formulas and locate their orbits."""
import argparse
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.spatial import cKDTree
import snapback
from analytic_bounds import iterate as analytic
from certify_cycles import iterate as jets
from cycle_seeds import orbit_eval
from reproduce import raw_parameters


def main(args):
    mp.mp.dps=mp.iv.dps=160
    rows=[]
    for p in sorted((args.source/'snapback').glob('*.json')):
        if p.stem.startswith('proposal') or p.stem=='audit':continue
        doc=json.loads(p.read_text())
        if 'blocks' not in doc:continue
        case=doc['case'];rec=json.loads((args.source/f'case-{case}.json').read_text())
        checks={}
        for name,model in [('jet',jets),('analytic',analytic)]:
            snapback.enclosure_model=model
            proof=snapback.close_branches(doc,rec)
            if proof['uniform_return_steps']!=doc['proof']['uniform_return_steps']:
                raise ArithmeticError('Different entropy denominators')
            checks[name]={k:proof[k] for k in ('local_contraction_upper_display','largest_inverse_box_fraction_display','preimage_chart_distance_upper_display','local_branch_lipschitz_upper_display','excursion_lipschitz_upper_display')}
        separated_sinks=[]
        sink_packet=args.source/'sinks'/f'cycles-{case}.json'
        if sink_packet.exists():
            _,basis_inverse=snapback.exact_basis(doc['basis'])
            root=snapback.block([mp.iv.mpf(v) for v in doc['fixed']['center']],mp.iv.mpf(doc['fixed']['root_radius']))
            small=mp.iv.mpf(doc['small_chart_radius'])
            for j,sc in enumerate(json.loads(sink_packet.read_text())['cycles']):
                if sc['proof']['unstable_dimension']!=0:continue
                bounds={}
                for name,model in [('jet',jets),('analytic',analytic)]:
                    snapback.verify(sc,rec,sc['proof']['radius'],model)
                    seed=snapback.block([mp.iv.mpf(v) for v in sc['center']],mp.iv.mpf(sc['proof']['radius']))
                    _,_,orbit=model(seed,sc['period'],rec,mp.iv)
                    distances=[max(abs(v).a for v in snapback.matvec(basis_inverse,[u-w for u,w in zip(x,root)])) for x in orbit]
                    minimum=min(distances)
                    if not minimum>small:raise ArithmeticError('Sink intersects the entropy chart')
                    bounds[name]=float(minimum.a)
                separated_sinks.append(dict(index=j,least_period=sc['period'],distance_lower_display=bounds))
        raw=raw_parameters(rec)
        _,_,periodic=orbit_eval(np.array(doc['fixed']['center'],float),raw,doc['fixed']['period'])
        excursion=[]
        for block in doc['blocks']:
            _,_,arc=orbit_eval(np.array(block['center'],float),raw,block['period'])
            excursion.extend(arc)
        with np.load(rec['cloud']) as f:z=f['points'].astype(float)
        rad=np.sqrt(np.var(z,axis=0).sum());half=len(z)//2
        tree=cKDTree(z)
        periodic_max=float(tree.query(periodic)[0].max()/rad)
        excursion_dist=tree.query(excursion)[0]/rad
        temporal=cKDTree(z[:half]).query(z[half::11])[0]/rad
        rows.append(dict(case=case,packet=p.name,period=doc['fixed']['period'],blocks=len(doc['blocks']),return_steps=doc['proof']['uniform_return_steps'],W=doc['local']['radius'],V=doc['small_chart_radius'],checks=checks,disjoint_certified_sinks=separated_sinks,sample_proximity=dict(units='RMS radius',periodic_orbit_max=periodic_max,excursion_max=float(excursion_dist.max()),excursion_median=float(np.median(excursion_dist)),late_to_early_nearest_q50_q99=np.quantile(temporal,[.5,.99]).tolist())))
    out=args.source/'snapback'/'audit.json'
    out.write_text(json.dumps(dict(precision_decimal=160,mpmath=mp.__version__,cases=rows),indent=2)+'\n',encoding='utf8')
    for row in rows:print(row['case'],row['return_steps'],row['sample_proximity'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True,type=Path)
    main(p.parse_args())
