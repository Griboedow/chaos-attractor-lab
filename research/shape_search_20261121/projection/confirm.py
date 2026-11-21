"""Long binary64 orbit, QR spectrum, independent starts and derivative controls."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from model import replay,derivative_audit
from screen import KINDS

LABELS={
    370279:'Ribbed fans',220329:'Folded whorl',301027:'Nested scrolls',
    516308:'Domed veil',428554:'Crested folds',72726:'Three curled bands',
    253746:'Pleated sail',454181:'Arched comb',418637:'Paired arches',
    449458:'Looped shell',351879:'Spiral wheel',354353:'Three folded sheets',
    198237:'Split crescents',127600:'Lenticular rim',
    1455969:'Tapered fern rows',1147837:'Curled fronds',1391489:'Paired combs',
    1282003:'Repeated folded crescents',1167605:'Crossed sails',
    1212574:'Coiled lobes',1006801:'Stepped curtain',1518001:'Rim with inner folds',
    1451361:'Parallel arches',1496296:'Coiled arc and fine train',
    1523981:'Pleated shield',1229838:'Offset joined sheets',
    346817:'Bent folding ramp',10238:'Fluted open shell',367791:'Folded wedge and tail',
    317370:'Terraced walls',303603:'Stepped arches',394727:'Four curled fins',
    175369:'Open scalloped bowl',475471:'Forked shell',367853:'Rows beneath an arch',
    186800:'Banded panel',121936:'Three tapered panels',119143:'Window grid',
    261538:'Folded canopy',45458:'Stacked pleats',
    1291273:'Layered cusps',1322937:'Slanted frames',1181714:'Staggered terraces',
    1079549:'Spatial lattice',1364145:'Crossed shelves',1161974:'Diagonal blades',
    2503686:'Coiled edge',2361186:'Paired spiral buds',2006357:'Crossed loops',
    2087905:'Folded arches',2304011:'Crescent with inner spiral',
    2008466:'Pleated dome',2189166:'Helical chain',2220593:'Folded wings',
    2425586:'Crossed rims',2508549:'Hooked funnel',2158223:'Stepped coil',
    2417165:'Strong two-rate blend',2427630:'Second high-rate blend',
    2039634:'Blended spatial fold',2028066:'Blended folded layers',
}


def save_json(path,doc):
    tmp=path.with_suffix(f'.{os.getpid()}.tmp.json')
    tmp.write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n',encoding='utf8')
    os.replace(tmp,path)


def periods(z):
    span=max(float(np.linalg.norm(np.ptp(z,axis=0))),1.)
    ref=z[-768:]
    residuals=[float(np.max(np.linalg.norm(ref-z[-768-k:-k],axis=1)))/span for k in range(1,513)]
    period=next((k+1 for k,v in enumerate(residuals) if v<1e-8),None)
    return period,min(residuals)


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    args.cache.mkdir(parents=True,exist_ok=True)
    for case in args.cases:
        record=args.out/f'case-{case}.json'
        if record.exists():
            continue
        with np.load(args.inspection/f'case-{case}.npz') as f:
            raw=f['data'].astype(float)
            kind=int(f['kind'])
            prior=f['rates'][:,0].tolist()
        linked=raw[1,3]>0
        smooth=raw[1,3]==2
        initial=raw[1,:3].copy()
        z,rates,mins,rad,branch=replay(raw,initial,args.burn,args.window,args.windows,args.tail,True)
        oracle,fd=derivative_audit(raw,z[np.linspace(0,len(z)-1,16,dtype=int)])
        cohort_seed=821711+case
        rand=np.random.default_rng(cohort_seed)
        starts=np.vstack([initial,rand.uniform(-2.,2.,size=(args.starts-1,3))])
        cohort=[]
        for state in starts:
            tail,rs,_,_,visits=replay(raw,state,args.start_burn,args.start_window,4,2048,args.cohort_qr)
            p,resid=periods(tail)
            cohort.append(dict(initial=state.tolist(),windows=rs[:,0].tolist(),period=p,minimum_period_residual=resid,branch_counts=visits.tolist()))
            if args.cohort_qr:
                cohort[-1]['qr_windows']=rs.tolist()
        p,resid=periods(z)
        cloud=args.cache/f'case-{case}-points.npz'
        tmp=cloud.with_suffix('.tmp.npz')
        np.savez_compressed(tmp,points=z.astype(np.float32),case=case,initial=initial)
        os.replace(tmp,cloud)
        stem='Linked-torus posterior reflection' if linked else 'Ring-torus posterior reflection'
        if smooth:
            stem='Blended linked-torus reflection'
        n_positive=sum(min(e['windows'])>.015 and e['period'] is None for e in cohort)
        strong=bool(np.min(rates[:,0])>.015 and p is None and n_positive>=int(np.ceil(len(cohort)/2)) and fd<1e-5)
        max_site=float(np.linalg.norm(raw[2:2+int(raw[0,3]),:3],axis=1).max())
        beta=float(raw[0,0])
        project_bound=1+float(raw[0,2])+int(linked)
        radius_bound=max(np.linalg.norm(initial),(beta*project_bound+abs(1-2*beta)*max_site)/(1-beta))
        result=dict(case=case,family=stem,label=LABELS.get(case,f'Case {case}'),kind=KINDS[kind],passes=strong,cloud=str(cloud.resolve()),cloud_sha256=hashlib.sha256(cloud.read_bytes()).hexdigest(),parameters=dict(major_radius=1.,tube_radius=float(raw[0,2]),beta=beta,inv_tau=float(raw[0,1]),tau=1/float(raw[0,1]),second_ring=bool(linked),sites=raw[2:2+int(raw[0,3]),:3].tolist()),initial=initial.tolist(),inspect_leading_windows=prior,burn=args.burn,windows=args.windows,window_length=args.window,points_saved=args.tail,scope='Numerical finite-time evidence; no trapping, topological equivalence or historical novelty certificate.',rates=rates.mean(axis=0).tolist(),qr_windows=rates.tolist(),minimum_leading_window=float(rates[:,0].min()),minimum_second_window=float(rates[:,1].min()),sampled_minimum_axis=float(mins[0]),sampled_minimum_centerline=float(mins[1]),sampled_minimum_branch_margin=float(mins[2]) if linked else None,branch_counts=branch.tolist(),max_norm=rad,analytic_radius_bound=radius_bound,period=p,minimum_period_residual=resid,nested_map_max_abs=oracle,central_difference_relative_jacobian=fd,cohort_seed=cohort_seed,cohort_burn=args.start_burn,cohort_window_length=args.start_window,cohort_positive=n_positive,cohort=cohort)
        if smooth:
            result['parameters']['mix_temperature']=float(raw[2,3])
        save_json(record,result)
        print(case,strong,np.round(rates.mean(axis=0),4),n_positive,'/',len(cohort),flush=True)
        manifest=[]
        for saved in sorted(args.out.glob('case-*.json')):
            r=json.loads(saved.read_text())
            if r['passes']:
                manifest.append(dict(case=r['case'],family=r['family'],label=r['label'],cloud=r['cloud'],path=r['cloud'],cloud_relative='projection/final/'+Path(r['cloud']).name,record=saved.name,rates=r['rates'],positive_starts=r['cohort_positive'],total_starts=len(r['cohort']),source='research/shape_search_20261121/projection/reproduce.py',points_saved=r['points_saved'],burn=r['burn'],observed=r['windows']*r['window_length'],cloud_sha256=r['cloud_sha256']))
        save_json(args.out/'manifest.json',manifest)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--inspection',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cases',type=int,nargs='+',required=True)
    p.add_argument('--burn',type=int,default=100000)
    p.add_argument('--windows',type=int,default=10)
    p.add_argument('--window',type=int,default=500000)
    p.add_argument('--tail',type=int,default=250000)
    p.add_argument('--starts',type=int,default=17)
    p.add_argument('--start-burn',type=int,default=25000)
    p.add_argument('--start-window',type=int,default=16384)
    p.add_argument('--cohort-qr',action='store_true')
    main(p.parse_args())
