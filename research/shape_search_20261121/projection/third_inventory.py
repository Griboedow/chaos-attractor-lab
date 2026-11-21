"""Build the three-rate evidence index with portable clouds and empirical exclusions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

LABELS={
    3008791:'Cupped fold',3011275:'Flared sail',3014945:'Folded sail with a second basin',
    3061571:'Raised sail',3064303:'Curved sail',3092055:'Separated patches',
    3143582:'Triangular veil',3231505:'Folded sail with a period-four basin',
    3265947:'Layered sail',3275794:'Folded sail with a period-two basin',
    3295921:'Curled bow',3333965:'Slanted veil',3336225:'Broad folded sail',3386440:'Horned rim',
    6044998:'Long pointed veil',6071984:'Arched veil',6174331:'Double edge',
    6227178:'Curled cup',6263583:'Two pointed bands',6268324:'Triangular canopy',
    6271020:'Split lower sail',6438193:'Arched canopy',6560098:'Diagonal folded sail',
    6583985:'Horned cup',6616162:'Forked sail',6653105:'Broad creased veil',
    6679894:'Long triangular sail',6746662:'Crescent-edged sail',6758635:'Deep canopy',
    6802681:'Layered triangular sail',6858107:'Hooked veil',7033904:'Folded double edge',
}


def read(path):return json.loads(path.read_text())


def main(args):
    data=[];rows=[]
    certificates={}
    audit=args.source/'snapback'/'audit.json'
    if audit.exists():certificates={v['case']:v for v in read(audit)['cases']}
    for p in sorted(args.source.glob('case-*.json')):
        rec=read(p);case=rec['case']
        rec['label']=LABELS.get(case,'Folded sail')
        p.write_text(json.dumps(rec,indent=2)+'\n',encoding='utf8')
        cloud=Path(rec['cloud'])
        assert hashlib.sha256(cloud.read_bytes()).hexdigest()==rec['cloud_sha256']
        with np.load(cloud) as f:
            z=f['points'].astype(float)
        assert z.shape==(250000,3) and np.isfinite(z).all()
        eig=np.linalg.eigvalsh(np.cov(z.T));eig/=eig.sum()
        audits=[read(path) for path in sorted((args.source/'precision').glob(f'mpfr-three-{case}-*.json'))]
        assert any(a['bits']==192 for a in audits)
        precision_min=min(v[2] for a in audits for v in a['windows'])
        empirical=rec['three_positive_passes'] and precision_min>.001
        entropy=certificates.get(case)
        cycles=args.source/'cycles'/f'cycles-{case}.json'
        sinks=args.source/'sinks'/f'cycles-{case}.json'
        row=dict(case=case,label=rec['label'],parent=rec['parent_case'],phase=rec['third_growth_stream'],sites=len(rec['parameters']['sites']),lambda1=rec['rates'][0],lambda2=rec['rates'][1],lambda3=rec['rates'][2],min_third=rec['minimum_third_window'],three_positive_starts=rec['cohort_three_positive'],mpfr_bits='/'.join(str(a['bits']) for a in audits),mpfr_min_third=precision_min,pca_smallest=eig[0],pca_middle=eig[1],pca_largest=eig[2],retained=empirical,entropy_lower_denominator=entropy['return_steps'] if entropy else None,certified_cycles=len(read(cycles)['cycles']) if cycles.exists() else 0,certified_sinks=len(read(sinks)['cycles']) if sinks.exists() else 0)
        rows.append(row)
        if not empirical:continue
        data.append(dict(case=case,family=rec['family'],label=rec['label'],study='three-direction growth',cloud=str(cloud),path=str(cloud),cloud_relative='projection/final/'+cloud.name,record=p.name,rates=rec['rates'],positive_starts=rec['cohort_positive'],total_starts=len(rec['cohort']),three_positive_starts=rec['cohort_three_positive'],source='research/shape_search_20261121/projection/reproduce.py',points_saved=rec['points_saved'],burn=rec['burn'],observed=rec['window_length']*rec['windows'],cloud_sha256=rec['cloud_sha256'],mpfr_bits=[a['bits'] for a in audits],min_mpfr_third=precision_min,certified_cycles=row['certified_cycles'],certified_sinks=row['certified_sinks'],entropy_certificate='snapback/'+entropy['packet'] if entropy else None,entropy_lower_denominator=row['entropy_lower_denominator']))
    (args.source/'manifest.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8')
    with (args.source/'cases.csv').open('w',newline='',encoding='utf8') as f:
        out=csv.DictWriter(f,fieldnames=list(rows[0]));out.writeheader();out.writerows(rows)
    print('long',len(rows),'retained',len(data),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    main(p.parse_args())
