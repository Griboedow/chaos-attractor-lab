"""Recheck stored cycles and sinks against both interval model implementations."""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
from analytic_bounds import iterate as analytic
from certify_cycles import iterate as jets,verify


def main(args):
    mp.mp.dps=mp.iv.dps=110
    rows=[]
    for group in ('cycles','sinks'):
        for file in sorted((args.source/group).glob('cycles-*.json')):
            doc=json.loads(file.read_text());case=doc['case']
            rec=json.loads((args.source/f'case-{case}.json').read_text())
            for i,c in enumerate(doc['cycles']):
                results={}
                for name,model in [('jet',jets),('analytic',analytic)]:
                    result=verify(c,rec,c['proof']['radius'],model)
                    if result['routh_signs']!=c['proof']['routh_signs']:raise ArithmeticError('Inconsistent spectral signs')
                    results[name]=result['contraction_upper_display']
                rows.append(dict(case=case,group=group,index=i,period=c['period'],unstable=c['proof']['unstable_dimension'],contraction_display=results))
    args.out.write_text(json.dumps(dict(mpmath=mp.__version__,precision_decimal=110,checked=len(rows),checks=rows),indent=2)+'\n',encoding='utf8')
    print('verified',len(rows),'cycles and sinks with both interval maps',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    main(p.parse_args())
