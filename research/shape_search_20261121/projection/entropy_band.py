"""Verify a frozen entropy itinerary uniformly on a closed interval of relaxations."""
import argparse
import json
from decimal import Decimal
from pathlib import Path
import mpmath as mp

import snapback
from analytic_bounds import iterate as analytic
from certify_cycles import iterate as jets


def check(packet,case,beta_digits,root_digits):
    record=dict(case,parameters=dict(case['parameters'],beta_radius=f'1e-{beta_digits}'))
    # Retain the original centers and inverse preconditioners. The root box now
    # encloses the continued periodic root for every beta in the parameter band.
    doc=dict(packet,fixed=dict(packet['fixed'],root_radius=f'1e-{root_digits}'))
    result={}
    for name,model in [('jet',jets),('analytic',analytic)]:
        snapback.enclosure_model=model
        proof=snapback.close_branches(doc,record)
        if proof['uniform_return_steps']!=packet['proof']['uniform_return_steps']:raise ArithmeticError('Unexpected return time')
        result[name]=dict(root_inclusion=proof['periodic_root_check']['max_krawczyk_fraction_display'],inverse_box_fraction=proof['largest_inverse_box_fraction_display'],local_inverse_lipschitz=proof['local_contraction_upper_display'])
    return dict(case=case['case'],beta_binary64=case['parameters']['beta'],beta_radius=f'1e-{beta_digits}',root_radius=f'1e-{root_digits}',return_steps=packet['proof']['uniform_return_steps'],checks=result)


def main(args):
    mp.mp.dps=mp.iv.dps=160
    c=json.loads(args.record.read_text());d=json.loads(args.packet.read_text())
    if args.verify:
        saved=json.loads(args.out.read_text());beta=-Decimal(saved['beta_radius']).adjusted();root=-Decimal(saved['root_radius']).adjusted()
        result=check(d,c,beta,root);print('verified',result,flush=True);return
    target=-Decimal(d['target_radius']).adjusted()
    # Search for usable widths. The final packet states only a verified interval.
    for width in (target+4,target+7,target+10,target+14,target+20,target+28,target+38,target+52,target+61):
        for root in (target+2,width-2,width-4,width-8,width-21):
            if root<=target:continue
            try:result=check(d,c,width,root)
            except (ArithmeticError,ValueError,ZeroDivisionError):continue
            result['packet']=args.packet.name
            args.out.parent.mkdir(parents=True,exist_ok=True)
            args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
            print('band',result,flush=True);return
    raise ArithmeticError('No tested relaxation band passed')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--record',type=Path,required=True);p.add_argument('--packet',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--verify',action='store_true')
    main(p.parse_args())
