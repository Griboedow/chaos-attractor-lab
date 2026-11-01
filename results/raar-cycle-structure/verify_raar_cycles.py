"""Audit frozen outcomes, interval inequalities and source provenance."""
import argparse,json
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np
from verify_rk23_controller import read,sha


def interval(s):
    a,b=s.strip('[]').split(',');return Decimal(a.strip()),Decimal(b.strip())


def inclusion(record):
    box=[interval(v) for v in record['box']];K=[interval(v) for v in record['krawczyk']]
    assert all(b[0]<k[0]<=k[1]<b[1] for b,k in zip(box,K))
    c=interval(record.get('contraction',record.get('contraction_bound')));assert 0<=c[0]<=c[1]<1


def main(args):
    root=args.root;manifest=json.loads((root/'manifest.json').read_text());assert sha(args.base/'manifest.json')==manifest['base_manifest_sha256']
    paths={r['path'] for r in manifest['files']};actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts};assert actual==paths
    for row in manifest['files']:
        p=root/row['path'];assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],p
    assert manifest['count']==len(paths) and manifest['bytes']==sum(r['bytes'] for r in manifest['files'])
    for row in json.loads((root/'copy-provenance.json').read_text())['copies']:assert sha(root/row['path'])==row['sha256']
    for row in json.loads((root/'assets.json').read_text())['files']:assert sha(args.assets/row['path'])==sha(root/row['archive_path'])==row['sha256']
    baseline=read(args.base/'scout/screen.npz')['params'][1051];summary=json.loads((root/'continuation/report.json').read_text());assert len(summary['branches'])==4
    for row in summary['branches']:
        d=read(root/f'continuation/period-{row["period"]}-branch.npz');assert len(d['beta'])==401 and d['complete'].all() and (d['residual']<2e-10).all() and (d['least_period_separation']>1e-5).all()
        assert row['completed']==row['least_period_preserved']==401
    cycles=json.loads((root/'cycle-proof/report.json').read_text());assert cycles['status']=='passed' and len(cycles['cases'])==5
    assert sha(root/'cycle-proof/certify_raar_cycles.py')==cycles['source_sha256'] and sha(root/'cycle-proof/input.json')==cycles['input_sha256']
    for r in cycles['cases']:
        inclusion(r);assert np.array_equal(r['parameters'][1:],baseline[1:]) and r['minimal_period_certified']
        assert [d['divisor'] for d in r['proper_divisor_exclusions']]==[d for d in range(1,r['period']) if r['period']%d==0]
        for d in r['proper_divisor_exclusions']:assert any(a>0 or b<0 for a,b in map(interval,d['coordinate_difference']))
        assert all(interval(v)[0]>0 for v in r['axis_distances'])
        if r['label'].startswith('stable-'):assert r['stable_certified'] and all(interval(v)[0]>0 for v in r['schur_jury_margins'])
        else:assert r['unstable_certified'] and interval(r['characteristic_at_minus_one'])[0]>0 and r['parameters'][0]==baseline[0]
    events=json.loads((root/'multiplier-proof/report.json').read_text());assert events['status']=='passed' and [r['period'] for r in events['cases']]==[7,14,28,56]
    assert sha(root/'multiplier-proof/certify_raar_multiplier_events.py')==events['source_sha256'] and sha(root/'multiplier-proof/input.json')==events['input_sha256']
    for r in events['cases']:
        inclusion(r['uniform_branch_certificate']);a,b=interval(r['uniform_branch_certificate']['q_derivative_at_minus_one']);assert a>0 or b<0
        assert r['uniform_branch_certificate']['other_two_inside_at_event_certified'] and all(interval(v)[0]>0 for v in r['uniform_branch_certificate']['remaining_multipliers_at_event_schur_margins'])
        assert np.array_equal(r['base_parameters'][1:],baseline[1:]) and r['simple_multiplier_certified'] and r['minus_one_event_exists']
        assert r['endpoint_characteristic_signs']==[-1,1]
        for endpoint,sign,dyadic in zip(r['endpoints'],[-1,1],r['beta_bracket_exact_dyadics']):
            inclusion(endpoint);a,b=interval(endpoint['q_minus_one']);assert (a>0 if sign==1 else b<0) and endpoint['inside_uniform_branch_box']
            assert endpoint['beta_exact_dyadic']==dyadic
            if sign==-1:assert endpoint['stable_certified'] and all(interval(v)[0]>0 for v in endpoint['schur_jury_margins'])
        with localcontext() as ctx:
            ctx.prec=800;ends=[Decimal(r['signed_mantissa'])*Decimal(2)**r['binary_exponent'] for r in r['beta_bracket_exact_dyadics']];assert ends[0]<ends[1] and Decimal('1.9e-300')<ends[1]-ends[0]<Decimal('2.1e-300')
    fine=read(root/'fine-slice/scan.npz');assert np.array_equal(fine['params'][:,0],np.linspace(.9608,.9664,1025)) and np.array_equal(fine['params'][:,1:],np.tile(baseline[1:],(1025,1))) and (fine['counters'][:,0]==0).all()
    computed=np.zeros(1025,dtype=int)
    for i,z in enumerate(fine['tails']):
        for lag in range(1,513):
            if np.max(np.abs(z[lag:]-z[:-lag]))<1e-9:computed[i]=lag;break
    assert np.array_equal(computed,fine['period']);flags=(fine['rates'][:,:,0].min(1)>.001)&(computed==0);assert np.array_equal(flags,fine['flags'])
    report=json.loads((root/'fine-slice/report.json').read_text());assert int(flags.sum())==report['positive']==90 and np.sum(computed==0)==92
    assert report['period_counts']=={str(k):int(np.sum(computed==k)) for k in np.unique(computed)}
    assert np.max(np.abs(fine['rates'].sum(2)-fine['detlogs']))<1e-10
    structure=json.loads((root/'structure-audit/report.json').read_text());assert structure['status']=='passed' and structure['unstable_real_multipliers']==1 and structure['stable_real_multipliers']==2
    assert interval(structure['determinant'])[1]<0 and interval(structure['characteristic_at_plus_one'])[0]>0
    assert sha(root/'structure-audit/input.json')==sha(root/'cycle-proof/base-unstable-period-7.json')==structure['input_sha256']
    for r in structure['multiplier_brackets']:
        a,b=map(interval,r['endpoint_values']);assert (a[1]<0 and b[0]>0) or (a[0]>0 and b[1]<0)
    result=dict(status='passed',manifest_sha256=sha(root/'manifest.json'),files=len(paths),continuation_settings=1604,isolated_certificates=5,multiplier_event_certificates=4,fine_settings=1025,leading_positive=90,periodic=int(np.sum(computed>0)),unclassified=int(np.sum((computed==0)&~flags)),scope='Packet and certificate-output audit. Re-running the two independent interval drivers reproduces the mathematical inclusions; this audit alone is not an independent interval arithmetic implementation.')
    args.output.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('results/raar-cycle-structure'));p.add_argument('--base',type=Path,default=Path('results/raar-circles'));p.add_argument('--assets',type=Path,default=Path('docs/assets'));p.add_argument('--output',type=Path,required=True);main(p.parse_args())
