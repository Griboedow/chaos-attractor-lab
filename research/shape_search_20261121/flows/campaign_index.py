"""Count actual distinct screened coefficient vectors and separate controls and geometries."""
import argparse
import csv
import functools
import gzip
import hashlib
import json
import os
from collections import Counter,defaultdict
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from field import make_case,make_ring,make_rotor,make_planes,make_crossed

GENERATORS=dict(cartesian=make_case,ring=make_ring,rotor=make_rotor,planes=make_planes,crossed=make_crossed)


def fingerprint(p):
    return hashlib.sha256(np.asarray(p,dtype='<f8').ravel().tobytes()).hexdigest()


@functools.lru_cache(maxsize=90000)
def generated(family,seed,case):
    return fingerprint(GENERATORS[family](seed,case)[0])


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--flip',required=True,type=Path)
    a=ap.parse_args();b=a.cache;r=a.out
    stages=defaultdict(lambda:defaultdict(set));counts=defaultdict(Counter)
    def add(stage,family,dig):
        stages[family][stage].add(dig);counts[family][stage]+=1
    archive_agreements=Counter()
    for file in sorted(b.glob('*/screen.csv')):
        folder=file.parent;spec=read(folder/'run.json');family=spec.get('family','cartesian');seed=spec['seed']
        with file.open() as f:rows=list(csv.DictReader(f))
        if not rows: continue
        if 'status' in rows[0]:
            tag='period_screen'
        elif spec.get('ids'):
            tag='fixed_step_replay'
        else:
            tag='broad_chaos_screen'
        overrides=read(Path(spec['inputs'])) if spec.get('inputs') else {}
        for row in rows:
            case=int(row['case'])
            dig=fingerprint(overrides[str(case)]['parameters']) if overrides else generated(family,seed,case)
            add(tag,family,dig)
        # Independently check available saved parameter arrays from each primary screening batch.
        if tag=='broad_chaos_screen':
            for f in folder.glob('case-*.npz'):
                case=int(f.stem[5:])
                with np.load(f) as z: assert fingerprint(z['parameters'])==generated(family,seed,case)
                archive_agreements[family]+=1
    phase=b/'band-scan-a'
    base_folder=Path(read(phase/'settings.json')['source'])
    p_records={}
    with (phase/'scan.csv').open() as inp:
        for row in csv.DictReader(inp):
            base=int(row['base']);i=int(row['parameter']);t=float(row['fraction'])
            if base not in p_records:p_records[base]=np.asarray(read(base_folder/f'case-{base}.json')['parameters'])
            p=p_records[base].copy();p[i]+=t*max(.45,abs(p[i]))
            add('single_coefficient_screen','ring',fingerprint(p))
    m=read(r/'manifest.json');branch=read(a.flip/'report.json')
    if not branch['trials'] or branch['trials'][-1]['fraction']!=1:
        raise ValueError('The sampled period-doubling branch has not reached the chaotic endpoint')
    for row in branch['trials']:
        p=p_records[28756].copy();p[8]=row['e']
        add('period_doubling_grid','ring',fingerprint(p))
    for row in branch['bisection_steps']:
        p=p_records[28756].copy();p[8]=row['e'];add('diagnostic_bisection','ring',fingerprint(p))
    prim=['broad_chaos_screen','period_screen','single_coefficient_screen','period_doubling_grid']
    total=0;table={}
    for family in GENERATORS:
        group=stages[family];unique=set().union(*(group[k] for k in prim));total+=len(unique)
        # The hash in each row includes every distinct little-endian binary64 coefficient sequence.
        table[family]=dict(rows=counts[family],unique_per_stage={k:len(v) for k,v in group.items()},
            unique_screened_equations=len(unique),extra_from_period_screen=len(group['period_screen']-group['broad_chaos_screen']),
            additional_control_only=len(set().union(*group.values())-unique),
            equation_set_sha256=hashlib.sha256(''.join(sorted(unique)).encode()).hexdigest(),
            raw_coefficient_archives_crosschecked=archive_agreements[family])
    stable=read(r/'periodic-accounting.json');controls=read(r/'controls.json')
    recov=read(r/'recovery.json');orbit=read(r/'periodic-orbits.json');links=read(r/'link-checks.json')
    outcomes=read(r/'geometry.json')
    count=lambda z:dict(Counter(row['family'] for row in z))
    per_family={f:dict(chaotic_selected=count(m['cases']).get(f,0),stable_knots=stable['families'].get(f,0),
        linked_sets=count(m['linked_cycles']).get(f,0),unstable_cores=count(m['knotted_periodic_skeletons']).get(f,0),
        tested_chaos_starts=sum(c['starts_tested'] for c in m['cases'] if c['family']==f),
        positive_chaos_starts=sum(c['starts_passed'] for c in m['cases'] if c['family']==f)) for f in GENERATORS}
    broad={key:dict(rows=sum(c['starts_tested'] for c in m['cases']),passing=sum(c['starts_passed'] for c in m['cases']))
           for key in ['standard_fresh_start_growth']}
    broad['local_recapture']=dict(cases=len(recov),trials=sum(sum(t['starts'] for t in groups) for groups in recov.values()),
                                  passed=sum(sum(t['passed'] for t in groups) for groups in recov.values()))
    broad['chaotic_adaptive_spectra']=len(controls)
    broad['sliced_geometry_cases']=len(outcomes)
    broad['unstable_short_returns_in_chaotic_clouds']=dict(cases=sum(bool(o['orbits']) for o in orbit.values()),
                                                           periodic_orbits=sum(len(o['orbits']) for o in orbit.values()))
    broad['stable_cycle_perturbations']=6*stable['distinct_accepted_parameter_cases']
    broad['linked_basin_perturbations']=dict(trials=sum(len(l['basins']['trials']) for l in links.values()),
             captured=sum(t['captured'] for l in links.values() for t in l['basins']['trials']))
    broad['selected_geometry_index']=dict(chaotic_clouds=len(m['cases']),stable_orbits=stable['distinct_accepted_parameter_cases'],
             composite_linked_sets=len(m['linked_cycles']),unstable_knotted_cores=len(m['knotted_periodic_skeletons']))
    if (r/'nested-skeletons.json').exists():
        broad['supplemental_same_parameter_cores']=len(read(r/'nested-skeletons.json')['cores'])
    if (r/'doubling-pairs.json').exists():
        broad['supplemental_parent_daughter_pairs']=len(read(r/'doubling-pairs.json'))
    references={}
    for file in ('archive-comparison.json.gz','ode-reference-comparison.json.gz'):
        if (r/file).exists():
            z=json.load(gzip.open(r/file,'rt'))
            references[file]=dict(total=len(z),per_collection=dict(Counter(k.split('/')[0] for k in z)),
                 reference_counts=sorted({item['reference_count'] for x in z.values() for item in x['compared'].values()}))
    summary=dict(unique_screened_equations=total,by_family=table,selected_by_family=per_family,
                 selected_control_counts=broad,archive_cohorts=references,
                 interpretation=['Independent screened equations are exact binary64 parameter vectors per family in a primary survey, a periodic screen, a one-coefficient sweep or the 36-point bifurcation grid.',
                    'Initial conditions, time steps, solver calls and component symmetries do not create new equations. Root-finding and diagnostic bisection probes are excluded from the primary total.',
                    'Equation counts do not quotient by symmetries, change of units or dynamical conjugacy.',
                    'One linked set contains several coexisting cycles. Each unstable core uses the same coefficients as its chaotic cloud. A cataloged geometry is a reproducible point sample, not a distinct mechanism.'])
    (r/'index-summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(total,{f:data['unique_screened_equations'] for f,data in table.items()},'references',references)


if __name__=='__main__':
    main()
