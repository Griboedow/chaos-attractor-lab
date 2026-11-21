"""Assemble portable case parameters and compact verification records."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--selections',type=Path,required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plan = read(args.selections)
    cases = []
    controls = {}
    periodic = {}
    recovery = {}
    for bundle in plan:
        for case, label in bundle['cases']:
            source = args.cache/bundle['dense']/f'case-{case}.npz'
            with np.load(source) as f:
                p, initial, rates = f['parameters'],f['initial'],f['rates']
                shape = list(f['points'].shape)
            valid = read(args.cache/bundle['adaptive']/f'case-{case}.json')
            checks = []
            for name in bundle['controls']:
                archive = args.cache/name/f'case-{case}.npz'
                with np.load(archive) as f:
                    checks.append(dict(run=name,initial=f['initial'].tolist(),rates=f['rates'].tolist()))
            controls[str(case)] = dict(family=bundle['family'],steps=checks,adaptive=valid)
            for name in bundle.get('local_controls', []):
                path = args.cache/name/f'case-{case}.npz'
                if path.exists():
                    with np.load(path) as f:
                        controls[str(case)].setdefault('nearby_initial_states',[]).append(
                            dict(run=name,initial=f['initial'].tolist(),rates=f['rates'].tolist()))
            orbits = args.cache/bundle.get('orbits','no-orbit-run')/f'case-{case}.json'
            if orbits.exists():
                periodic[str(case)] = read(orbits)
            returned = args.cache/bundle['dense'].replace('dense-','return-')/f'case-{case}.json'
            if returned.exists():
                groups = read(returned)['groups']
                recovery[str(case)] = [
                    {**{k:v for k,v in group.items() if k!='trials'},
                     'minimum_rate_each_trial':[min(r['rates']) for r in group['trials']],
                     'distance_p95_each_trial':[r['distance_p95_rms_units'] for r in group['trials']]}
                    for group in groups]
            cases.append(dict(case=case,family=bundle['family'],label=label,path=str(source.resolve()),
                              cloud_relative=source.resolve().relative_to(args.cache.resolve().parent).as_posix(),
                              point_key='points',shape=shape,sha256=sha(source),parameters=p.tolist(),
                              initial=initial.tolist(),dense_run=bundle['dense'],rates=rates.tolist(),
                              adaptive_lambda1=valid['mean_spectrum'][0],
                              starts_passed=sum(min(r['rates'])>.008 for r in checks),starts_tested=len(checks)))
    names = sorted({b[k] for b in plan for k in ('dense','adaptive')}|
                   {n for b in plan for n in b['controls']+b.get('local_controls',[])})
    survey = {}
    for path in sorted(args.cache.glob('*/screen.csv')):
        spec = read(path.parent/'run.json')
        if 'count' not in spec or 'ids' not in spec:
            continue
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        expected = len(spec['ids']) if spec['ids'] else spec['count']
        if len(rows) == expected:
            survey[path.parent.name] = dict(requested=expected,evaluated=len(rows),positive_and_spatial=sum(int(r['keep']) for r in rows),
                                     nonfinite=sum(not np.isfinite(float(r['window_0'])) for r in rows),
                                     csv_sha256=sha(path),run_settings={k:v for k,v in spec.items() if k!='out'})
    protocols = {n:{k:v for k,v in read(args.cache/n/'run.json').items() if k!='out'} for n in names}
    selected_ids = {case['case'] for case in cases}
    with (args.out/'long-screen.csv').open('w',encoding='utf-8',newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['family','case','run','window_length','count','minimum_lambda1',
                         'mean_lambda1','spatial_growth_filter','selected'])
        for run, record in survey.items():
            if 'long' not in run:
                continue
            settings = record['run_settings']
            with (args.cache/run/'screen.csv').open(encoding='utf-8') as source:
                for row in csv.DictReader(source):
                    values = [float(row[f'window_{i}']) for i in range(settings['blocks'])]
                    writer.writerow([settings.get('family','cartesian'),row['case'],run,settings['window'],
                                     len(values),min(values),round(float(np.mean(values)),8),row['keep'],
                                     int(int(row['case']) in selected_ids)])
    stable = []
    knot_controls = {}
    knot_starts = []
    knot_plan = args.out/'knot-selections.json'
    if knot_plan.exists():
        witnesses = args.out/'knot-witnesses'
        witnesses.mkdir(exist_ok=True)
        for group in read(knot_plan):
            hits = {h['case']:h for h in (json.loads(line) for line in
                    (args.cache/group['screen']/'hits.jsonl').read_text().splitlines())}
            for case,label in group['cases']:
                detail = read(args.cache/group['validated']/f'case-{case}.json')
                source = args.cache/group['validated']/f'case-{case}.npz'
                witness = args.cache/group['witness']/f'case-{case}.json'
                registered = witnesses/witness.name
                registered.write_bytes(witness.read_bytes())
                knot_starts.append(hits[case])
                knot_controls[str(case)] = detail
                stable.append(dict(case=case,family=detail['family'],label=label,
                    evidence='Numerically attracting periodic orbit with exact integer-polygon witness',
                    path=str(source.resolve()),cloud_relative=source.resolve().relative_to(args.cache.resolve().parent).as_posix(),
                    point_key='points',shape=[32768,3],sha256=sha(source),parameters=detail['parameters'],
                    initial=detail['refined_start'],period=detail['period'],determinant=detail['determinant'],
                    alexander_coefficients=detail['alexander_coefficients'],
                    max_transverse_multiplier_abs=detail['max_transverse_multiplier_abs'],
                    perturbations_returned=sum(v['p95_distance']<.02 for v in detail['recovery']),
                    perturbations_tested=len(detail['recovery']),
                    polygon_witness=registered.relative_to(args.out).as_posix()))
    linked=[]
    link_plan=args.out/'link-selections.json'
    if link_plan.exists():
        target=args.out/'link-witnesses';target.mkdir(exist_ok=True)
        link_checks={}
        for group in read(link_plan):
            for case,label in group['cases']:
                file=f'case-{case}.json'
                source=args.cache/group['source']/f'case-{case}.npz'
                witness=args.cache/group['source']/file
                base=read(witness);fine=read(args.cache/group['fine']/file)
                assert all(v['linking']==base['views'][0]['linking'] for v in base['views']+fine['views'])
                (target/file).write_bytes(witness.read_bytes())
                checks=read(args.cache/group['basins']/file)
                sections=read(args.cache/group['sections']/file)
                link_checks[str(case)]=dict(section=sections,basins=checks,fine_witness_sha256=sha(args.cache/group['fine']/file),
                        fine_vertices_per_component=len(fine['integer_polygons'][0]),fine_projected_views=fine['views'])
                linked.append(dict(case=case,label=label,family=base['family'],period=base['period'],
                      component_count=base['component_count'],shape=[32768*base['component_count'],3],
                      point_key='points',component_key='component',offsets_key='offsets',path=str(source.resolve()),
                      cloud_relative=source.resolve().relative_to(args.cache.resolve().parent).as_posix(),
                      sha256=sha(source),row_vector_transforms=base['transforms'],linking_matrix=base['views'][0]['linking'],
                      polygon_witness=(target/file).relative_to(args.out).as_posix(),representative_case=case))
        write(args.out/'link-checks.json',link_checks)
    skeletons=[]
    band_plan=args.out/'band-selections.json'
    if band_plan.exists():
        target=args.out/'band-witnesses';target.mkdir(exist_ok=True)
        bands={}
        band_audit={}
        for group in read(band_plan):
            for case in group['cases']:
                file=f'case-{case}.json'
                source=args.cache/group['source']/f'case-{case}.npz'
                b=read(args.cache/group['source']/file)
                assert b['completed'] and b['same_alexander_polynomial']
                cert=b.pop('projection_witness')
                cert['views']=[dict(integer_matrix=v['matrix'],crossings=v['crossings'],coefficients=v['coefficients'])
                               for v in b['final_exact_views']]
                cert['case']=case;cert['quantization_per_rms']=10**9
                (target/file).write_text(json.dumps(cert,separators=(',',':'))+'\n',encoding='utf-8')
                sections=args.cache/group['sections']/file
                if sections.exists(): b['finite_difference_section']=read(sections)
                if group.get('audited_case')==case:
                    tracks=[]
                    for folder in group.get('long_tracks',[]):
                        with np.load(args.cache/folder/f'case-{case}.npz') as f:
                            spec=read(args.cache/folder/'run.json')
                            tracks.append(dict(initial=f['initial'].tolist(),rates=f['rates'].tolist(),
                               dt=spec['dt'],burn=spec['burn'],window=spec['window'],blocks=spec['blocks']))
                    band_audit[str(case)]=dict(geometry=read(args.cache/group['audit']/'geometry.json'),
                               tight_section=read(args.cache/group['tight_section']/file),long_runs=tracks)
                bands[str(case)]=b
                skeletons.append(dict(case=case,base_case=b['base'],family='ring',shape=[32768,3],point_key='points',
                       path=str(source.resolve()),cloud_relative=source.resolve().relative_to(args.cache.resolve().parent).as_posix(),
                       sha256=sha(source),refined_start=b['refined_start'],period=b['final']['period'],
                       dominant_multiplier=b['final']['dominant'],alexander_coefficients=b['final_exact_views'][0]['coefficients'],
                       polygon_witness=(target/file).relative_to(args.out).as_posix(),
                       orbit_to_cloud_distance_p95=b['orbit_to_cloud_distance_p95'],cloud_to_orbit_distance_p95=b['cloud_to_orbit_distance_p95']))
        write(args.out/'band-continuation.json',bands)
        if band_audit:write(args.out/'band-audit.json',band_audit)
    write(args.out/'manifest.json',dict(point_columns=['x','y','z'],point_dtype='float64',cases=cases,
              stable_cycles=stable,linked_cycles=linked,knotted_periodic_skeletons=skeletons))
    if stable:
        write(args.out/'knot-controls.json',knot_controls)
        (args.out/'knot-starts.jsonl').write_text(''.join(json.dumps(h)+'\n' for h in knot_starts),encoding='utf-8')
    write(args.out/'controls.json',controls)
    write(args.out/'accounting.json',survey)
    write(args.out/'run-settings.json',protocols)
    write(args.out/'periodic-orbits.json',periodic)
    write(args.out/'recovery.json',recovery)


if __name__ == '__main__':
    main()
