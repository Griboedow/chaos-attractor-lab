"""Freeze completed cluster, return, destination and initial-condition controls."""
import argparse,datetime,hashlib,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/motor-activity-continuation'
SCOPES=['cluster-spectrum','cluster-return-check','cluster-return-figure-review','probe-destinations',
        'cluster-strict','cluster-strict-completed','tuple15-excluded','tuple15-ensemble','tuple15-ensemble-audit']
SOURCES=['full_motor_activity_cluster_spectrum.py','check_motor_cluster_returns.py','plot_motor_cluster_returns.py',
    'check_motor_probe_destinations.py','check_motor_probe_confirmation.py','check_motor_probe_symmetry.py',
    'freeze_motor_cluster_checkpoint.py','resume_motor_cluster_control.py','classify_motor_tuple15_excluded.py',
    'check_motor_activity_initial_ensemble.py','report_motor_activity_initial_ensemble.py']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding='utf-8'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output.resolve()
    if out.exists():
        raise FileExistsError('Choose a new freeze folder')
    checks=[]
    def check(label,passed,**details):
        checks.append(dict(check=label,passed=bool(passed),**details))
    for name in ['cluster-spectrum','tuple15-excluded']:
        m=read(BASE/name/'manifest-final.json')
        check('child_manifest_passed',m['status']=='passed',scope=name)
        for path,entry in m['files'].items():
            check('child_manifest_hash',sha(BASE/name/path)==entry['sha256'],scope=name,path=path)
    m=read(BASE/'probe-destinations/manifest.json')
    for entry in m['files']:
        check('destination_packet_hash',sha(ROOT/entry['path'])==entry['sha256'],path=entry['path'])
    for name in ['checkpoint-02','checkpoint-03']:
        m=read(BASE/'cluster-strict'/name/'checkpoint.json')
        check('complete_block_checkpoint',m['status']=='passed' and all(m['checks'].values()),checkpoint=name)
        for path,expected in m['hashes'].items():
            check('checkpoint_hash',sha(BASE/'cluster-strict'/name/path)==expected,checkpoint=name,path=path)
    folder=BASE/'cluster-strict-completed';strict=read(folder/'analysis.json')
    check('strict_10000_complete',strict['status']=='complete' and strict['observed_tau']==10000. and strict['elapsed_tau']==10200.)
    check('strict_short_control',strict['short_step_control_maximum_state_and_tangent_error']<1e-10)
    check('strict_history_hash',sha(folder/'history.npz')==strict['history_sha256'])
    for name,digest in strict['source_hashes'].items():
        check('strict_source_hash',sha(folder/name)==digest,source=name)
    for name,digest in strict['input_hashes'].items():
        check('strict_input_hash',sha(folder/('input-'+name))==digest,source=name)
    inherited=np.load(folder/'input-history.npz');final=np.load(folder/'history.npz')
    original=np.load(BASE/'cluster-strict/lane-46/history.npz')
    for name in ['time','states']:
        check('strict_inherited_prefix',np.array_equal(inherited[name],final[name][:len(inherited[name])]),array=name)
        check('strict_reproduces_original_partial_trajectory',np.array_equal(original[name],final[name][:len(original[name])]),array=name)
    old=read(BASE/'cluster-strict/lane-46/analysis.json')
    check('original_partial_scope_preserved',old['status']=='wall budget incomplete' and old['elapsed_tau']>old['observed_tau']+200.)
    check('all20strict_windows',len(strict['windows'])==20 and all(w['duration']==500. for w in strict['windows']))
    for a,b in zip(old['windows'],strict['windows']):
        check('strict_existing_window',a==b,end_tau=a['end_tau'])
    endpoint=final['final_state'];check('strict_final_state',np.array_equal(endpoint,final['states'][-1]))
    q=np.array([w['quotient_leading_ftle'] for w in strict['windows']])
    c=np.array([w['contrast_ftle'] for w in strict['windows']])
    strict_summary=dict(quotient_leading=strict['quotient_leading_ftle'],contrast_leading=np.array(strict['contrast_ftle'])[:,0].tolist(),
        minimum500_tau_quotient=float(q.min()),maximum500_tau_contrast_leaders=c[:,:,0].max(axis=0).tolist(),
        short_control_max_error=strict['short_step_control_maximum_state_and_tangent_error'])
    returns=read(BASE/'cluster-return-check/analysis.json');sample=np.load(BASE/'cluster-return-check/sample.npz')
    check('return_check_complete',returns['status']=='complete' and returns['duration']==1000.)
    check('return_sample_hash',sha(BASE/'cluster-return-check/sample.npz')==returns['sample_sha256'])
    seed=np.load(BASE/'cluster-spectrum/lane-46/history.npz')
    check('return_exact_saved_state',np.array_equal(sample['initial_state'],seed['final_state']))
    check('return_section_order',len(sample['section_time'])==124 and np.all(np.diff(sample['section_time'])>0))
    check('return_section_residual',np.max(abs(sample['section_states'][:,16]))<1e-10)
    figure=read(BASE/'cluster-return-figure-review/analysis.json')
    check('return_figure_data',figure['sample_sha256']==returns['sample_sha256'])
    check('return_figure_hash',sha(BASE/'cluster-return-figure-review/returns.png')==figure['image_sha256'])
    ensemble=read(BASE/'tuple15-ensemble-audit/analysis.json')
    check('ensemble_passed',ensemble['status']=='passed' and all(c['passed'] for c in ensemble['checks']))
    for name,digest in ensemble['input_hashes'].items():
        check('ensemble_input_hash',sha(BASE/'tuple15-ensemble'/name)==digest,source=name)
    check('ensemble_reference_hash',sha(Path(ensemble['reference_path']))==ensemble['reference_sha256'])
    for name,digest in ensemble['source_hashes'].items():
        check('ensemble_source_hash',sha(BASE/'tuple15-ensemble-audit'/name)==digest,source=name)
    positive=[x for x in ensemble['outcomes'] if x['screen_positive']]
    other=[x for x in ensemble['outcomes'] if not x['screen_positive']]
    ensemble_summary=dict(declared=16,positive_screen=len(positive),subthreshold=len(other),
        subthreshold_maximum_saddle_distance=max(x['saddle_full_row_rms_distance'] for x in other),
        subthreshold_maximum_one_period_closure=max(x['one_period_full_state_closure_error'] for x in other),
        original_seed46_case=ensemble['outcomes'][1])
    if not all(c['passed'] for c in checks):
        raise ArithmeticError([c for c in checks if not c['passed']])
    out.mkdir(parents=True);(out/'sources').mkdir()
    for name in SOURCES:
        (out/'sources'/name).write_bytes((ROOT/'research'/name).read_bytes())
    (out/'audit-source.py').write_bytes(Path(__file__).read_bytes())
    for name in ['motor-activity-continuation.md','motor-lane46-spectrum.md','motor-tuple15-excluded.md']:
        (out/name).write_bytes((ROOT/'docs'/name).read_bytes())
    files=[]
    for folder in [BASE/name for name in SCOPES]+[out]:
        files.extend(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    entries=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
    result=dict(status='passed',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,
        strict_summary=strict_summary,ensemble_summary=ensemble_summary,scope_directories=SCOPES,
        file_count=len(entries),total_bytes=sum(e['bytes'] for e in entries),files=entries,
        local_figure_review='Full cluster spectrum, return intervals and excluded-lane orbit-distance PNGs inspected; browser rendering unavailable.',
        scope='Completed finite numerical follow-ups only. Separate interval proof, branch continuation and logarithmic source-start studies are excluded.')
    (out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('passed',len(checks),'checks',len(entries),'files',strict_summary,ensemble_summary['positive_screen'])


if __name__=='__main__':
    main()
