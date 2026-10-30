"""Freeze the completed activity continuation without including later live runs."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/motor-activity-continuation'
SCOPES=['cpu','gpu-main','gpu-strict','gpu-audit','gpu-audit-recheck','lane49',
        'cluster-controls','full-space-perturbations','perturbation-audit','publication-01']
SOURCES=['continue_motor_activity_cpu.py','continue_motor_activity_gpu.py',
         'check_motor_activity_lane49.py','check_motor_activity_restart.py',
         'check_motor_activity_clusters.py','check_motor_activity_perturbations.py',
         'report_motor_activity_gpu.py','report_motor_activity_perturbations.py']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output=args.output.resolve()
    if args.output.exists():
        raise FileExistsError('Choose a new freeze directory')
    checks=[]
    def check(label,passed,**details):
        checks.append(dict(check=label,passed=bool(passed),**details))
    def snapshots(folder,record,key='source_snapshot_sha256'):
        for name,expected in record[key].items():
            check('archived_source_hash',sha(folder/name)==expected,file=(folder/name).relative_to(ROOT).as_posix())
    cpu=read(BASE/'cpu/manifest.json')
    check('cpu_manifest_passed',cpu['status']=='passed' and all(cpu['checks'].values()))
    for name,entry in cpu['files'].items():
        check('cpu_frozen_file_hash',sha(BASE/'cpu'/name)==entry['sha256'],file=name)
    for name in ['gpu-audit','gpu-audit-recheck','perturbation-audit']:
        data=read(BASE/name/'analysis.json')
        check('independent_analysis_passed',data['status']=='passed' and all(c['passed'] for c in data['checks']),folder=name)
        for filename,expected in data['input_sha256'].items():
            p=ROOT/filename if name.startswith('gpu') else BASE/'full-space-perturbations'/filename
            check('independent_analysis_input_hash',sha(p)==expected,file=p.relative_to(ROOT).as_posix())
        check('analysis_source_hash',sha(BASE/name/'source.py')==data['source_sha256'],folder=name)
    for setting in ['gpu-main','gpu-strict']:
        folder=BASE/setting;record=read(folder/'summary.json')
        check('cohort_complete',record['status']=='complete' and len(record['blocks'])==40,setting=setting)
        check('cohort_input_hash',sha(ROOT/record['source_path'])==record['source_sha256'],setting=setting)
        check('cohort_prior_block_hash',sha(ROOT/record['source_final_block'])==record['source_final_block_sha256'],setting=setting)
        snapshots(folder,record)
    for lane in [36,42,46]:
        folder=BASE/'cluster-controls'/f'lane-{lane}';record=read(folder/'analysis.json')
        check('conditional_complete',record['status']=='complete' and record['observed_tau']==2000.,lane=lane)
        check('conditional_input_hash',sha(Path(record['source_path']))==record['source_sha256'],lane=lane)
        check('conditional_embedding',record['transverse_full_jacobian_error']<1e-10,lane=lane)
        snapshots(folder,record)
    sample=read(BASE/'lane49/sample.json')
    check('lane49_sample_found_return',sample['status']=='periodic_return_found')
    snapshots(BASE/'lane49',sample,'source_sha256_by_name')
    cycle=read(BASE/'lane49/periodic/analysis.json')
    check('lane49_periodic_completed',cycle['status']=='complete')
    check('lane49_cycle_input',sha(ROOT/cycle['source'])==cycle['source_sha256'])
    check('lane49_cycle_source',sha(BASE/'lane49/periodic/source.py')==cycle['code_sha256'])
    restart=read(BASE/'lane49/restart-control/analysis.json')
    check('lane49_restart_complete',restart['status']=='complete')
    snapshots(BASE/'lane49/restart-control',restart)
    transverse=read(BASE/'lane49/transverse-control-provenance.json')
    check('lane49_transverse_output_hash',sha(BASE/'lane49/transverse-control.json')==transverse['output_sha256'])
    for name,expected in transverse['sources'].items():
        check('lane49_transverse_source',sha(ROOT/name)==expected,file=name)
    protocol=read(BASE/'full-space-perturbations/protocol.json')
    snapshots(BASE/'full-space-perturbations',protocol)
    failed=[c for c in checks if not c['passed']]
    if failed:
        raise ArithmeticError(failed)
    args.output.mkdir(parents=True)
    (args.output/'sources').mkdir()
    for name in SOURCES:
        (args.output/'sources'/name).write_bytes((ROOT/'research'/name).read_bytes())
    (args.output/'audit-source.py').write_bytes(Path(__file__).read_bytes())
    (args.output/'report.md').write_bytes((ROOT/'docs/motor-activity-continuation.md').read_bytes())
    files=[]
    for folder in [BASE/name for name in SCOPES]+[args.output]:
        files.extend(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    entries=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed',checks=checks,
        scope_directories=SCOPES,excludes_later_live_runs=True,file_count=len(entries),total_bytes=sum(e['bytes'] for e in entries),files=entries,
        review=dict(local_figures=['cpu/analysis/cpu-spectrum.png','gpu-audit-recheck/cohort-continuation.png','perturbation-audit/perturbation-outcomes.png'],
            figures_inspected=True,browser_rendering='Unavailable in this session; local Markdown paths, math fences and PNGs checked.',
            independent_source_review='Weighted restriction, eliminated displacement, contrast embeddings, column QR convention and physical perturbation norm reviewed read-only; no defect found.',
            report_snapshot='The current report remains editable; this directory preserves its publication text.'),
        limitations=['Finite-time numerical evidence; no new strange-attractor or priority theorem.',
            'Both source snapshots actually executed by each run and current reproduction interfaces are preserved.',
            'This inventory excludes itself and does not certify future extensions.'])
    (args.output/'manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result['status'],len(checks),'checks',len(entries),'files',result['total_bytes'],'bytes')


if __name__=='__main__':
    main()
