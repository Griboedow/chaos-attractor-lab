"""Check saved postcritical evidence and emit an exact-byte manifest.

This does not rerun integrations or upgrade numerical evidence to a proof.
"""
import argparse
import datetime
import hashlib
import json
import re
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/motor-postcritical'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    if output.exists():raise FileExistsError(output)
    files=sorted(p for p in BASE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p!=output)
    failures=[];checks=[]
    sources=list((ROOT/'research').glob('*.py'))+list((ROOT/'research').glob('*.cl'))+[p for p in files if p.suffix=='.py']
    source_hashes={}
    for p in sources:
        data=p.read_bytes()
        for value in [hashlib.sha256(data).hexdigest(),hashlib.sha256(data.replace(b'\r\n',b'\n')).hexdigest()]:
            source_hashes.setdefault(value,[]).append(p.relative_to(ROOT).as_posix())
    source_references=[]
    for p in files:
        if p.suffix=='.npz':
            with np.load(p,allow_pickle=False) as data:
                valid=all(np.isfinite(data[name]).all() for name in data.files if np.issubdtype(data[name].dtype,np.number))
                if 'status' in data.files:valid=valid and bool(np.all(data['status']==0))
            checks.append(dict(file=p.relative_to(ROOT).as_posix(),kind='finite arrays and completed GPU lanes',passed=bool(valid)))
            if not valid:failures.append(p.relative_to(ROOT).as_posix()+': nonfinite/incomplete arrays')
        if p.suffix!='.json':continue
        value=json.loads(p.read_text())
        for key in ['script_sha256','code_sha256','renderer_sha256']:
            expected=value.get(key)
            if isinstance(expected,str):
                found=source_hashes.get(expected,[])
                source_references.append(dict(file=p.relative_to(ROOT).as_posix(),key=key,sha256=expected,matches=found))
                if not found:failures.append(p.relative_to(ROOT).as_posix()+': missing executed source '+expected)
            elif isinstance(expected,dict):
                for name,h in expected.items():
                    found=source_hashes.get(h,[])
                    source_references.append(dict(file=p.relative_to(ROOT).as_posix(),key=key+':'+name,sha256=h,matches=found))
                    if not found:failures.append(p.relative_to(ROOT).as_posix()+': missing executed source '+name)
        if isinstance(value.get('source_sha256'),str):
            name=value.get('source_file',value.get('source'))
            if isinstance(name,str):
                source=ROOT/name
                if source.is_dir():
                    candidates=[source/'cycle.npz',source/'analysis.json']
                    source=next((candidate for candidate in candidates if candidate.is_file() and digest(candidate)==value['source_sha256']),candidates[0])
                passed=source.is_file() and digest(source)==value['source_sha256']
                checks.append(dict(file=p.relative_to(ROOT).as_posix(),kind='upstream exact source bytes',passed=passed,source=source.relative_to(ROOT).as_posix()))
                if not passed:failures.append(p.relative_to(ROOT).as_posix()+': upstream source mismatch')
    for manifest in [BASE/'figure-manifest-v2.json',BASE/'figure-preview-initial/archive-render-manifest.json']:
        values=json.loads(manifest.read_text())
        for name,h in values['images'].items():
            passed=digest(ROOT/name)==h
            checks.append(dict(file=name,kind='rendered image hash',passed=passed))
            if not passed:failures.append(name+': figure hash mismatch')
    report=ROOT/'docs/motor-postcritical.md'
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
        if '://' in target or target.startswith('#'):continue
        path=(report.parent/target.split('#')[0]).resolve()
        if path!=output and not path.exists():failures.append('Missing report link: '+target)
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status='passed' if not failures else 'failed',failures=failures,checks=checks,
        executed_source_references=source_references,
        files=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in files],
        code_files={p.relative_to(ROOT).as_posix():digest(p) for p in sorted((ROOT/'research').glob('*motor_postcritical*.py'))},
        report=dict(path=report.relative_to(ROOT).as_posix(),sha256=digest(report)),
        historical_notes=['figure-manifest.json identifies the superseded preview; exact original PNGs are preserved under figure-preview-initial.',
                          'initial-survey/samples/source-snapshot-note.json identifies the executed driver despite a CLI repair during sampling.',
                          'final-six-cluster-rho-3p268/interpretation.json excludes a stale-group projected shooting result from the connection evidence.'],
        scientific_scope='Provenance, data-completion and link checks only. No additional scientific integration and no interval theorem.')
    output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=record['status'],files=len(files),checks=len(checks),source_references=len(source_references),failures=failures),indent=2))
    if failures:raise SystemExit(2)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=BASE/'manifest.json')
    a=parser.parse_args();main(a.output.resolve())
