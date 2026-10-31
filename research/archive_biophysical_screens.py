"""Archive completed source checks and nonretained numerical scouts."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from crossfeeding import field_jac as crossfeeding_field
from rho_adsi_kinetics import field_jac as rho_field


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);inventory=[]
    names=['rho-adsi-pilot','rho-adsi-focus','rho-adsi-kinetics','rho-adsi-checks','rho-adsi-checks02',
        'rho-adsi-return-audit','rho-adsi-focus-returns','rho-adsi-focus-returns02','rho-adsi-kinetics-returns',
        'crossfeeding-source','crossfeeding-pilot','crossfeeding-recheck','crossfeeding-returns',
        'crossfeeding-recheck-returns','motor-clutch-source']
    scout_names={'rho-adsi-pilot','rho-adsi-focus','rho-adsi-kinetics','crossfeeding-pilot','crossfeeding-recheck'}
    for name in names:
        source=args.external/name;target=out/name;target.mkdir()
        for path in sorted(source.iterdir()):
            if path.suffix in ['.json','.py','.cl']:shutil.copyfile(path,target/path.name)
            elif path.suffix=='.npz':
                keep=name in scout_names or name=='crossfeeding-source'
                if keep:shutil.copyfile(path,target/path.name)
                inventory.append(dict(run=name,name=path.name,bytes=path.stat().st_size,sha256=digest(path),included=keep))
    author=args.external/'ross-author-code/manifest.json';shutil.copyfile(author,out/'crossfeeding-author-source-manifest.json')
    summaries=[]
    expected={'rho-adsi-pilot':(2050,1715),'rho-adsi-focus':(12860,12756),'rho-adsi-kinetics':(8016,7883),
        'crossfeeding-pilot':(8194,5552),'crossfeeding-recheck':(8194,7318)}
    for name,(total,completed) in expected.items():
        path=out/name/'complete-or-budgeted.npz';d=np.load(path)
        assert len(d['status'])==total and (d['status']==0).sum()==completed
        primary=(d['status']==0)&(d['windows'].min(axis=1)>.005);assert not primary.any()
        summaries.append(dict(run=name,total=total,completed=completed,primary_retained=int(primary.sum()),sha256=digest(path),
            statuses={str(int(s)):int((d['status']==s).sum()) for s in np.unique(d['status'])}))
    alignments=[]
    for name,dimension,field in [('rho-adsi-kinetics',3,rho_field),('crossfeeding-recheck',5,crossfeeding_field)]:
        d=np.load(out/name/'complete-or-budgeted.npz');ids=np.flatnonzero((d['status']==0)&(d['windows'].min(axis=1)>.001));rows=[]
        for i in ids:
            f=field(d['states'][i,:dimension],d['parameters'][i])[0];v=d['states'][i,dimension:2*dimension]
            cosine=float(abs(f@v)/(np.linalg.norm(f)*np.linalg.norm(v)))
            rows.append(dict(id=int(i),one_minus_absolute_cosine=1-cosine))
        alignments.append(dict(run=name,weak_count=len(ids),maximum_one_minus_absolute_cosine=max(r['one_minus_absolute_cosine'] for r in rows),runs=rows))
    assert alignments[0]['weak_count']==97 and alignments[1]['weak_count']==44
    focus=[]
    for name in ['rho-adsi-focus-returns','rho-adsi-focus-returns02']:
        focus.extend(json.loads(p.read_text()) for p in sorted((out/name).glob('case-*.json')))
    assert len(focus)==32 and len({r['id'] for r in focus})==32
    assert sum(r['first_numerical_period_lag']==1 for r in focus)==30
    cross=[]
    for name in ['crossfeeding-returns','crossfeeding-recheck-returns']:
        runs=json.loads((out/name/'analysis.json').read_text())['runs'];assert len(runs)==4
        assert all(r['first_numerical_period_lag']==1 for r in runs);cross.extend(runs)
    write(out/'local-inventory.json',inventory)
    write(out/'audit.json',dict(scouts=summaries,endpoint_tangent_flow_alignments=alignments,
        rho_focus_return_checks=32,rho_focus_lag_one=30,rho_focus_insufficient_events=[r['id'] for r in focus if r['first_numerical_period_lag'] is None],
        crossfeeding_return_checks=len(cross),crossfeeding_lag_one=len(cross),
        scope='No retained chaos candidate. Finite-time screens and numerical recurrences are not global exclusions. Rechecks overlap original designs; do not add their counts as independent parameter coverage.'))
    entries=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write(out/'manifest.json',dict(count=len(entries),bytes=sum(e['bytes'] for e in entries),files=entries))
    print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),checks='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--external',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
