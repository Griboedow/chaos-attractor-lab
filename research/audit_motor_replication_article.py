"""Audit the article evidence chain and package a reviewed manuscript.

Read-only by default. --write-records creates review/manifest files once,
after visual inspection is explicitly recorded with --visual-review-complete.
This checks provenance and reported inequalities; it does not replace the
validated integrations or the analytic argument.
"""
import argparse
import datetime
from fractions import Fraction
import hashlib
import importlib.metadata
import itertools
import json
from pathlib import Path
import platform
import re

import numpy as np
from PIL import Image
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT/'papers/replication-instability'
DATA=ROOT/'results/motor-replication-article'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def bounds(record):
    def dyadic(pair):
        n,e=int(pair[0]),int(pair[1])
        return Fraction(n)*Fraction(2)**e
    m,r=dyadic(record['midpoint']),dyadic(record['radius'])
    return m-r,m+r


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-records',action='store_true')
    parser.add_argument('--visual-review-complete',action='store_true')
    args=parser.parse_args();checks=[]
    def check(name,passed):
        checks.append(dict(name=name,passed=bool(passed)))
    for packet in ['motor-five-cycle-proof','motor-five-replication-proof']:
        manifest=read(ROOT/'results'/packet/'manifest.json')
        for item in manifest['files']:
            path=ROOT/item['path']
            check('frozen:'+item['path'],path.exists() and sha(path)==item['sha256'])
    cycle_path=ROOT/'results/motor-five-cycle-proof/recheck-01.json'
    cycle=read(cycle_path)
    check('cycle exists and full-space spectral bounds passed',
          cycle['status']=='attracting_five_cluster_cycle_certified' and cycle['spectrum']['all_nonflow_moduli_below_point13'])
    check('reported Krawczyk contraction',bounds(cycle['krawczyk_contraction'])[1]<Fraction('1.960e-19'))
    check('reported inclusion ratio',max(bounds(x)[1] for x in cycle['krawczyk_ratios'])<Fraction('2.574e-12'))
    lo,hi=bounds(cycle['positive_period'])
    check('printed positive-period enclosure',Fraction('25.0016843758970971')<lo<hi<Fraction('25.0016843758970972'))
    for filename,source in [('exact-index.json',ROOT/'results/motor-five-replication-proof/singleton-blocks.json'),
                            ('exact-index-recheck.json',DATA/'reintegration-640.json')]:
        record=read(DATA/filename)
        check(filename+':input chain',record['input_sha256']==sha(source) and record['cycle_sha256']==sha(cycle_path))
        check(filename+':executed checker',record['source_sha256']==sha(ROOT/'research/certify_motor_replication_index.py')==sha((DATA/filename).with_suffix('.source.py')))
        check(filename+':status',record['status']=='exact_replication_index_certified')
        for block,value in record['blocks'].items():
            check(filename+':root '+block,bounds(value['polynomial_left'])[1]<0<bounds(value['polynomial_right'])[0]
                  and Fraction('-23.083469180')<bounds(value['root_left'])[0]
                  and bounds(value['root_right'])[1]<Fraction('-23.083469172'))
            check(filename+':quadratic '+block,all(bounds(x)[0]>Fraction('.9999') for x in value['sharper_jury']['margins']))
            check(filename+':complete block index '+block,value['count']==dict(unstable=1,neutral=0,stable=2))
    reintegration=read(DATA/'reintegration-640.json')
    algebra=read(DATA/'algebra-controls.json')
    check('fresh independent algebra controls',algebra['status']=='passed' and all(x['passed'] for x in algebra['checks']))
    check('fresh algebra executed source',sha(DATA/'algebra-controls.source.py')==sha(ROOT/'research/check_motor_five_proof.py'))
    check('new integration exact input',reintegration['input_sha256']==sha(cycle_path))
    check('new full47D embedding control',reintegration['independent_N16_embedding_control']['passed'])
    for name,expected in reintegration['source_sha256'].items():
        check('new integration source:'+name,sha(DATA/('reintegration-640.'+name))==expected==sha(ROOT/'research'/name))
    folder=DATA/'nonlinear';protocol=read(folder/'protocol.json');analysis=read(folder/'analysis.json')
    observed={(r['cluster'],r['amplitude'],r['sign']) for r in protocol['starts']}
    check('complete declared ensemble',len(protocol['starts'])==8 and observed==set(itertools.product((2,4),(1e-6,1e-8),(-1,1))))
    check('numerical source snapshots',sha(folder/'source.py')==protocol['source_sha256']==sha(ROOT/'research/illustrate_motor_replication.py')
          and sha(folder/'motor_rows.py')==protocol['field_sha256']==sha(ROOT/'research/motor_rows.py'))
    check('numerical status',analysis['status']=='passed' and len(analysis['runs'])==8)
    for start,record in zip(protocol['starts'],analysis['runs']):
        path=folder/f"trajectory-{record['index']:02}.npz";h=np.load(path)
        check('trajectory hash '+path.name,sha(path)==record['trajectory_sha256'])
        check('full-space shape and values '+path.name,h['states'].shape==(2401,47) and np.isfinite(h['states']).all()
              and np.array_equal(h['states'][0],start['start']))
        error=float(np.max(abs(h['main_states']-h['states'])))
        check('state agreement '+path.name,error==record['max_state_disagreement'] and error<2.69e-8)
        relative=h['observed']/h['prediction']-1
        check('return residual record '+path.name,np.array_equal(relative,np.array(record['relative_return_errors'])))
        if start['amplitude']==1e-8:
            check('reported small-signal agreement '+path.name,np.max(abs(relative[1:4]))<1e-6)
    radau=np.load(folder/'radau-control.npz')['states'];dop=np.load(folder/'trajectory-00.npz')['states']
    check('different solver control',sha(folder/'radau-control.npz')==analysis['radau_sha256']
          and float(np.max(abs(radau-dop)))==analysis['independent_radau_max_state_difference']<1.97e-8)
    visual=read(PAPER/'figures/manifest.json')
    check('renderer source',visual['renderer_sha256']==sha(ROOT/'research/render_motor_replication_article.py'))
    for key in ('input_sha256','outputs'):
        for name,expected in visual[key].items():
            check('figure '+key+':'+name,sha(ROOT/name)==expected)
    gif=Image.open(PAPER/'figures/cycle-rotation.gif')
    check('animation structure',gif.n_frames==96 and gif.info.get('loop')==0)
    pdf=PdfReader(PAPER/'manuscript.pdf');texts=[page.extract_text() for page in pdf.pages]
    check('PDF has expected complete document',len(pdf.pages)==10 and all(len(t)>300 for t in texts)
          and 'Theorem 1' in ''.join(texts) and 'References' in ''.join(texts))
    tex=(PAPER/'main.tex').read_text(encoding='utf-8')
    labels=set(re.findall(r'\\label\{([^}]+)\}',tex))
    check('LaTeX cross-references',set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',tex))<=labels)
    cites={x for group in re.findall(r'\\cite\{([^}]+)\}',tex) for x in group.split(',')}
    check('LaTeX bibliography',cites<=set(re.findall(r'\\bibitem\{([^}]+)\}',tex)))
    for path in [PAPER/'README.md',PAPER/'priority-and-positioning.md',PAPER/'research-programme.md',ROOT/'docs/motor-replication-index.md']:
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            check('local link:'+target,(path.parent/target.split('#')[0]).exists() or (args.write_records and target=='review.json'))
    if (PAPER/'manifest.json').exists():
        for item in read(PAPER/'manifest.json')['files']:
            check('article manifest:'+item['path'],sha(ROOT/item['path'])==item['sha256'])
    passed=all(item['passed'] for item in checks)
    report=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed' if passed else 'failed',
                check_count=len(checks),checks=checks,pdf_pages=len(pdf.pages),python=platform.python_version(),
                packages={name:importlib.metadata.version(name) for name in ['numpy','scipy','matplotlib','Pillow','pypdf']},
                proof_arithmetic='Python-FLINT 0.9.0 / Arb; 640-bit reintegration and 768-bit replay',
                visual_review=dict(completed=args.visual_review_complete,pages=list(range(1,len(pdf.pages)+1)) if args.visual_review_complete else [],
                                   scope='Assistant inspection of rendered PDF pages, all three figures and representative GIF frames; recorded separately from automated checks.'),
                scientific_scope='Complete periodic-orbit replication index; no least period, basin fraction, asymptotic chaos or historical-priority claim.')
    print(report['status'],len(checks),'checks',flush=True)
    for item in checks:
        if not item['passed']:print('FAILED',item['name'],flush=True)
    if not passed:raise SystemExit(2)
    if args.write_records:
        if not args.visual_review_complete:raise ValueError('Inspect the rendered document before packaging')
        for path in (PAPER/'review.json',PAPER/'manifest.json'):
            if path.exists():raise FileExistsError('Preserve the existing review packet')
        (PAPER/'review.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        paths=list(PAPER.rglob('*'))+list(DATA.rglob('*'))
        paths += [ROOT/'research'/name for name in ['certify_motor_replication_index.py','illustrate_motor_replication.py',
                                                 'render_motor_replication_article.py','audit_motor_replication_article.py']]
        paths += [ROOT/'docs/motor-replication-index.md',ROOT/'requirements-article.txt']
        paths=sorted({p for p in paths if p.is_file() and p.suffix not in ('.log','.aux','.out')})
        manifest=dict(created_utc=report['created_utc'],self_excluded=True,
                      scope='Research manuscript, new exact-index proof extension, full nonlinear ensemble and publication figures. Earlier packets remain separately frozen.',
                      files=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in paths])
        (PAPER/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
