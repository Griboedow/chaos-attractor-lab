"""Compare the reviewed RK45 manuscript with the retained finite evidence."""
import json
from pathlib import Path
import numpy as np
from pypdf import PdfReader
from verify_rk23_controller import read,sha


def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/rk45-controller';root=repo/'results/rk45-controller'
    js=lambda p:json.loads(p.read_text(encoding='utf-8'))
    manifest=js(paper/'manifest.json')
    for row in manifest['files']:
        p=paper/row['path'];assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],p
    assert {r['path'] for r in manifest['files']}=={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    build=js(paper/'build.json');review=js(paper/'review.json')
    assert sha(paper/'manuscript.pdf')==build['pdf_sha256']==review['pdf_sha256']
    assert sha(paper/'content.py')==build['content_sha256']
    assert sha(repo/'research/build_rk45_manuscript.py')==build['builder_sha256']
    assert sha(root/'manifest.json')==build['evidence_manifest_sha256']
    pages=PdfReader(paper/'manuscript.pdf').pages;text=' '.join(' '.join(p.extract_text() for p in pages).split())
    assert len(pages)==build['planned_pages']==len(review['visually_reviewed_pages'])
    assert all(len(p.extract_text())>1000 for p in pages) and '\ufffd' not in text
    for row in js(paper/'figures/provenance.json')['sources']:assert sha(repo/row['path'])==row['sha256']
    scout=read(root/'oscillator/scout/screen.npz');follow=read(root/'oscillator/follow/follow.npz')
    assert len(scout['params'])==2051 and scout['flags'].sum()==1061
    assert list(np.flatnonzero(scout['counters'][:,0]))==[802]
    assert np.array_equal(follow['ids'],np.flatnonzero(scout['flags']))
    for val in [len(scout['params']),scout['flags'].sum(),follow['flags'].sum()]:assert str(val) in text
    d=read(root/'oscillator/five-million/follow.npz')
    for ident in [335,448,1890,2006]:
        rates=d['rates'][list(d['ids']).index(ident)]
        assert np.all(rates[:,0]>0) and f'{rates[:,0].min():.6f}' in text
        for value in rates.mean(0):assert f'{value:.6f}' in text
    minima=[d['rates'][list(d['ids']).index(335),:,1].min()]
    for digits in [70,90]:minima.append(np.array(js(root/f'oscillator/precision/case-335-{digits}.json')['rates'])[:,1].min())
    minima.append(read(root/'oscillator/upstream/case-335-analysis.npz')['rates'][:,1].min())
    assert min(minima)>0 and all(f'{v:.6f}' in text for v in minima)
    for r in js(root/'oscillator/error-norm-comparison/report.json')['runs']:
        for v in r['means']:assert f'{v:.6f}' in text
    cohort=read(root/'oscillator/cohort/cohort.npz')
    for ident,expected in [(335,[24,6,26]),(448,[24,6,29]),(2006,[24,6,32]),(1890,[24,6,32])]:
        for group,count,total in zip(['start','step','neighbor'],expected,[24,6,32]):
            mask=(cohort['ids']==ident)&(cohort['groups']==group)
            assert mask.sum()==total and cohort['flags'][mask].sum()==count and f'{count} / {total}' in text
        mask=(cohort['ids']==ident)&(cohort['groups']=='neighbor')
        assert not np.any(cohort['flags'][mask]&(cohort['rates'][mask,:,1].min(1)>0))
    for ident in [335,448]:
        d=read(root/f'oscillator/approach/case-{ident}.npz')
        assert d['snapshots'].shape==(6,32768,3) and np.all(d['status']==0)
        assert f'{np.median(d["distances"][-1]):.5f}' in text
    audit=js(root/'source/structure/analysis.json');assert audit['status']=='passed'
    for r in audit['selected_cases']:
        assert r['sturm_roots']==0 and r['gain_condition']==(r['id']!=1890)
        assert f'{r["squared_gain_decimal"]:.6g}' in text
        if r['gain_condition']:assert f'{r["normalized_P_bound_constant"]:.6f}' in text
    dr=js(root/'oscillator/damping/report.json')
    assert [(r['positive'],r['two_positive'],r['periodic']) for r in dr['groups']]==[(211,4,41),(194,125,63)]
    for r in dr['groups']:
        for key in ['count','positive','two_positive','periodic']:assert str(r[key]) in text
    phase=next(r for r in js(root/'oscillator/phase/report.json')['runs'] if r['id']==335)
    assert f'{phase["distance_from_half_turn_quantiles"][2]:.6f}' in text
    assert 'not a priority certificate' in text and 'Historical uniqueness remains unresolved' in text
    assert 'not claimed as a new mechanism' in text and 'not a certified basin' in text
    assert 'No compact global attractor' in text and 'not declared a proven regular attractor' in text
    assert 'Stuart' in text and 'Christodoulou' in text
    print(f'PASS: {len(pages)} reviewed pages, publication manifest, figure provenance, displayed rates/counts, exact-certificate constants, norm controls, damping slices and explicit priority limits')


if __name__=='__main__':main()
