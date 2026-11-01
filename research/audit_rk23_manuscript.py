"""Check the reviewed publication copy and its displayed finite-data claims."""
import json
from pathlib import Path
import numpy as np
from pypdf import PdfReader
from verify_rk23_controller import read,sha


def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/rk23-controller';root=repo/'results/rk23-controller'
    js=lambda p:json.loads(p.read_text(encoding='utf-8'))
    manifest=js(paper/'manifest.json')
    for row in manifest['files']:
        p=paper/row['path'];assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],p
    assert {r['path'] for r in manifest['files']}=={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    build=js(paper/'build.json');review=js(paper/'review.json')
    assert sha(paper/'manuscript.pdf')==build['pdf_sha256']==review['pdf_sha256']
    assert sha(paper/'content.py')==build['content_sha256']
    assert sha(repo/'research/build_rk23_manuscript.py')==build['builder_sha256']
    pages=PdfReader(paper/'manuscript.pdf').pages;text='\n'.join(p.extract_text() for p in pages)
    assert len(pages)==11 and len(review['visually_reviewed_pages'])==11
    assert all(len(p.extract_text())>1400 for p in pages) and '\ufffd' not in text
    for row in js(paper/'figures/provenance.json')['sources']:assert sha(repo/row['path'])==row['sha256']
    scout=read(root/'oscillator/scout/screen.npz');follow=read(root/'oscillator/follow/follow.npz')
    assert len(scout['params'])==1027 and scout['flags'].sum()==396 and follow['flags'].sum()==393
    assert sorted(follow['period'][follow['period']>0])==[14,19,172]
    for folder,ids,direction in [('five-million',[381,354,883],0),('hyper-five-million',[37,526,689],1)]:
        d=read(root/'oscillator'/folder/'follow.npz')
        for ident in ids:
            rates=d['rates'][list(d['ids']).index(ident)]
            assert np.all(rates[:,direction]>0) and f'{rates[:,direction].min():.6f}' in text
            if direction==0:
                for value in rates.mean(0):assert f'{value:.6f}' in text
    for ident in [37,526,689]:
        minima=[]
        for digits in [70,90]:minima.append(np.array(js(root/f'oscillator/hyper-precision/case-{ident}-{digits}.json')['rates'])[:,1].min())
        minima.append(read(root/f'oscillator/hyper-upstream/case-{ident}-analysis.npz')['rates'][:,1].min())
        assert all(f'{v:.6f}' in text for v in minima)
        assert (min(minima)>0)==(ident==526)
    for r in js(root/'oscillator/error-norm-comparison/report.json')['runs']:
        assert max(r['means'])<0
        for v in r['means']:assert f'{v:.6f}' in text
    cohort=read(root/'oscillator/cohort/cohort.npz')
    for ident,expected in [(381,[24,6,16]),(354,[24,6,32]),(883,[24,6,24])]:
        for group,count,total in zip(['start','step','neighbor'],expected,[24,6,32]):
            mask=(cohort['ids']==ident)&(cohort['groups']==group)
            assert mask.sum()==total and cohort['flags'][mask].sum()==count and f'{count} / {total}' in text
    for ident in [381,883]:
        d=read(root/f'oscillator/approach/case-{ident}.npz')
        assert d['snapshots'].shape==(6,32768,3) and np.all(d['status']==0)
        assert f'{np.median(d["distances"][-1]):.5f}' in text
    audit=js(root/'source/rk23-structure-audit.json');special=audit['special_case']
    assert audit['status']=='passed' and special['sturm_roots_on_closed_interval']==0
    assert f'{special["squared_controller_gain_decimal"]:.9f}' in text
    assert f'{special["normalized_P_bound_constant"]:.6f}' in text
    assert 'historical uniqueness remain unresolved' in text and 'not a priority certificate' in text
    assert 'mixed' in text and 'not a certified basin' in text and 'no compact global attractor' in text
    print('PASS: 11 reviewed pages, complete publication manifest, exact figure inputs, displayed rates and counts, mixed independent outcomes, exact-certificate constants and explicit priority limits')


if __name__=='__main__':main()
