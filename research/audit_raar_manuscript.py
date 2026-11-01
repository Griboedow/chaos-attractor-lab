"""Compare the reviewed manuscript against the exact frozen numerical evidence."""
import json
from pathlib import Path
import numpy as np
from pypdf import PdfReader
from verify_rk23_controller import read,sha


def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-circles';root=repo/'results/raar-circles'
    js=lambda p:json.loads(p.read_text(encoding='utf-8'))
    manifest=js(paper/'manifest.json')
    for row in manifest['files']:
        p=paper/row['path'];assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],p
    assert {r['path'] for r in manifest['files']}=={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    build=js(paper/'build.json');review=js(paper/'review.json')
    assert sha(paper/'manuscript.pdf')==build['pdf_sha256']==review['pdf_sha256']
    assert sha(paper/'content.py')==build['content_sha256']
    assert sha(repo/'research/build_raar_manuscript.py')==build['builder_sha256']
    assert sha(repo/'research/build_raar_figures.py')==build['figure_builder_sha256']
    assert sha(root/'manifest.json')==build['evidence_manifest_sha256']
    pages=PdfReader(paper/'manuscript.pdf').pages;text=' '.join(' '.join(p.extract_text() for p in pages).split())
    assert len(pages)==build['planned_pages']==len(review['visually_reviewed_pages'])==12
    assert all(len(p.extract_text())>1000 for p in pages) and '\ufffd' not in text
    for row in js(paper/'manuscript-figures/sources.json')['sources']:assert sha(repo/row['path'])==row['sha256']
    d=read(root/'five-million/follow.npz')
    for i,rates in zip(d['ids'],d['rates']):
        assert f'R{i}' in text and f'{rates[:,0].min():.6f}' in text
        for value in rates.mean(0):assert f'{value:.6f}' in text
    series=[d['rates'][list(d['ids']).index(820)]]
    for digits in [70,90]:series.append(np.array(js(root/f'precision/case-820-{digits}.json')['rates']))
    series.append(read(root/'upstream/case-820.npz')['rates'])
    for rates in series:
        for value in [rates[:,1].mean(),rates[:,1].min()]:assert f'{value:.6f}' in text
    assert all(r[:,1].mean()>0 for r in series) and series[-1][:,1].min()<0
    for i in [1051,1483,820,75]:
        data=read(root/f'approach/case-{i}.npz')
        for value in [np.median(data['distances'][0]),np.median(data['distances'][-1]),np.median(data['reference_split_distances'])]:assert f'{value:.6g}' in text
    for r in js(root/'relaxation/report.json')['groups']:
        for key in ['count','positive','two_positive','periodic']:assert str(r[key]) in text
    for r in js(root/'fixed-proof/report.json')['cases']:assert r['number_of_fixed_points']==(2 if r['id']==1679 else 0) and r['all_fixed_points_unstable']
    for phrase in ['one negative','historical uniqueness remain unresolved','not proof that it is unprecedented','not a proof of structural stability or a certified bifurcation diagram','not a classification of the attracting cloud','not proof-assistant formalizations','2018','385–403']:assert phrase.lower() in text.lower(),phrase
    print('PASS: 12 reviewed pages, exact manifests and figure provenance, displayed rates/distances, fixed-point counts, parameter slices and explicit novelty limits')


if __name__=='__main__':main()
