"""Audit the publication copy, frozen numerical claims, and figure input hashes."""
import json
from pathlib import Path
import numpy as np
from pypdf import PdfReader
from verify_optimizer_followups import arrays,sha


def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/yogi-memory';root=repo/'results/yogi-tanh'
    read=lambda p:json.loads(p.read_text())
    manifest=read(paper/'manifest.json')
    for r in manifest['files']:
        p=paper/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256']
    expected={r['path'] for r in manifest['files']}
    assert expected=={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    build=read(paper/'build.json');review=read(paper/'review.json');pdf=paper/'manuscript.pdf'
    assert sha(pdf)==build['pdf_sha256']==review['pdf_sha256']
    assert sha(paper/'content.py')==build['content_sha256'] and sha(repo/'research/build_yogi_manuscript.py')==build['builder_sha256']
    reader=PdfReader(pdf);text='\n'.join(p.extract_text() for p in reader.pages)
    assert len(reader.pages)==11 and len(review['visually_reviewed_pages'])==11
    assert '\ufffd' not in text and all(len(p.extract_text())>1500 for p in reader.pages)
    for r in read(paper/'figures/provenance.json')['sources']:assert sha(repo/r['path'])==r['sha256']
    d=arrays(root/'five-million/histories.npz');precision=read(root/'precision/analysis.json')
    for ident in [39,7239]:
        w=d['windows'][list(d['ids']).index(ident)];assert np.all(w[:,:2]>0)
        for value in [*w.mean(0),w[:,1].min()]:assert f'{value:.6f}' in text
        for digits in [70,90]:
            row=next(r for r in precision['runs'] if r['id']==ident and r['digits']==digits)
            assert f"{row['minimum_windows'][1]:.6f}" in text
        folder='upstream-analysis' if ident==39 else 'upstream-strong-analysis'
        row=next(r for r in read(root/folder/'analysis.json')['runs'] if r['id']==ident)
        assert f"{row['window_minima'][1]:.6f}" in text
    co=read(root/'cohorts/analysis.json')
    for row in co['summary']:
        good=row['two'] if row['id'] in [39,7239] else row['leading']
        assert f"{good} / {row['total']}" in text
    controls=read(root/'structure/analysis.json')
    assert controls['completed']==1024 and controls['adam_epsilon_outside_local_derivative_identity'] and controls['derivative_equals_relative_matrix']
    assert 'not established' in text and 'known Adam' in text and 'settles' in text and 'research draft' in text
    print('PASS: 11 publication pages, reviewed PDF hash, every manuscript file, exact figure inputs, displayed rates/cohort counts, known-equivalence controls and explicit claim limits')


if __name__=='__main__':main()
