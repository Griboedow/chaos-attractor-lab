"""Check the reviewed companion against its frozen evidence and exact figure bytes."""
import json
from decimal import Decimal
from pathlib import Path
from pypdf import PdfReader
from verify_rk23_controller import sha


def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-cycle-structure';root=repo/'results/raar-cycle-structure';js=lambda p:json.loads(p.read_text(encoding='utf-8'))
    manifest=js(paper/'manifest.json')
    for r in manifest['files']:
        p=paper/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],p
    assert {r['path'] for r in manifest['files']}=={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    b=js(paper/'build.json');review=js(paper/'review.json');assert sha(paper/'manuscript.pdf')==b['pdf_sha256']==review['pdf_sha256']
    assert sha(paper/'content.py')==b['content_sha256'] and sha(repo/'research/build_raar_cycle_manuscript.py')==b['builder_sha256'] and sha(root/'manifest.json')==b['evidence_manifest_sha256']
    for p in (paper/'manuscript-figures').iterdir():assert sha(p)==sha(root/'figures'/p.name)
    for r in js(root/'figures/sources.json')['sources']:assert sha(repo/r['path'])==r['sha256']
    pages=PdfReader(paper/'manuscript.pdf').pages;assert len(pages)==b['planned_pages']==6 and review['visually_reviewed_pages']==list(range(1,7))
    text=' '.join(' '.join(p.extract_text() for p in pages).split());assert '\ufffd' not in text and all(len(p.extract_text())>1000 for p in pages)
    for r in js(root/'cycle-proof/report.json')['cases']:assert repr(r['parameters'][0]) in text
    for r in js(root/'multiplier-proof/report.json')['cases']:
        assert format(Decimal(r['beta_display']),'.18f') in text and r['endpoints'][0]['stable_certified']
        assert r['uniform_branch_certificate']['other_two_inside_at_event_certified']
    fine=js(root/'fine-slice/report.json');assert fine['complete']==1025 and fine['positive']==90 and sum(v for k,v in fine['period_counts'].items() if k!='0')==933
    for p in [7,14,28,56,112,224]:assert f'{p}:{fine["period_counts"][str(p)]}' in text
    structural=js(root/'structure-audit/report.json');assert structural['unstable_real_multipliers']==1 and structural['stable_real_multipliers']==2
    for phrase in ['not prove an infinite cascade','historical uniqueness remains unresolved','not proof-assistant formalizations','negative determinant','two remain unclassified','nondegenerate flip','four attracting cycles are at different','do not express physical measurement accuracy','time-dependent forcing','2014','2026']:assert phrase in text,phrase
    print('PASS: six reviewed pages, exact manifests and figure bytes, cycle parameters, rounded event centers, fine-slice counts, orientation exclusion and explicit novelty limits')


if __name__=='__main__':main()
