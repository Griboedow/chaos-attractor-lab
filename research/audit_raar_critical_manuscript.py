"""Check reviewed publication bytes and their completed evidence bindings."""
import json
from pathlib import Path
from pypdf import PdfReader
from audit_raar_coverings import sha


def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-critical-membership';packet=repo/'results/raar-critical-membership'
    read=lambda p:json.loads(p.read_text())
    rows=read(paper/'manifest.json')['files']
    actual={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    assert actual=={r['path'] for r in rows}
    for row in rows:
        p=paper/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes']
    build,review=read(paper/'build.json'),read(paper/'review.json')
    assert sha(paper/'manuscript.pdf')==build['pdf_sha256']==review['pdf_sha256']
    assert sha(paper/'content.py')==build['content_sha256']
    assert sha(repo/'research/build_raar_critical_manuscript.py')==build['builder_sha256']
    assert sha(packet/'manifest.json')==build['evidence_manifest_sha256']
    assert sha(packet/'theorem.md')==build['theorem_sha256']
    for field,name in [('membership_audit_sha256','raar-critical-membership-audit03'),('secondary_audit_sha256','raar-critical-membership-audit04')]:
        p=packet/'runs'/name/'report.json';assert sha(p)==build[field] and read(p)['status']=='passed'
    pages=PdfReader(paper/'manuscript.pdf').pages
    assert len(pages)==build['pages']==build['planned_pages']==5
    assert review['visually_reviewed_pages']==list(range(1,6))
    assert [len(p.extract_text()) for p in pages]==build['page_text_lengths']
    assert 'Hausdorff' in pages[2].extract_text() and '0.20984' in pages[2].extract_text()
    assert 'locally injective' in pages[4].extract_text() and 'historical' in pages[4].extract_text().lower()
    for ext in ['png','pdf']:
        name=f'raar-critical-membership.{ext}'
        assert sha(paper/'manuscript-figures'/name)==sha(packet/'runs/raar-critical-membership-render02'/name)
    print(json.dumps({'status':'passed','pages':len(pages),'pdf_sha256':build['pdf_sha256'],'evidence_manifest_sha256':build['evidence_manifest_sha256']}),flush=True)


if __name__=='__main__':main()
