"""Check the reviewed publication bytes and their frozen evidence bindings."""
import json
from pathlib import Path
from pypdf import PdfReader
from audit_raar_trapping_batch import file_hash as sha


def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-trapping';packet=repo/'results/raar-trapping'
    read=lambda p:json.loads(p.read_text(encoding='utf-8'))
    rows=read(paper/'manifest.json')['files']
    actual={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    assert actual=={r['path'] for r in rows}
    for row in rows:
        path=paper/row['path'];assert sha(path)==row['sha256'] and path.stat().st_size==row['bytes']
    build,review=read(paper/'build.json'),read(paper/'review.json')
    assert sha(paper/'manuscript.pdf')==build['pdf_sha256']==review['pdf_sha256']
    assert sha(paper/'content.py')==build['content_sha256']
    assert sha(repo/'research/build_raar_trapping_manuscript.py')==build['builder_sha256']
    assert sha(packet/'manifest.json')==build['evidence_manifest_sha256']
    assert sha(packet/'runs/raar-trapping-closure-audit01/report.json')==build['closure_report_sha256']
    pages=PdfReader(paper/'manuscript.pdf').pages
    assert len(pages)==build['pages']==build['planned_pages']==6
    assert review['visually_reviewed_pages']==list(range(1,7))
    assert [len(p.extract_text()) for p in pages]==build['page_text_lengths']
    assert all(len(p.extract_text())>1000 for p in pages)
    assert all(token in pages[3].extract_text() for token in ['21,405,802','23,055,541','235,677'])
    assert '0.01414' in pages[4].extract_text() and 'transitive' in pages[4].extract_text()
    for name in ['region','approach']:
        for ext in ['png','pdf']:
            file=f'raar-trapping-{name}.{ext}'
            assert sha(paper/'manuscript-figures'/file)==sha(packet/'runs/raar-trapping-render02'/file)
    print(json.dumps({'status':'passed','pages':len(pages),'pdf_sha256':build['pdf_sha256'],'evidence_manifest_sha256':build['evidence_manifest_sha256'],'scope':'Publication-byte, evidence-binding and reviewed-page accounting audit; visual review recorded separately.'}),flush=True)


if __name__=='__main__':main()
