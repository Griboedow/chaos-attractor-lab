"""Audit the reviewed article against frozen certificates and figure bytes."""
import json
from pathlib import Path
from pypdf import PdfReader
from archive_raar_symbolic import sha

def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-symbolic-dynamics';root=repo/'results/raar-symbolic-dynamics';js=lambda p:json.loads(p.read_text(encoding='utf-8'))
    rows=js(paper/'manifest.json')['files'];assert {r['path'] for r in rows}=={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    for r in rows:p=paper/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==r['bytes']
    b=js(paper/'build.json');review=js(paper/'review.json');assert sha(paper/'manuscript.pdf')==b['pdf_sha256']==review['pdf_sha256'];assert sha(paper/'content.py')==b['content_sha256'];assert sha(repo/'research/build_raar_symbolic_manuscript.py')==b['builder_sha256'];assert sha(root/'manifest.json')==b['evidence_manifest_sha256']
    pages=PdfReader(paper/'manuscript.pdf').pages;assert len(pages)==b['planned_pages']==7 and review['visually_reviewed_pages']==list(range(1,8));text=' '.join(' '.join(p.extract_text() for p in pages).split());assert '\ufffd' not in text and all(len(p.extract_text())>1500 for p in pages)
    for p in (paper/'manuscript-figures').iterdir():assert sha(p)==sha(root/'runs/raar-symbolic-render02'/p.name)
    raw=js(root/'runs/raar-covering-refined01/input.json')
    for v in raw['params']:assert repr(v) in text
    for phrase in ['16,467','9,449','0.01414','579.0979','579.0981','0.10367','Theorem 9','Theorem 16','not a wholly independent software stack','historical uniqueness remain unresolved','exact least period 49','ill-conditioned','not a validated invariant manifold','not individually root-certified','No third-party publication PDF is redistributed','Different equations or appearances do not exclude nonlinear conjugacy']:assert phrase in text,phrase
    assert js(root/'runs/raar-covering-audit-primary01/report.json')['complete_graph'];assert js(root/'runs/raar-covering-audit-ad01/report.json')['complete_graph']
    print('PASS: seven visually reviewed pages, exact evidence and figure bytes, parameter table, proof counts, entropy bound, located saddle and explicit attraction/priority limits')

if __name__=='__main__':main()
