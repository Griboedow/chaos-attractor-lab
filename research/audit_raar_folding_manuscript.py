"""Check the reviewed note's claims, exact source figures and paper manifest."""
import json
from pathlib import Path
from pypdf import PdfReader
from verify_raar_folding import sha

def main():
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-folding';root=repo/'results/raar-folding';load=lambda p:json.loads(p.read_text())
    files=load(paper/'manifest.json')['files']
    assert {r['path'] for r in files}=={p.relative_to(paper).as_posix() for p in paper.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    for r in files:p=paper/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==r['bytes']
    build=load(paper/'build.json');review=load(paper/'review.json');assert sha(paper/'manuscript.pdf')==build['pdf_sha256']==review['pdf_sha256'];assert sha(paper/'content.py')==build['content_sha256'];assert sha(repo/'research/build_raar_folding_manuscript.py')==build['builder_sha256'];assert sha(root/'manifest.json')==build['evidence_manifest_sha256']
    for name,folder in [('raar-cycle-unfolding','render'),('raar-cycle-seed-controls','render'),('raar-four-preimages','inverse-render')]:
        for ext in ['.pdf','.png']:assert sha(paper/'manuscript-figures'/(name+ext))==sha(root/folder/(name+ext))
    pages=PdfReader(paper/'manuscript.pdf').pages;assert len(pages)==build['planned_pages']==5;assert review['visually_reviewed_pages']==list(range(1,6));text=' '.join(' '.join(p.extract_text() for p in pages).split());assert '\ufffd' not in text
    for r in load(root/'seeds/report.json')['runs']:assert format(r['distance_quantiles'][-1][0],'.8f') in text
    four=[int(r['label'].split('-')[1]) for r in load(root/'inverse-proof/report.json')['cases'] if r['number_of_preimages']==4]
    assert len(four)==7 and all(str(i) in text for i in four)
    for phrase in ['61,440','64 preimages','18 have two','seven have four','degree eight','not a validated chart','Unequal sample densities','Historical uniqueness remains unresolved','not a proof-assistant formalization','supplement remains unreviewed','invertible restriction to a smaller invariant set remains possible','does not represent physical measurement accuracy']:assert phrase in text,phrase
    print('PASS: five reviewed pages, exact figures/manifests, all six control medians, exact four-preimage target list and explicit scientific limits')
if __name__=='__main__':main()
