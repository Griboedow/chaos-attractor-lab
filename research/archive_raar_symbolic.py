"""Freeze every declared proposal, failed pilot, completed covering and figure."""
import hashlib,json,shutil
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    repo=Path(__file__).resolve().parents[1];base=Path('D:/Codex/chaos-attractor-search-20261031');root=repo/'results/raar-symbolic-dynamics';root.mkdir(exist_ok=False)
    names=['raar-return-search01','raar-covering-proposal01','raar-covering-refined01','raar-covering-proof-pilot01']+[f'raar-covering-proof-edge{i}-01' for i in range(1,8)]+['raar-covering-ad-edge0-01']+[f'raar-covering-ad-edge{i}-02' for i in range(8)]+['raar-covering-audit-partial01','raar-covering-audit-primary01','raar-covering-audit-ad01','raar-symbolic-orbits01','raar-symbolic-cycle-proof01','raar-symbolic-cycle-audit01','raar-symbolic-render01','raar-symbolic-render02']
    provenance=[]
    for name in names:
        for p in sorted((base/name).rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts:continue
            q=root/'runs'/name/p.relative_to(base/name);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q);assert sha(p)==sha(q);provenance.append(dict(source=str(p),path=q.relative_to(root).as_posix(),sha256=sha(q)))
    source=root/'source';source.mkdir()
    for name in ['search_raar_returns','propose_raar_coverings','refine_raar_coverings','certify_raar_coverings','check_raar_coverings_ad','audit_raar_coverings','find_raar_symbolic_orbits','certify_raar_cycles','audit_raar_symbolic_cycle','continue_raar_cycles','trace_raar_cycle_seeds','search_raar_circles','render_twist','render_raar_symbolic','archive_raar_symbolic']:
        shutil.copy2(repo/'research'/(name+'.py'),source/(name+'.py'))
    (root/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    review=dict(date='2026-11-01',method_source=dict(authors=['Piotr Zgliczynski','Marian Gidea'],title='Covering relations for multidimensional dynamical systems',journal='Journal of Differential Equations',year=2004,volume=202,pages='32-58',doi='10.1016/j.jde.2004.03.013',url='https://ww2.ii.uj.edu.pl/~zgliczyn/papers/wts/cov.pdf',pdf_sha256=sha(base/'raar-covering-literature/covering-2004.pdf'),review_scope='Author preprint dated 20 April 2004: Theorem 9, page 11 (closed chains yield periodic points), and Theorem 16, page 17 (one exit direction), text and formulas visually inspected. Corollary 12 and Section 5.2 read for symbolic context. Third-party full text not redistributed.'),queries=['"RAAR" "topological entropy"','"Douglas Rachford" "covering relations"','"relaxed averaged alternating reflections" "symbolic dynamics"','"Douglas-Rachford" horseshoe chaos'],outcome='No exact spatial-circle covering result located in the returned results of this bounded search. Not a priority certificate. RAAR, projection-search chaos and the covering method are established prior work.',open_review_gaps=['Elser 2007 supplementary information','Elser 2025 book Chapter 5 full text'],other_reviews=['results/raar-circles/priority/review.json','results/raar-cycle-structure/priority/review.json','results/raar-folding/priority/review.json','results/raar-priority-extension/review.json'])
    (root/'priority.json').write_text(json.dumps(review,indent=2)+'\n',encoding='utf-8')
    for name in ['geometry','graph','itineraries']:shutil.copy2(root/'runs/raar-symbolic-render02'/f'raar-symbolic-{name}.png',repo/'docs/assets'/f'raar-symbolic-{name}.png')
    (root/'README.md').write_text('''# Certified symbolic dynamics in spatial-circle RAAR

The eight edges of a branching covering graph pass two interval map implementations. Standard covering theory gives a compact invariant set with topological entropy at least log(2)/49 per original iteration. This is not a certificate of an attracting chaotic set or historical uniqueness.

All runs are preserved under their original directory names in `runs`; source paths inside reports record provenance. The primary certificates are `raar-covering-proof-pilot01` and `raar-covering-proof-edge1-01` through `edge7-01`. The second implementation is `raar-covering-ad-edge0-02` through `edge7-02`. Both use the identical exact input from `raar-covering-refined01/input.json`.

The initial 36 numerical window proposals all failed their nonlinear sample tests. The refined proposal succeeds numerically and is then certified. The interrupted AD pilot `raar-covering-ad-edge0-01` is incomplete; its gzip stream may lack a footer and is never counted as a certificate. The partial structural audit is retained separately from both complete audits.

The two map implementations use analytic vector projectors versus scalar nested reflectors with automatic derivatives, and different matrix inversion methods. They share the subdivision driver and mpmath interval library. These are computer-assisted interval proofs, not proof-assistant formalizations.

The located least-period-49 saddle has a separate 200-digit existence/least-period certificate, a certified window itinerary, and exactly one unstable multiplier. Longer coded multiple-shooting examples are numerical illustrations. Their binary64 monodromy products are ill-conditioned: small eigenvalues must not be used as stability evidence.

The return scout's legacy `max_excursion` field measures distance from one selected phase of the seven-cycle, not distance from the whole cycle. It is not used as a geometric excursion estimate. The final itinerary figure uses nearest distance to all seven cycle points.

Final figure files are in `runs/raar-symbolic-render02`. Earlier renderings are retained. Run `python research/verify_raar_symbolic.py` from the repository root. See [the report](../../docs/raar-symbolic-dynamics.md) and [manuscript](../../docs/raar-symbolic-dynamics/manuscript.md).
''',encoding='utf-8')
    files=[dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(root.rglob('*')) if p.is_file()]
    (root/'manifest.json').write_text(json.dumps(dict(files=files),indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(files=len(files),bytes=sum(r['bytes'] for r in files),manifest_sha256=sha(root/'manifest.json'))))

if __name__=='__main__':main()
