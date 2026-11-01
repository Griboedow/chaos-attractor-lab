"""Hash the completed checkpoint; no numerical calculations or source mutations."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/unexpected-search'
paths=[p for p in BASE.rglob('*') if p.is_file() and p.name!='manifest.json']
paths += [p for p in (ROOT/'research/unexpected_search').rglob('*') if p.is_file()]
paths += [ROOT/p for p in ['research/verify_unexpected_search.py','research/archive_unexpected_search.py','docs/tetrahedral-flow.md','docs/anderson-memory.md']]
paths += [ROOT/'docs/assets'/name for name in ['tetrahedral-flow.png','tetrahedral-flow-oblique.png','tetrahedral-vortex-565.png','tetrahedral-flow.gif','anderson-memory.png']]
items=[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(paths)]
(BASE/'manifest.json').write_text(json.dumps({'schema':1,'date':'2026-11-01','files':items},indent=2)+'\n',encoding='utf-8')
print(json.dumps({'files':len(items),'bytes':sum(i['bytes'] for i in items)}))
