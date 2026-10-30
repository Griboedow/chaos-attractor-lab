"""Verify the frozen network bundles and documentation links without edits."""
import ast
import datetime
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'manifest-final-02.json'
if OUTPUT.exists():
    raise FileExistsError('Preserve prior manifest')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

checks = []
for folder in (HERE, HERE / 'figure-layout-review'):
    record = json.loads((folder / 'analysis.json').read_text())
    checks.append(dict(check='run status', path=str(folder.relative_to(ROOT)),
                       passed=record['status'] == 'passed'))
    checks.append(dict(check='executed source', path=str(folder.relative_to(ROOT)),
                       passed=sha(folder / 'source.py') == record['source_sha256']))
    for name, digest in record['outputs'].items():
        checks.append(dict(check='output SHA256', path=str((folder/name).relative_to(ROOT)),
                           passed=sha(folder/name) == digest))
    for name, digest in record['inputs'].items():
        checks.append(dict(check='input SHA256', path=name, passed=sha(ROOT/name) == digest))
for name in ('exact-graph.json', 'phase-audit.json'):
    checks.append(dict(check='figure revision preserves science bytes', path=name,
                       passed=(HERE/name).read_bytes() == (HERE/'figure-layout-review'/name).read_bytes()))
script = ROOT/'research/motor_network_combinatorics.py'
ast.parse(script.read_text(encoding='utf-8'))
checks.append(dict(check='current script equals latest executed source',
                   passed=script.read_bytes() == (HERE/'figure-layout-review/source.py').read_bytes()))
report = ROOT/'docs/motor-network-combinatorics.md'
for link in re.findall(r'\]\(([^)]+)\)', report.read_text(encoding='utf-8')):
    target = (report.parent / link).resolve()
    checks.append(dict(check='report link', path=link, passed=target.exists()))
files = sorted([script, report] + [p for p in HERE.rglob('*') if p.is_file()])
record = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              status='passed' if all(x['passed'] for x in checks) else 'failed',
              checks=checks, files={str(p.relative_to(ROOT)):dict(sha256=sha(p), bytes=p.stat().st_size) for p in files},
              total_bytes=sum(p.stat().st_size for p in files),
              scope='Publication byte/link audit. Exact graph and numerical phase audit have separate scopes.')
OUTPUT.write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
print(record['status'], 'checks', len(checks), 'files', len(files), 'bytes', record['total_bytes'])
if record['status'] != 'passed':
    raise SystemExit(2)
