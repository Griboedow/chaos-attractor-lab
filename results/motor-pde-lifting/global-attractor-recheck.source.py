"""Repeat the algebra and link the frozen prerequisite used by the first note."""
import hashlib
import json
from pathlib import Path
import runpy
import sys

HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'global-attractor-algebra-recheck.json'
if OUTPUT.exists():
    raise FileExistsError('Preserve the previous recheck')
snapshot = HERE/'global-attractor-prerequisite.md'
original = json.loads((HERE/'global-attractor-algebra.json').read_text())
snapshot_hash = hashlib.sha256(snapshot.read_bytes()).hexdigest()
assert snapshot_hash == original['prerequisite_note_sha256']
sys.argv = [str(HERE/'global-attractor-algebra.py'), '--output', str(OUTPUT)]
runpy.run_path(sys.argv[0], run_name='__main__')
record = json.loads(OUTPUT.read_text())
record.update(original_prerequisite_snapshot='results/motor-pde-lifting/global-attractor-prerequisite.md',
    original_prerequisite_snapshot_sha256=snapshot_hash,
    original_record_sha256=hashlib.sha256((HERE/'global-attractor-algebra.json').read_bytes()).hexdigest(),
    recheck_driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    prerequisite_note='The original record references the frozen snapshot. This recheck also records the current prerequisite note hash; later editorial changes do not replace the original evidence.')
OUTPUT.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
