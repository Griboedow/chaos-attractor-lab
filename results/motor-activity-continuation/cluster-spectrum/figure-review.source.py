"""Improve adjacent tick readability while preserving the first frozen figure."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
out = root / 'figure-review'
if out.exists():
    raise FileExistsError('Preserve existing figure review')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
audit = json.loads((root / 'analysis/analysis.json').read_text())
record = json.loads((root / 'lane-46/analysis.json').read_text())
assert audit['status'] == 'passed' and sha(root / 'lane-46/analysis.json') == audit['run_record_sha256']
q = np.array([w['quotient_spectrum'] for w in record['windows']])
c = np.array([w['contrast_spectra'] for w in record['windows']])
times = np.arange(1, len(q) + 1) * 500
spectrum = np.r_[record['quotient_spectrum'], np.array(record['contrast_spectra']).ravel()]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
axes[0].plot(times, q[:, 0], 'o-', label='Quotient QR 1', color='#A85636')
axes[0].plot(times, q[:, 1], 'o-', label='Quotient QR 2', color='#747B85')
for k, pair in enumerate(([1, 2], [3, 6], [4, 7])):
    axes[0].plot(times, c[:, k, 0], '.-', label=f'Pair {pair[0]},{pair[1]} leading')
axes[0].axhline(0, color='black', linewidth=.8)
axes[0].set(xlabel='Observed tau (500-tau windows)', ylabel='Finite-window exponent',
            title='Expansion within clusters; transverse contraction')
axes[0].legend(fontsize=8, ncol=2)
axes[0].grid(alpha=.18)
axes[1].bar(np.arange(1, 24), spectrum, color=['#A85636'] * 14 + ['#477D99'] * 9)
axes[1].axhline(0, color='black', linewidth=.8)
axes[1].axvline(14.5, color='#777777', linewidth=.8, linestyle='--')
axes[1].set(xlabel='14 quotient + 9 pair-contrast directions', ylabel='Cumulative finite-time exponent',
            title='All 23 directions after 4,000 observed tau')
axes[1].set_xticks([1, 7, 14, 18, 23])
axes[1].grid(axis='y', alpha=.18)
fig.suptitle('Lane 46: declared 1+2+2+2+1 cluster restriction', fontsize=12)
fig.tight_layout()
out.mkdir()
fig.savefig(out / 'spectrum.png', dpi=170)
plt.close(fig)
review = dict(status='passed', created_utc=datetime.now(timezone.utc).isoformat(), source_sha256=sha(Path(__file__)),
              input_analysis_sha256=sha(root / 'analysis/analysis.json'), figure_sha256=sha(out / 'spectrum.png'),
              change='Remove adjacent 14/15 tick-label crowding. All plotted values unchanged; original figure and inventory retained.')
(out / 'analysis.json').write_text(json.dumps(review, indent=2) + '\n', encoding='utf-8', newline='\n')
inventory = {p.relative_to(root).as_posix(): dict(bytes=p.stat().st_size, sha256=sha(p))
             for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
(root / 'manifest-final.json').write_text(json.dumps(dict(status='passed', files=inventory,
    file_count=len(inventory), checks=audit['checks']), indent=2) + '\n', encoding='utf-8', newline='\n')
print(json.dumps(review))
