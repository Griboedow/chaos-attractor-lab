"""Freeze a bounded RAAR comparison without changing the upstream reference packet."""
import argparse,json,shutil
from pathlib import Path
from compare_raar_catalog import sha

def main(args):
    repo=Path(__file__).resolve().parents[1];out=repo/'results/raar-catalog-comparison';out.mkdir(exist_ok=False)
    for src,name in [(args.comparison,'comparison'),(args.initial_render,'render-initial'),(args.render,'render')]:shutil.copytree(src,out/name)
    (out/'README.md').write_text('# Bounded RAAR geometry comparison\n\nThe existing million-state R1051 sample is compared with the 108 frozen regimes in results/attractor-novelty/reference. No reference system is reintegrated. Both shape-preserving and covariance-whitened representations retain full rankings, four resampling controls, temporal-half comparisons and three local alignments. The plotted points are actual stored states. The initial rendering and corrected caption are both retained.\n\nThis is a finite empirical comparison, not a proof of dynamical inequivalence, historical novelty or chaotic attraction. Local alignment does not certify an optimal residual. Run research/verify_raar_catalog.py. Reproduce with research/compare_raar_catalog.py --out NEW_DIRECTORY, then research/render_raar_catalog.py --source NEW_DIRECTORY --out NEW_RENDER_DIRECTORY.\n',encoding='utf-8')
    shutil.copy2(out/'render/raar-catalog-comparison.png',repo/'docs/assets/raar-catalog-comparison.png')
    files=[dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    (out/'manifest.json').write_text(json.dumps(dict(files=files),indent=2)+'\n',encoding='utf-8');print(len(files),'files frozen')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--comparison',type=Path,required=True);p.add_argument('--initial-render',type=Path,required=True);p.add_argument('--render',type=Path,required=True);main(p.parse_args())
