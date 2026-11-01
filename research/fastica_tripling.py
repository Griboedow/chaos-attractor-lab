"""Reproduce the structured FastICA study in an explicit external directory."""
from pathlib import Path
import argparse,os,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
STAGES={'pilot':'model.py','exact':'exact_check.py','controls':'controls.py',
        'long':'long_check.py','precision':'precision_check.py','library':'library_check.py'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--stage',choices=['all',*STAGES],default='all')
    args=parser.parse_args();out=args.output.resolve()
    if out.is_relative_to(ROOT):parser.error('Use a directory outside the checkout, normally on D:.')
    out.mkdir(parents=True,exist_ok=True)
    for source in (ROOT/'results/fastica-tripling/code').glob('*.py'):
        dest=out/source.name
        if dest.exists() and dest.read_bytes()!=source.read_bytes():
            raise ValueError(f'Refusing to replace changed source: {dest}')
        shutil.copyfile(source,dest)
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    env.setdefault('NUMBA_CACHE_DIR',str(out/'numba'))
    env.setdefault('NUMBA_NUM_THREADS','6')
    env.setdefault('OPENBLAS_NUM_THREADS','1');env.setdefault('MKL_NUM_THREADS','1')
    selected=STAGES if args.stage=='all' else {args.stage:STAGES[args.stage]}
    for name,source in selected.items():
        print(f'Running {name} in {out}',flush=True)
        subprocess.run([sys.executable,str(out/source)],cwd=out,env=env,check=True)

if __name__=='__main__':main()
