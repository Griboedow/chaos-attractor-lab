"""Resolve run storage and initialize a replay from the published exact inputs."""
from pathlib import Path
import os,shutil
ROOT=Path(__file__).resolve().parents[2]
def output_root():
    return Path(os.environ.get('CHAOS_FASTICA_OUTPUT','D:/Codex/fastica-reproduction'))
if __name__=='__main__':
    out=output_root();out.mkdir(parents=True,exist_ok=False)
    (out/'follow-01').mkdir()
    for p in (ROOT/'results/fastica/follow-01').glob('*-input.json'):
        shutil.copyfile(p,out/'follow-01'/p.name)
    print(out)
