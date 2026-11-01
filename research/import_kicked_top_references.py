"""Include the repository's frozen recomputations of the 1988 prior model."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np


def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=root/'results/attractor-novelty/reference')
    out=parser.parse_args().out;(out/'trajectories').mkdir(parents=True,exist_ok=True)
    path=root/'results/twist-followup/prior-comparison/spectra.npz'
    source=np.load(path)
    rows=[]
    for index,name in [(1,'KickedTop-A35-b0.1'),(2,'KickedTop-A15-b0.5')]:
        np.savez_compressed(out/'trajectories'/f'{name}.npz',points=source['tails'][index])
        rows.append({'name':name,'source_index':index,'parameters_in_twist_coordinates':source['parameters'][index].tolist(),
                     'start_in_twist_coordinates':source['starts'][index].tolist(),
                     'warm_steps':100000,'observed_steps':1000000,'retained_final_points':20000,
                     'spectrum':source['windows'][index].mean(0).tolist()})
    info={'source_path':str(path.relative_to(root)).replace('\\','/'),
          'source_sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest(),
          'source_commit':'689f1b75d2583010de732c4151550275c9be20b1',
          'report':'docs/twist-priority-and-mechanism.md',
          'coordinate_relation':'Uniform scale, rotation and translation conjugate the prior model to this twist-map slice.',
          'excluded':{'source_index':0,'A':-5.,'b':.5,'reason':'Fixed-point outcome; singular point cloud, not usable for shape normalization. Earlier follow-up retains the OCR-caption reproduction gap and all 33 fixed-point runs.'},
          'runs':rows}
    (out/'kicked-top-runs.json').write_text(json.dumps(info,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
