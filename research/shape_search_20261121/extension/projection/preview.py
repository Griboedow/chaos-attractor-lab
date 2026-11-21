"""Cache-only contact sheets in the earlier campaign's canonical PCA view."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'projection'))
from inspect_cases import picture


def main(args):
    records=[json.loads(p.read_text()) for p in args.records.glob('case-*.json')]
    keep=sorted([r for r in records if r['retained']],key=lambda r:-r['normalized_descriptor_distance_to_nearest_old_mean'])
    args.out.mkdir(parents=True,exist_ok=True)
    for n in range(0,len(keep),20):
        batch=keep[n:n+20]
        im=Image.new('RGB',(1000,166*((len(batch)+3)//4)),(8,18,26));d=ImageDraw.Draw(im)
        for j,row in enumerate(batch):
            path=args.cache_root/row['cloud_relative']
            with np.load(path) as f:z=f['points']
            x=j%4*250;y=j//4*166;im.paste(picture(z,246,145),(x,y+18))
            d.text((x+4,y+2),f"{row['case']}   d={row['normalized_descriptor_distance_to_nearest_old_mean']:.2f}   L3={row['rates'][2]:.3f}",fill='white')
        im.save(args.out/f'sheet-{n//20:02d}.jpg',quality=88)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--records',required=True,type=Path);p.add_argument('--cache-root',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    main(p.parse_args())
