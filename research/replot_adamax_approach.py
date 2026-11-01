"""Improve spacing using frozen approach data; no trajectories are recomputed."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from render_twist import density,colorize


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    meta=json.loads((args.source/'analysis.json').read_text());rows=meta['runs'];sources={}
    canvas=Image.new('RGB',(2760,1220),'#071320');draw=ImageDraw.Draw(canvas)
    large=ImageFont.truetype(findfont('DejaVu Sans'),43);small=ImageFont.truetype(findfont('DejaVu Sans'),25)
    draw.text((35,20),'FROM A VOLUME TO FOLDED SHEETS',font=large,fill='#e6edf5')
    draw.text((35,80),'32,768 initial states per row / same points followed through time / no interpolation',font=small,fill='#9eb5c9')
    for i,row in enumerate(rows):
        path=args.source/f"case-{row['id']}.npz";d=np.load(path);sources[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        q=(d['snapshots']-row['center'])@np.array(row['basis']);top=200+i*500
        title='CROWN / 847' if row['id']==847 else 'SERRATED WINGS / 4099'
        draw.text((35,top-48),title,font=small,fill='#e6edf5')
        for j,t in enumerate(d['times']):
            h=density(q[j],np.array(row['projection']),430,row['halfwidth']);assert int(h.sum())==len(d['starts'])
            im,_=colorize(h,vmax=np.log1p(70));left=25+j*455;canvas.paste(im,(left,top))
            draw.text((left+12,top+12),f'n = {t:,}',font=small,fill='#cbd9e5')
    draw.text((35,1170),'Each row uses one fixed spatial scale at every time. Final clusters remain discrete iterates, not drawn connecting curves.',font=small,fill='#9eb5c9')
    canvas.save(out/'adamax-approach.png')
    (out/'metadata.json').write_text(json.dumps(dict(source=str(args.source.resolve()),source_hashes=sources,
        change='Increase vertical spacing between rows by 40 pixels. Same point data, cameras, density mapping and scale.'),indent=2)+'\n')
    for name in ['replot_adamax_approach.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
