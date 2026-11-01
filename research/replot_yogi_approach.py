"""Replot all approach states with explicit final zooms and competing outcomes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from render_twist import density,colorize


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    analysis=json.loads((args.source/'analysis.json').read_text());previews=json.loads((args.preview/'analysis.json').read_text())['runs']
    canvas=Image.new('RGB',(3420,1550),'#071320');draw=ImageDraw.Draw(canvas)
    large=ImageFont.truetype(findfont('DejaVu Sans'),49);small=ImageFont.truetype(findfont('DejaVu Sans'),25);tiny=ImageFont.truetype(findfont('DejaVu Sans'),23)
    draw.text((40,24),'YOGI / A VOLUME APPROACHES TWO COMPETING OUTCOMES',font=large,fill='#e8f0f5')
    draw.text((40,103),'32,768 initial states per row / columns 1-5 share one spatial scale / last column is an explicit zoom',font=small,fill='#9ab0c1')
    records=[]
    for i,row in enumerate(analysis['runs']):
        ident=row['grid_id'];d=np.load(args.source/f'case-{ident}.npz');z=d['snapshots'];fixed=d['fixed_tail'];q=(z-np.array(row['center']))@np.array(row['basis']);projection=np.array(row['projection'])
        preview=next(r for r in previews if r['id']==ident);zoom=float(preview['halfwidth']);full=row['halfwidth'];top=230+650*i
        draw.text((40,top-48),f"G{ident:04d} / {int(fixed.sum()):,} detected fixed tails, {len(fixed)-int(fixed.sum()):,} other endpoints",font=small,fill='#e8f0f5')
        for j in range(6):
            k=min(j,4);halfwidth=zoom if j==5 else full;points=q[k];h=density(points,projection,535,halfwidth);count=int(h.sum())
            if j<5:assert count==len(points)
            im,_=colorize(h,vmax=np.log1p(45));left=30+565*j;canvas.paste(im,(left,top))
            draw.text((left+13,top+12),f'n = {int(d["times"][k]):,}'+(' / ZOOM' if j==5 else ''),font=small,fill='#e8f0f5')
            if j>=4:
                xy=points[fixed]@projection.T
                for a,b in xy:
                    px=int(np.floor((a/halfwidth+1)*535/2));py=534-int(np.floor((b/halfwidth+1)*535/2))
                    if 0<=px<535 and 0<=py<535:draw.point((left+px,top+py),fill='#69d9e5')
            if j==5:
                draw.text((left+13,top+53),f'{full/zoom:.1f}x magnification',font=tiny,fill='#9ab0c1')
                draw.text((left+13,top+535-38),f'{len(points)-count:,} endpoints outside zoom',font=tiny,fill='#9ab0c1')
            records.append(dict(grid_id=ident,column=j,time=int(d['times'][k]),halfwidth=halfwidth,points=len(points),included=count,excluded=len(points)-count,zoom=j==5))
    draw.text((40,1462),'Cyan marks detected fixed tails at the final time; all endpoints are retained in the five full-view panels.',font=small,fill='#9ab0c1')
    draw.text((40,1502),'Finite ensemble experiment: no basin probability or infinite-time attraction theorem. No connecting curves or interpolated motion.',font=tiny,fill='#9ab0c1')
    canvas.save(out/'yogi-grid-approach.png')
    (out/'analysis.json').write_text(json.dumps(dict(source=str(args.source),panels=records,
        source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.source.glob('*.npz')},
        correction='Layout-only follow-up: increase row-label clearance, add explicitly scaled final zooms, overlay numerical fixed-tail outcomes. No dynamics rerun and no full-view point discarded.'),indent=2)+'\n')
    for n in ['replot_yogi_approach.py','render_twist.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['source','preview','output']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
