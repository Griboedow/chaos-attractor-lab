"""Color subsequently classified regular starts in the frozen approach figure."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from render_twist import density


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    rows=json.loads((args.source/'analysis.json').read_text())['runs'];follow=json.loads(args.followup.read_text())
    regular=[r['index'] for r in follow['runs'] if r['kind']=='far_endpoint' and r['tail_period_up_to256']==68]
    assert len(regular)==436 and len(regular)==len(set(regular))
    canvas=Image.new('RGB',(2760,1270),'#071320');draw=ImageDraw.Draw(canvas)
    font=lambda size:ImageFont.truetype(findfont('DejaVu Sans'),size)
    draw.text((35,20),'FROM A VOLUME TO FOLDED SHEETS',font=font(43),fill='#e6edf5')
    draw.text((35,80),'32,768 initial states per row / same points followed through time / no interpolation',font=font(25),fill='#9eb5c9')
    histograms=[];regular_histograms=[]
    for i,row in enumerate(rows):
        d=np.load(args.source/f"case-{row['id']}.npz");q=(d['snapshots']-row['center'])@np.array(row['basis']);top=200+i*500
        mask=np.zeros(len(d['starts']),dtype=bool)
        if row['id']==4099:mask[regular]=True
        title='CROWN / 847' if row['id']==847 else 'SERRATED WINGS / 4099'
        draw.text((35,top-48),title,font=font(25),fill='#e6edf5')
        for j,t in enumerate(d['times']):
            h=density(q[j],np.array(row['projection']),430,row['halfwidth'])
            hr=density(q[j,mask],np.array(row['projection']),430,row['halfwidth'])
            assert int(h.sum())==len(d['starts']) and int(hr.sum())==int(mask.sum())
            intensity=np.minimum(np.log1p(h)/np.log1p(70),1)**.72
            fraction=np.divide(hr,h,out=np.zeros(h.shape),where=h>0)
            rgb=(np.array([82,206,240])[None,None,:]*(1-fraction[:,:,None])+np.array([255,167,75])[None,None,:]*fraction[:,:,None])*intensity[:,:,None]
            im=Image.fromarray(np.uint8(rgb));left=25+j*455;canvas.paste(im,(left,top));histograms.append(h);regular_histograms.append(hr)
            draw.text((left+12,top+12),f'n = {t:,}',font=font(25),fill='#cbd9e5')
    draw.text((35,1160),'Orange: 436 wing-box starts subsequently reach a checked 68-cycle or its reflected copy.',font=font(26),fill='#ffa74b')
    draw.text((35,1203),'Cyan: remaining starts, not individually certified chaotic. Each row keeps one camera and spatial scale.',font=font(25),fill='#9eb5c9')
    canvas.save(out/'adamax-approach.png');np.savez_compressed(out/'histograms.npz',total=histograms,regular=regular_histograms,regular_start_indices=regular)
    (out/'metadata.json').write_text(json.dumps(dict(source=str(args.source.resolve()),followup_sha256=hashlib.sha256(args.followup.read_bytes()).hexdigest(),
        regular_starts=436,regular_scope='All436far endpoints, followed for another50000warm+500000observed updates, meet the68-step closure criterion. The additional weak-volume record repeats one of these original starts.',
        colors=dict(remaining=[82,206,240],regular=[255,167,75]),mixing='Pixel-count-weighted mixture; intensity from total count, log1p with vmaxlog71,gamma.72.',
        point_data_unchanged=True,spatially_clipped=0),indent=2)+'\n')
    for name in ['render_adamax_approach_outcomes.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--followup',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
