"""Explain the crown by grouping actual states by time since maximum-memory reset."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from render_twist import camera,density,colorize


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    z=np.load(args.source/'case-847.npy',mmap_mode='r');meta=json.loads((args.source/'analysis.json').read_text())['runs']
    row=next(r for r in meta if r['id']==847);b,c,e=row['parameters'];ages=np.full(len(z),-1,dtype=np.int16)
    age=-1
    for i in range(1,len(z)):
        if abs(z[i-1,0])+e>=c*z[i-1,2]:age=0
        elif age>=0:age+=1
        ages[i]=age
    q=(z-row['center'])@np.array(row['pca_basis_columns']);projection=camera(np.pi/4,.25)
    groups=[ages==0,ages==1,ages>=2];labels=['RESET / AGE 0','ONE MEMORY STEP','TWO TO FOUR MEMORY STEPS']
    canvas=Image.new('RGB',(2100,980),'#071320');draw=ImageDraw.Draw(canvas)
    font=lambda size:ImageFont.truetype(findfont('DejaVu Sans'),size)
    draw.text((40,22),'THE MAXIMUM MEMORY CREATES SHEETS',font=font(46),fill='#e6edf5')
    draw.text((40,92),'Reset identity:  M = sign(x) U (U - e) - X U',font=font(33),fill='#ffc07d')
    draw.text((40,146),'Case 847 / actual states / one rigid camera and spatial scale / no connecting lines',font=font(27),fill='#9eb5c9')
    histograms=[];counts=[]
    for k,(mask,label) in enumerate(zip(groups,labels)):
        left=35+690*k;h=density(q[mask],projection,650,row['halfwidth']);count=int(mask.sum());assert int(h.sum())==count
        im,_=colorize(h,vmax=np.log1p(100));canvas.paste(im,(left,270));counts.append(count);histograms.append(h)
        draw.text((left+10,210),label,font=font(27),fill='#e6edf5')
        draw.text((left+10,248),f'{count:,} states',font=font(22),fill='#9eb5c9')
    draw.text((40,937),'After each reset, memory-only updates transport the same two quadratic surfaces into additional layers.',font=font(25),fill='#9eb5c9')
    canvas.save(out/'adamax-sheets.png');np.savez_compressed(out/'histograms.npz',histograms=histograms,counts=counts)
    result=dict(id=847,source_sha256=hashlib.sha256((args.source/'case-847.npy').read_bytes()).hexdigest(),
        source_metadata_sha256=hashlib.sha256((args.source/'analysis.json').read_bytes()).hexdigest(),
        counts=counts,unknown_initial_age_count=int((ages<0).sum()),maximum_observed_age=int(ages.max()),
        projection=projection.tolist(),halfwidth=row['halfwidth'],center=row['center'],basis=row['pca_basis_columns'],
        resolution=650,palette='magma',log_density_vmax=float(np.log1p(100)),gamma=.72,spatially_clipped=0)
    assert sum(counts)+result['unknown_initial_age_count']==len(z)
    (out/'metadata.json').write_text(json.dumps(result,indent=2)+'\n')
    for name in ['render_adamax_sheets.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    print(result['counts'],'unknown initial age',result['unknown_initial_age_count'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
