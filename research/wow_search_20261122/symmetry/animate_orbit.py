"""Show consecutive map iterates against a faint native-coordinate density background."""
import argparse
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib import font_manager


def render(cloud,path):
    with np.load(cloud) as data:
        xy=data['points']
    n=560
    pad=30
    s=(n-2*pad)/(2*abs(xy).max())
    x=np.rint(n/2+s*xy[:,0]).astype(int)
    y=np.rint(n/2-s*xy[:,1]).astype(int)
    count=np.bincount(y*n+x,minlength=n*n).reshape(n,n)
    v=np.clip(np.log1p(count)/np.log1p(np.percentile(count[count>0],99.9)),0,1)
    bg=np.array([7.,14.,25.])
    cyan=np.array([98.,191.,193.])
    im=Image.fromarray((bg+v[:,:,None]**.72*(cyan-bg)).astype(np.uint8))
    regular=ImageFont.truetype(font_manager.findfont('DejaVu Sans'),17)
    heading=ImageFont.truetype(font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold')),24)
    frames=[]
    for k in range(108):
        f=Image.new('RGB',(n,n+114),tuple(bg.astype(int)));f.paste(im,(0,60));d=ImageDraw.Draw(f)
        d.text((23,14),'One orbit, two alternating rings',font=heading,fill='#eaf1f2')
        for j in range(max(0,k-5),k+1):
            age=k-j
            radius=7 if age==0 else 3
            px,py=int(x[j]),int(y[j])+60
            col=('#ffd281' if age==0 else '#a39169')
            d.ellipse((px-radius,py-radius,px+radius,py+radius),fill=col,outline='#fff0d0' if age==0 else None,width=2)
        r=np.linalg.norm(xy[k])
        d.text((25,n+66),'Case 13745',font=regular,fill='#9fb1c6')
        d.text((205,n+66),f'Step {k+1:03d}',font=regular,fill='#f5e3c0')
        d.text((355,n+66),f'Radius {r:.3f}',font=regular,fill='#9fb1c6')
        frames.append(f)
    frames[0].save(path,format='WEBP',save_all=True,append_images=frames[1:],duration=125,loop=0,lossless=True,method=5)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cloud',type=Path);parser.add_argument('output',type=Path)
    opts=parser.parse_args();render(opts.cloud,opts.output)
