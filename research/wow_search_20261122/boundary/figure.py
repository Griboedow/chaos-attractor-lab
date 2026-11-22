"""Render the finite-window forward basins and independently sampled inverse boundary."""
import argparse,json,math
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from scipy.ndimage import maximum_filter,gaussian_filter
from search import plane,roots,RAW,OUTPUT
from inverse import refined
from selected import D,K,BOUNDS
# Darker, cool root colors; a warm sequential palette marks the nine cycle phases.
COLORS=['#186b7a','#398f81','#637db2','#56558d','#7986a7','#237e97',
        '#9f4e54','#bd5350','#d3674c','#e27e4e','#e79955','#eab45f','#eacb77','#d99287','#be7290']
def rgb(h):return tuple(int(h[a:a+2],16) for a in (1,3,5))
def shaded(labels,steps,edge=True):
    pal=np.array([rgb('#111e2b')]+[rgb(v) for v in COLORS])
    col=pal[labels+1].astype(float)
    gain=(.63+.43*np.exp(-np.maximum(steps-3,0)/32.))[:,:,None]
    q=col*gain
    if edge:
        ed=np.zeros(labels.shape,bool)
        ed[:,1:]|=labels[:,1:]!=labels[:,:-1];ed[1:,:]|=labels[1:,:]!=labels[:-1,:]
        q[ed]=q[ed]*.44+np.array([15,23,35])*.56
    return Image.fromarray(np.clip(q,0,255).astype('uint8'))
def preview():
    from selected import APPROX
    a=refined();r=roots(D,a)
    choices=[(.734,.009,.45),(.734,.009,.16),(.734,.009,.035),(.42,.03,.7),
      (.03,.0,.68),(.58,.18,.5),(.37,.38,.48),(-.27,.34,.55),
      (-.57,-.23,.7),(-.76,.15,.35),(.77,-.10,.21),(.54,.12,.25)]
    im=Image.new('RGB',(4*225,3*200),'#16202b');g=ImageDraw.Draw(im)
    for i,(cx,cy,size) in enumerate(choices):
        box=(cx-size/2,cx+size/2,cy-size*.40,cy+size*.40);l,c=plane(D,a,K,r,*box,225,180,300)
        x=(i%4)*225;y=i//4*200;im.paste(shaded(l,c),(x,y));g.text((x+4,y+182),f'{i}: {cx} {cy} width {size}',fill='white')
    im.save(RAW/'zoom-review.jpg',quality=86)

def font(size,face='regular'):
    choices={'serif':['C:/Windows/Fonts/georgia.ttf','/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf'],
             'regular':['C:/Windows/Fonts/segoeui.ttf','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'],
             'semibold':['C:/Windows/Fonts/seguisb.ttf','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf']}
    for path in choices[face]:
        try:return ImageFont.truetype(path,size)
        except OSError:continue
    return ImageFont.load_default(size=size)

def paint():
    from pathlib import Path
    bg='#f3f0e9';ink='#1d2731';muted='#475460';rule='#c7c5bb'
    image=Image.new('RGB',(1632,898),bg);pen=ImageDraw.Draw(image)
    pen.text((54,25),'Fifteen basins, one boundary',font=font(55,'serif'),fill=ink)
    pen.text((57,106),'Newton iteration for one sextic polynomial, observed every ninth step',font=font(23),fill=muted)
    pen.line((56,151,1576,151),fill=rule,width=2)
    x1,x2,yt=56,904,227
    pen.text((x1,170),'The starting plane',font=font(26,'semibold'),fill=ink)
    pen.text((x2,170),'A window around one pole',font=font(26,'semibold'),fill=ink)
    pen.text((x1,205),'Width 3.70 in the complex plane',font=font(16),fill=muted)
    pen.text((x2,205),'Width 2 × 10⁻³⁰ in the same coordinates',font=font(16),fill=muted)
    w1,h1=820,532;w2,h2=665,532
    with np.load(RAW/'case-6902-forward.npz') as z:
        l=z['labels'];s=z['steps'];first=shaded(l,s-s.min()).resize((w1,h1),Image.Resampling.LANCZOS)
    with np.load(RAW/'case-6902-microscopic.npz') as z:
        l=z['labels'];s=z['steps'];second=shaded(l,s-s.min()).resize((w2,h2),Image.Resampling.LANCZOS)
    image.paste(first,(x1,yt));image.paste(second,(x2,yt))
    cert=json.loads((OUTPUT/'certificate.json').read_text(encoding='utf-8'));px,py=cert['approximate_pole']
    lo,hi,bt,tp=BOUNDS;cx=x1+(px-lo)/(hi-lo)*w1;cy=yt+(tp-py)/(tp-bt)*h1
    for xx,yy,r in [(cx,cy,7),(x2+w2/2,yt+h2/2,.03*w2/2)]:
        pen.ellipse((xx-r,yy-r,xx+r,yy+r),outline='#17222c',width=4)
        pen.ellipse((xx-r,yy-r,xx+r,yy+r),outline='#fcf8e9',width=2)
    pen.text((x1,772),'Six roots',font=font(19,'semibold'),fill=ink)
    for i in range(D):pen.rounded_rectangle((x1+92+i*27,779,x1+113+i*27,793),2,fill=COLORS[i])
    pen.text((x1+281,772),'Nine cycle phases',font=font(19,'semibold'),fill=ink)
    for i in range(K):pen.rounded_rectangle((x1+454+i*27,779,x1+475+i*27,793),2,fill=COLORS[D+i])
    pen.text((x1,812),'Darker shades mark slower capture within each view.',font=font(17),fill=muted)
    pen.text((x2,773),'The marked pole maps to infinity.',font=font(19,'semibold'),fill=ink)
    pen.text((x2,807),'All 15 certified starting boxes lie inside the right-hand ring',font=font(17),fill=muted)
    pen.text((x2,833),'of radius 3 × 10⁻³². Every color has the same Julia boundary.',font=font(17),fill=muted)
    out=OUTPUT/'fifteen-newton-basins.png'
    # Palette quantization is presentation-only; both full label arrays are retained externally.
    image.quantize(colors=256,method=Image.Quantize.MEDIANCUT,dither=Image.Dither.NONE).save(out,optimize=True)
    print(out,out.stat().st_size,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--preview',action='store_true')
    if ap.parse_args().preview:preview()
    else:paint()
