"""Build a small external contact sheet used to select a certified parameter."""
import colorsys,json
import numpy as np
from PIL import Image,ImageDraw
from search import RAW,roots,plane

BASE=['#daab61','#6ca0a7','#9885bd','#be7b7c','#87a77e','#91a8ca']
PHASE=['#55ab92','#669d7b','#91bb7d','#acc66a','#d6d584','#ddab80','#c47774','#9c6d9d','#6f82ab']
def shade(labels,counts,d,k):
    def rgb(hex):return np.array([int(hex[i:i+2],16) for i in (1,3,5)])
    col=np.array([[36,32,37]]+[rgb(v) for v in BASE[:d]]+[rgb(PHASE[round(t*8/max(1,k-1))]) for t in range(k)])
    y=col[labels+1]
    brightness=(.39+.68*np.exp(-np.asarray(counts)/24))[:,:,None]
    return np.clip(y*brightness,0,255).astype(np.uint8)

def main():
    rows=json.loads((RAW/'candidate-search.json').read_text())['rows']
    subset=[]
    for d in (4,5,6):
        for k in (3,5,7,9):
            best=sorted([r for r in rows if r['d']==d and r['k']==k],key=lambda r:r['score'],reverse=True)
            if best:subset+=best[:2]
    im=Image.new('RGB',(4*240,((len(subset)+3)//4)*206),'#182027');draw=ImageDraw.Draw(im)
    for j,row in enumerate(subset):
        d,k=row['d'],row['k'];a=complex(*row['a']);ext=row['extent'];h=ext*.78
        l,n=plane(d,a,k,roots(d,a),-ext,ext,-h,h,240,180,160)
        x=(j%4)*240;y=(j//4)*206;im.paste(Image.fromarray(shade(l,n,d,k)),(x,y))
        draw.text((x+2,y+184),f"{j}: {d}+{k}, {a.real:.3f}+{a.imag:.3f}i, {row['cycle_fraction']:.2f}",fill='white')
    im.save(RAW/'contact.jpg',quality=84);(RAW/'contact-cases.json').write_text(json.dumps(subset,indent=2)+'\n')
if __name__=='__main__':main()
