"""Draw two certified first-image cores sharing the Alexander polynomial of T(3,4)."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from PIL import Image,ImageDraw

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from flow_figures import PAPER,INK,DIM,PALETTE,coords,fit,font,tube


def plot(rows,out):
    im=Image.new("RGB",(1460,800),PAPER)
    d=ImageDraw.Draw(im)
    d.text((40,23),"Two knots with the same Alexander polynomial",font=font(38,"bold"),fill=INK)
    d.text((43,88),"First trapped-tube cores from explicitly certified Fourier embeddings",font=font(22),fill=DIM)
    cases=((91625,"Torus knot T(3,4)","A three-strand braid makes four turns.",0),
           (152655,"Non-torus knot","The two invariants together exclude every torus knot.",1))
    for j,(n,title,caption,idx) in enumerate(cases):
        row=next(x for x in rows if x["case"]==n)
        if n==91625:
            assert row["exact_word"]==[2,1]*4
        t=np.arange(32768)*(2*np.pi/32768)
        z=sum((x+1j*y)*np.exp(1j*k*t) for k,x,y in row["harmonics"])/row["scale"]
        angle=row["degree"]*t+row["base_ripple"]*np.sin(row["degree"]*t)+row["base_offset"]
        xyz=np.column_stack(((1.15+z.real)*np.cos(angle),(1.15+z.real)*np.sin(angle),z.imag))
        w,h=680*2,442*2
        q=coords(xyz,(63,-35))
        q=fit(q,w,h)(q)
        rgb=np.empty((h,w,3),dtype=np.uint8);rgb[:]=PAPER
        depth=np.full((h,w),-np.inf)
        tube(rgb,depth,q,PALETTE[idx],w/183,True)
        left=36+j*720
        im.paste(Image.fromarray(rgb).resize((680,442),Image.Resampling.LANCZOS),(left,172))
        d.text((left+10,142),f"CASE {n}",font=font(20,"bold"),fill=DIM)
        d.text((left+10,626),title,font=font(29,"bold"),fill=INK)
        d.text((left+10,676),caption,font=font(20),fill=DIM)
    d.line((730,149,730,714),fill=(206,214,211),width=2)
    d.text((45,750),"Both cores have determinant 3. Each panel uses the same camera and preserves equal physical units.",font=font(19),fill=DIM)
    im.save(out,optimize=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificates",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    plot(json.loads(args.certificates.read_text())["rows"],args.out)
