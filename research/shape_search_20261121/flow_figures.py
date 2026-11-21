"""Draw labeled components and a periodic skeleton from saved native spatial states."""
import argparse
from pathlib import Path
import sys

import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from numba import njit

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
from render_portraits import camera

PAPER=(249,248,244)
INK="#172e38"
DIM="#47616c"
PALETTE=np.array([[3,130,150],[211,94,51],[125,88,183],[93,145,66]],dtype=float)


def font(size,weight="regular"):
    from matplotlib.font_manager import FontProperties
    return ImageFont.truetype(findfont(FontProperties(family="DejaVu Sans",weight=weight)),size)


def coords(z,view):
    p=camera(*view)
    return z@np.vstack([p,np.cross(*p)]).T


@njit(cache=True)
def tube(rgb,depth,points,color,r,close):
    n=len(points)
    for i in range(n if close else n-1):
        a,b=points[i],points[(i+1)%n]
        steps=int(np.ceil(max(abs(a[0]-b[0]),abs(a[1]-b[1]))/.85))+1
        for step in range(steps):
            t=step/steps
            u=(1-t)*a+t*b
            for iy in range(max(0,int(u[1]-r)),min(len(rgb),int(u[1]+r)+1)):
                for ix in range(max(0,int(u[0]-r)),min(rgb.shape[1],int(u[0]+r)+1)):
                    dx,dy=(ix-u[0])/r,(iy-u[1])/r
                    mag=1-dx*dx-dy*dy
                    if mag<=0:
                        continue
                    nz=np.sqrt(mag)
                    z=u[2]+r*nz
                    if z<depth[iy,ix]:
                        continue
                    depth[iy,ix]=z
                    light=.78+.29*(.64*nz-.45*dx-.45*dy)
                    for k in range(3):
                        rgb[iy,ix,k]=min(255,max(0,color[k]*light))


def fit(points,w,h):
    low=np.min(points,axis=0)
    high=np.max(points,axis=0)
    scale=.90*min(w/(high[0]-low[0]),h/(high[1]-low[1]))
    def project(p):
        p=(p-(low+high)/2)*scale
        p[:,0]+=w/2
        p[:,1]=h/2-p[:,1]
        return p
    return project


def read(path):
    with np.load(path) as data:
        return data["points"]


def link_tile(cache,number,view,w,h):
    path=cache/"flows/symmetry-links-a"/f"case-{number}.npz"
    with np.load(path) as f:
        z,offset=f["points"],f["offsets"]
    xyz=coords(z,view)
    project=fit(xyz,w,h)
    xyz=project(xyz)
    rgb=np.empty((h,w,3),dtype=np.uint8);rgb[:]=PAPER
    depth=np.full((h,w),-np.inf)
    for j,(a,b) in enumerate(zip(offset,offset[1:])):
        tube(rgb,depth,xyz[a:b:4],PALETTE[j],w/285,True)
    return Image.fromarray(rgb)


def links(cache,out):
    panels=((2479,(35.26438968,45),"Three interlocked orbits","Every pair has |Lk| = 30."),
            (26471,(76,-15),"Four interlocked orbits","All six pairs have |Lk| = 6."))
    im=Image.new("RGB",(1536,842),PAPER)
    draw=ImageDraw.Draw(im)
    draw.text((44,26),"Coexisting periodic orbits",font=font(46,"bold"),fill=INK)
    draw.text((46,93),"Each color traces a separate attracting cycle of the same flow.",font=font(24),fill=DIM)
    for j,(case,view,title,caption) in enumerate(panels):
        left=36+758*j
        im.paste(link_tile(cache,case,view,714*2,510*2).resize((714,510),Image.Resampling.LANCZOS),(left,169))
        draw.text((left+12,143),f"CASE {case}",font=font(20,"bold"),fill=DIM)
        draw.text((left+12,689),title,font=font(30,"bold"),fill=INK)
        draw.text((left+13,737),caption,font=font(23),fill=DIM)
    draw.line((767,152,767,781),fill=(206,214,211),width=2)
    draw.text((47,800),"Native spatial coordinates. Linking numbers are computed on closed polygonal orbits.",font=font(19),fill=DIM)
    im.save(out,optimize=True)


def overlay_tile(cloud,cycle,view,w,h):
    q=coords(cloud,view)
    v=coords(cycle,view)
    project=fit(np.vstack([q,v]),w,h)
    z=project(q)
    counts=np.histogram2d(z[:,1],z[:,0],bins=(h,w),range=((0,h),(0,w)))[0]
    density=np.log1p(counts)
    limit=np.quantile(density[counts>0],.997)
    alpha=np.minimum(1,density/limit)**.55
    bg=np.array(PAPER)
    color=np.array([9,87,110])
    rgb=np.rint(bg*(1-alpha[:,:,None])+alpha[:,:,None]*color).astype(np.uint8)
    im=Image.fromarray(rgb)
    draw=ImageDraw.Draw(im)
    points=[tuple(v[:2]) for v in project(v)[::3]]
    draw.line(points+[points[0]],fill=(203,61,31),width=2,joint="curve")
    return im


def ribbon(cache,out):
    cloud=read(cache/"flows/dense-ring-knot/case-180844.npz")
    cycle=read(cache/"flows/knot-continuation-a/case-180844.npz")
    im=Image.new("RGB",(1536,810),PAPER)
    d=ImageDraw.Draw(im)
    d.text((43,22),"A chaotic band around a periodic knot",font=font(45,"bold"),fill=INK)
    d.text((46,90),"Case 180844  ·  The same trajectories seen from two directions",font=font(24),fill=DIM)
    for j,(label,view) in enumerate((("Oblique view",(27,-45)),("Near-axial view",(70,-15)))):
        left=34+j*760
        tile=overlay_tile(cloud,cycle,view,726,505)
        im.paste(tile,(left,174))
        d.text((left+15,143),label,font=font(23,"bold"),fill=DIM)
    d.line((768,146,768,670),fill=(206,214,211),width=2)
    d.line((49,716,98,716),fill=(9,87,110),width=7)
    d.text((114,697),"Sampled attracting band",font=font(24),fill=INK)
    d.line((518,716,567,716),fill=(203,61,31),width=4)
    d.text((584,697),"Unstable periodic orbit",font=font(24),fill=INK)
    d.text((47,759),"The trajectories share physical coordinates and scale within each view. No registration was fitted.",font=font(19),fill=DIM)
    im.save(out,optimize=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("kind",choices=["links","ribbon"])
    ap.add_argument("--cache",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    (links if args.kind=="links" else ribbon)(args.cache,args.out)
