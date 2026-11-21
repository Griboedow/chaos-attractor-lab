"""Continue the cyclic-cloud shortlist and inspect phase-resolved spatial motion."""
import argparse
import json
from pathlib import Path
from search import CACHE, START, sample
from cyclic import phase_metrics, scan_order
from preview import projection
from PIL import Image, ImageDraw
import numpy as np
from scipy.ndimage import maximum_filter


def summarize(z,m):
    groups=z[:len(z)-len(z)%m].reshape(-1,m,3)
    lo=groups.min(axis=0);hi=groups.max(axis=0)
    diag=np.linalg.norm(np.ptp(z,axis=0))
    gaps=[]
    for i in range(m):
        for j in range(i):
            gaps.append(np.max(np.maximum(lo[i]-hi[j],lo[j]-hi[i])))
    info=phase_metrics(z,m)
    info["minimum_box_gap_fraction"]=float(min(gaps)/diag)
    info["median_component_radius_fraction"]=float(np.median(info["phase_radii_95"])/diag)
    info["maximum_component_radius_fraction"]=float(max(info["phase_radii_95"])/diag)
    return info


def draw(z,m,label):
    n=176
    x=projection(z)
    lo,hi=x.min(axis=0),x.max(axis=0)
    pix=np.clip(((x-(hi+lo)/2)/max(np.max(hi-lo),1e-15)*.93*n+n/2).astype(int),0,n-1)
    key=pix[:,0]+n*(n-1-pix[:,1])
    rgb=np.zeros((n,n,3),dtype=float)
    for i in range(m):
        den=np.bincount(key[i::m],minlength=n*n).reshape(n,n)
        v=np.log1p(den)
        v=np.minimum(1,v/max(np.quantile(v[v>0],.90),.001))**.6
        v=np.maximum(v,.29*maximum_filter(v,2))
        hue=2*np.pi*i/m
        color=150+85*np.cos(hue+np.array([0.,-2.1,2.1]))
        rgb=np.maximum(rgb,v[:,:,None]*color)
    im=Image.new("RGB",(n,197))
    im.paste(Image.fromarray(rgb.astype("uint8")),(0,0))
    ImageDraw.Draw(im).text((2,n),label,fill="white")
    return im


def inspect(cache,batches,out_name,per_order):
    out=cache/out_name
    out.mkdir(parents=True,exist_ok=True)
    hits=[]
    for b in batches:
        hits.extend(json.loads((cache/f"cycle-{b:03}"/"screen.json").read_text())["hits"])
    selected=[]
    for m in range(2,37):
        group=[r for r in hits if r["order"]==m]
        group.sort(key=lambda r:r["median_nearest_ratio"])
        # Cover broad and tiny components at every order instead of keeping
        # only chains whose component clouds appear as almost single points.
        if len(group)>per_order:
            group=[group[j] for j in np.linspace(0,len(group)-1,per_order,dtype=int)]
        selected.extend(group)
    rows=[]
    tiles=[]
    for begin in range(0,len(selected),16):
        chunk=selected[begin:begin+16]
        params=np.array([r["params"] for r in chunk])
        w,t=sample(params,np.tile(START,(len(chunk),1)),80_000,160_000,8,32_768)
        for j,r in enumerate(chunk):
            z=t[j];m=int(r["order"])
            if w[j].min()<=.007 or np.max(np.ptp(z,axis=0))<1e-8:
                continue
            m2,v=scan_order(z,36)
            if m!=m2 or v<.60:
                continue
            info=summarize(z,m)
            if info["median_nearest_ratio"]<.5:
                continue
            rows.append(dict(case=r["case"],params=r["params"],rates=w[j].tolist(),**info))
            tiles.append(draw(z,m,f"{r['case']}  m{m}  {w[j].mean():.3f}"))
    for begin in range(0,len(tiles),24):
        sheet=Image.new("RGB",(176*6,197*4))
        for j,item in enumerate(tiles[begin:begin+24]):
            sheet.paste(item,(j%6*176,j//6*197))
        sheet.save(out/f"review-{begin//24:02}.jpg",quality=92)
    (out/"preview.json").write_text(json.dumps(rows)+"\n")
    print(dict(hits=len(hits),continued=len(selected),retained=len(rows)),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--first",type=int,default=8)
    p.add_argument("--batches",type=int,default=20)
    p.add_argument("--name",default="review-cyclic")
    p.add_argument("--per-order",type=int,default=32)
    a=p.parse_args()
    inspect(a.cache,range(a.first,a.first+a.batches),a.name,a.per_order)
