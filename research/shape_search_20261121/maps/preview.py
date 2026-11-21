"""Select geometrically diverse candidates and build disposable review sheets."""
import argparse
import json
from pathlib import Path
from search import CACHE, START, aligned, sample, shape_features
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.spatial.distance import cdist


def load_scouts(cache,batches):
    scouts = [dict(np.load(cache/f"batch-{i:03}"/"scout.npz")) for i in batches]
    return {k:np.concatenate([s[k] for s in scouts]) for k in scouts[0]}


def selection(a,limit):
    f = a["features"]
    ok = np.isfinite(f).all(axis=1)&(f[:,1]>.0015)&(f[:,0]>.025)&(f[:,14]<.64)
    ix = np.where(ok)[0]
    # Frequency, coarse dimensions, occupancy and skewness all affect this screen.
    desc = np.column_stack((f[ix],a["rates"][ix].mean(axis=1)))
    lo,hi = np.quantile(desc,[.03,.97],axis=0)
    g = np.clip((desc-lo)/np.maximum(hi-lo,.01),-.5,1.5)
    selected = [np.argmin(np.sum((g-np.median(g,axis=0))**2,axis=1))]
    distance = np.full(len(ix),np.inf)
    for _ in range(min(limit,len(ix))-1):
        distance = np.minimum(distance,np.sum((g-g[selected[-1]])**2,axis=1))
        distance[selected] = -1.
        selected.append(int(np.argmax(distance)))
    return ix[selected]


def projection(z):
    x,_ = aligned(z)
    a,b = .35,.42
    out = np.stack((x[:,0]*np.cos(a)+x[:,2]*np.sin(a),
                    x[:,1]*np.cos(b)+(x[:,0]*np.sin(a)-x[:,2]*np.cos(a))*np.sin(b)),axis=1)
    return out


def image_grid(z,n):
    x = projection(z)
    low,high = x.min(axis=0),x.max(axis=0)
    x = (x - (high+low)/2)/max(np.max(high-low),1e-15)*.94+.5
    box = np.clip((x*n).astype(int),0,n-1)
    return np.bincount(box[:,0]+n*(n-1-box[:,1]),minlength=n*n).reshape(n,n)


def shape_screen(z):
    g = image_grid(z,96)
    filled = ndimage.binary_closing(g>0,iterations=1)
    parts,count = ndimage.label(ndimage.binary_fill_holes(filled) & ~filled)
    areas = np.bincount(parts.ravel(),minlength=count+1)[1:]
    largest = int(np.max(areas,initial=0))
    return [largest,int((areas>=12).sum()),int((g>0).sum())]


def tile(g,case,rates,metric):
    v = np.log1p(g)
    v = np.minimum(1,v/max(np.quantile(v[v>0],.95),.001))**.69
    rgb = np.stack((245*v,199*v,125*v),axis=-1).astype("uint8")
    im = Image.new("RGB",(176,197),(0,0,0))
    im.paste(Image.fromarray(rgb),(0,0))
    ImageDraw.Draw(im).text((3,177),f"{case}   {rates.mean():.3f}   H{metric[1]}",fill="white")
    return im


def main(cache,batches,limit,name):
    out=cache/name
    out.mkdir(parents=True,exist_ok=True)
    a=load_scouts(cache,batches)
    idx=selection(a,limit)
    rows=[]
    for begin in range(0,len(idx),24):
        ids=idx[begin:begin+24]
        w,tails=sample(a["params"][ids],np.tile(START,(len(ids),1)),12_000,48_000,6,36_000)
        sheet=Image.new("RGB",(176*6,197*4))
        for j,i in enumerate(ids):
            case=int(a["case"][i])
            z=tails[j]
            try:
                metrics=shape_screen(z)
                k=image_grid(z,176)
                f=shape_features(z).tolist()
            except ValueError:
                metrics=[0,0,0];k=np.ones((176,176));f=[]
            rows.append(dict(case=case,params=a["params"][i].tolist(),
                             rates=w[j].tolist(),shape=metrics,features=f))
            sheet.paste(tile(k,case,w[j],metrics),(j%6*176,j//6*197))
        sheet.save(out/f"review-{begin//24:02}.jpg",quality=91)
    (out/"preview.json").write_text(json.dumps(rows,indent=1)+"\n")
    print(json.dumps(dict(scouts=len(a["case"]),previews=len(rows),directory=str(out))),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=CACHE)
    p.add_argument("--batches",type=int,nargs="+",required=True)
    p.add_argument("--limit",type=int,default=192)
    p.add_argument("--name",default="review")
    a=p.parse_args()
    main(a.cache,a.batches,a.limit,a.name)
