"""Replay selected parameter sets in binary64 and draw private inspection sheets."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw
from model import replay


def picture(points,width=245,height=178,view=0):
    z=points-points.mean(axis=0)
    _,v=np.linalg.eigh(z.T@z)
    z=z@v[:,::-1]
    for j in range(3):
        if np.mean(z[:,j]**3)<0: z[:,j]*=-1
    if view==0:
        coords=np.column_stack((.85*z[:,0]+.50*z[:,2],-.25*z[:,0]+.87*z[:,1]+.44*z[:,2]))
    elif view==1:
        coords=z[:,:2]
    else:
        coords=z[:,1:]
    lo,hi=np.quantile(coords,[.0005,.9995],axis=0)
    coords-=.5*(lo+hi)
    scale=.91*min(width/max(hi[0]-lo[0],1e-15),height/max(hi[1]-lo[1],1e-15))
    pixels=(coords*scale+np.array([width,height])/2).astype(int)
    good=(pixels[:,0]>=0)&(pixels[:,0]<width)&(pixels[:,1]>=0)&(pixels[:,1]<height)
    pixels=pixels[good]
    count=np.bincount(pixels[:,1]*width+pixels[:,0],minlength=height*width).reshape(height,width)
    image=np.log1p(count)**.65
    image=np.minimum(image/max(np.quantile(image[image>0],.985),1e-6),1)
    bg=np.array([8,18,26])
    fg=np.array([121,233,228])
    return Image.fromarray((bg+(fg-bg)*image[:,:,None]).astype('uint8')[::-1])


def make_sheets(paths,out,view):
    for start in range(0,len(paths),24):
        page=Image.new('RGB',(245*4,201*6),(8,18,26))
        draw=ImageDraw.Draw(page)
        for k,path in enumerate(paths[start:start+24]):
            with np.load(path) as f:
                img=picture(f['points'],view=view)
                j=int(f['case'])
                rates=f['rates'][:,0]
            x=k%4*245; y=k//4*201
            page.paste(img,(x,y+18))
            draw.text((x+6,y+3),f'{start+k:03d}  case {j}   {rates.mean():+.3f}',fill='white')
        page.save(out/f'sheet-v{view}-{start//24:02d}.jpg',quality=89)


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    with np.load(args.selection) as f:
        ids=f['ids']; data=f['data'].astype(np.float64); kinds=f['kinds']
    paths=[]
    for i,(case,raw,k) in enumerate(zip(ids,data,kinds)):
        path=args.out/f'case-{case}.npz'
        if not path.exists():
            cloud,rates,minimum,maxima,branches=replay(raw,raw[1,:3].copy(),args.burn,args.window,4,args.tail,False)
            tmp=path.with_suffix('.tmp.npz')
            np.savez(tmp,points=cloud.astype('f4'),rates=rates,minima=minimum,max_norm=maxima,case=case,data=raw,kind=k,branches=branches)
            os.replace(tmp,path)
        paths.append(path)
        if i%24==23 or i==len(ids)-1:
            print(i+1,flush=True)
    for v in args.views: make_sheets(paths,args.out,v)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--selection',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--burn',type=int,default=16384)
    p.add_argument('--window',type=int,default=8192)
    p.add_argument('--tail',type=int,default=24576)
    p.add_argument('--views',nargs='+',type=int,default=[0])
    main(p.parse_args())
