"""Render every retained native-coordinate point cloud into compact study plates."""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont


def font(size):
    for name in ('DejaVuSans.ttf','arial.ttf'):
        try:return ImageFont.truetype(name,size)
        except OSError:pass
    return ImageFont.load_default(size=size)


def density(cloud,width,height):
    z=cloud.astype(float);z-=z.mean(0)
    _,v=np.linalg.eigh(z.T@z);v=v[:,::-1]
    signed=z@v;signed*=np.where(np.mean(signed**3,axis=0)>=0,1.,-1.)
    right=np.array([np.sqrt(3)/2,0,.5]);up=np.array([-.25,np.sqrt(3)/2,np.sqrt(3)/4])
    xy=signed@np.column_stack((right,up))
    low,high=np.quantile(xy,[.0005,.9995],axis=0)
    xy-=(low+high)*.5
    scale=.92*min(width/(high[0]-low[0]),height/(high[1]-low[1]))
    pixels=np.floor(xy*scale+np.array([width,height])/2).astype(int)
    valid=((pixels>=0)&(pixels<np.array([width,height]))).all(1)
    pix=pixels[valid]
    counts=np.bincount(pix[:,1]*width+pix[:,0],minlength=height*width).reshape(height,width)[::-1]
    intensity=np.log1p(counts)**.72
    intensity=np.minimum(1,intensity/np.quantile(intensity[counts>0],.992))
    dark=np.array([12,23,30]);light=np.array([128,230,215])
    return Image.fromarray((dark[None,None,:]+(light-dark)[None,None,:]*intensity[:,:,None]).astype('u1'))


def main(args):
    rows=json.loads((args.source/'manifest.json').read_text())
    groups=[('Single ring torus','single',[r for r in rows if r['study']=='base survey' and r['family'].startswith('Ring')]),
            ('Nearest linked torus','linked',[r for r in rows if r['study']=='base survey' and r['family'].startswith('Linked')]),
            ('Blended linked tori','blended',[r for r in rows if r['study']=='base survey' and r['family'].startswith('Blended')]),
            ('Three positive numerical exponents','three',[r for r in rows if r['study']=='three-direction growth'])]
    args.out.mkdir(parents=True,exist_ok=True)
    index=[f'# Portraits of {len(rows)} retained parameter sets','','Each panel uses the last 250,000 native states of its own orbit. The camera is fixed in signed principal-axis coordinates. All axes share one scale within each panel. Density uses logarithmic compression. Labels identify a parameter set and a visual description.']
    total=0;plates=[]
    for title,slug,group in groups:
        index+=['',f'## {title}','']
        for j in range(0,len(group),12):
            items=group[j:j+12];columns=3;tilew,tileh=325,210
            height=((len(items)+columns-1)//columns)*tileh+52
            page=Image.new('RGB',(columns*tilew,height),(12,23,30));draw=ImageDraw.Draw(page)
            draw.text((15,12),f'{title}  ·  {j+1}–{j+len(items)} / {len(group)}',font=font(23),fill=(234,240,233))
            for n,row in enumerate(items):
                with np.load(row['cloud']) as f:img=density(f['points'],tilew-10,tileh-40)
                x=n%columns*tilew+5;y=n//columns*tileh+50
                page.paste(img,(x,y));draw.text((x+8,y+tileh-37),f'{row["case"]}   {row["label"]}',font=font(13),fill=(214,224,219))
            name=f'{slug}-{j//12+1}.webp'
            page.save(args.out/name,'WEBP',quality=84,method=6)
            if (args.out/name).stat().st_size>=262144:raise ValueError('Plate exceeds the repository asset limit')
            plates.append(dict(file=name,cache_relative='projection/atlas/'+name,sha256=hashlib.sha256((args.out/name).read_bytes()).hexdigest(),cases=[row['case'] for row in items]))
            index+=['',f'![{title}, cases {j+1} through {j+len(items)}]({name})']
            total+=len(items)
    assert total==len(rows)
    (args.out/'README.md').write_text('\n'.join(index).rstrip()+'\n',encoding='utf8')
    if args.manifest:
        info=dict(count=total,format='WebP quality 84',camera='Signed descending PCA, equal physical axis scale, right=(sqrt(3)/2,0,1/2), up=(-1/4,sqrt(3)/2,sqrt(3)/4). The 0.05 and 99.95 percentiles set the bounds.',plates=plates)
        args.manifest.write_text(json.dumps(info,indent=2)+'\n',encoding='utf8')
    print('rendered',total,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--manifest',type=Path)
    main(p.parse_args())
