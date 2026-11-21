"""Render the positive x meridian section of the five Case 180844 orbit polygons."""
import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[3]
COLORS = {1:'#7E4698',2:'#BF4928',4:'#257449',8:'#2062A1',16:'#A9770A'}
MARKERS = {1:'circle',2:'square',4:'triangle',8:'diamond',16:'cross'}


def font(size, bold=False):
    for family in (('arialbd.ttf','DejaVuSans-Bold.ttf') if bold else ('arial.ttf','DejaVuSans.ttf')):
        try:
            return ImageFont.truetype(family,size)
        except OSError:
            continue
    return ImageFont.load_default(size=size)


def section(path, closed):
    with np.load(path) as saved:
        points=saved['points'].astype(float)
    end=np.roll(points,-1,axis=0) if closed else points[1:]
    start=points if closed else points[:-1]
    keep=(start[:,1]<=0)&(end[:,1]>0)
    a,b=start[keep],end[keep]
    hits=a+(b-a)*(-a[:,1]/(b[:,1]-a[:,1]))[:,None]
    return hits[hits[:,0]>0][:,[0,2]]


def mark(draw,x,y,order,size):
    color=COLORS[order];k=MARKERS[order];w=3
    if k=='circle':
        draw.ellipse((x-size,y-size,x+size,y+size),outline=color,width=w)
    elif k=='square':
        draw.rectangle((x-size,y-size,x+size,y+size),outline=color,width=w)
    elif k=='triangle':
        draw.line([(x,y-size),(x+size,y+size),(x-size,y+size),(x,y-size)],fill=color,width=w)
    elif k=='diamond':
        draw.line([(x-size,y),(x,y-size),(x+size,y),(x,y+size),(x-size,y)],fill=color,width=w)
    else:
        draw.line([(x-size,y),(x+size,y)],fill=color,width=w)
        draw.line([(x,y-size),(x,y+size)],fill=color,width=w)


def panel(draw,box,xr,yr,xt,yt,precision,title,cloud,hits,detail=None):
    left,top,right,bottom=box
    def xy(x,y):
        return left+(x-xr[0])/(xr[1]-xr[0])*(right-left),bottom-(y-yr[0])/(yr[1]-yr[0])*(bottom-top)
    fg='#253342';muted='#505D68';grid='#E3E9EC'
    draw.text((left,top-42),title,font=font(22,True),fill=fg)
    for n in xt:
        x,_=xy(n,yr[0]);draw.line((x,top,x,bottom),fill=grid)
        draw.text((x,bottom+8),f'{n:.{precision}f}',font=font(15),fill=muted,anchor='mt')
    for n in yt:
        _,y=xy(xr[0],n);draw.line((left,y,right,y),fill=grid)
        draw.text((left-8,y),f'{n:.{precision}f}',font=font(15),fill=muted,anchor='rm')
    draw.line([(left,top),(left,bottom),(right,bottom)],fill=muted,width=2)
    draw.text(((left+right)/2,bottom+33),'Radius r',font=font(18),fill=fg,anchor='mt')
    draw.text((left-17,top-12),'z',font=font(18),fill=fg,anchor='rm')
    if cloud is not None:
        for a,b in cloud:
            x,y=xy(a,b);draw.ellipse((x-2,y-2,x+2,y+2),fill='#9CA8AF')
    for order in (1,2,4,8,16):
        for a,b in hits[order]:
            if xr[0]<=a<=xr[1] and yr[0]<=b<=yr[1]:
                x,y=xy(a,b);mark(draw,x,y,order,{1:10,2:9,4:8,8:7,16:5}[order] if cloud is not None else {1:13,2:11,4:10,8:9,16:6}[order])
    if detail:
        x,y=xy(detail[0][1],detail[1][0]);x0,y0=xy(detail[0][0],detail[1][1])
        draw.line([(x0-3,y0-3),(x+4,y0-3),(x+4,y+4),(x0-3,y+4),(x0-3,y0-3)],fill='#26313B',width=1)
        draw.line([(x+5,y0-1),(x+25,y0-19)],fill='#26313B',width=1)
        draw.text((x+30,y0-30),'detail',font=font(15),fill=fg)


def main(args):
    records=ROOT/'results/shape_search_20261121/flows'
    rows=json.loads((records/'section-crosscheck.json').read_text())['rows']
    hits={r['period_order']:section(args.cache_root/r['cache_relative'],True) for r in rows}
    for row in rows:
        assert len(hits[row['period_order']])==row['distinct_primary_section_points']
    band=section(args.cache_root/'flows/dense-ring-knot/case-180844.npz',False)
    im=Image.new('RGB',(1440,715),'white');draw=ImageDraw.Draw(im)
    draw.text((59,28),'Meridian returns of five periodic orbits',font=font(37,True),fill='#24323E')
    draw.text((61,79),'Case 180844   ·   y = 0, x > 0   ·   native radius and height',font=font(21),fill='#53616D')
    for j,order in enumerate((1,2,4,8,16)):
        x=72+j*199;mark(draw,x,134,order,7)
        draw.text((x+16,134),f'Order {order}  ({len(hits[order])} hits)',font=font(17),fill='#34424E',anchor='lm')
    draw.ellipse((1093,130,1100,137),fill='#9CA8AF')
    draw.text((1108,134),f'Band  ({len(band)} hits)',font=font(17),fill='#52606B',anchor='lm')
    zoom=((.8609,.8621),(.92200,.92280))
    panel(draw,(96,230,686,565),(.854,1.002),(.851,.931),
          [.86,.89,.92,.95,.98,1.00],[.86,.88,.90,.92],2,'Full section',band,hits,zoom)
    panel(draw,(823,230,1380,565),*zoom,[.8610,.8613,.8616,.8619],
          [.9220,.9222,.9224,.9226,.9228],4,'Upper branch detail · periodic orbits',None,hits)
    draw.text((60,650),'Each symbol is a linearly interpolated crossing. Orbit polygons close once. The sampled band is open and consecutive hits are unconnected.',font=font(17),fill='#4E5B66')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    im.save(args.out,optimize=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--cache-root',type=Path,required=True)
    p.add_argument('--out',type=Path,default=ROOT/'docs/assets/shape-search-20261121-meridian.png')
    main(p.parse_args())
