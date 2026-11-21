"""Make a compact cross-family figure from native samples and labeled cycle components."""
import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image,ImageDraw

from flow_figures import font,link_tile
from inspect_views import tile
from render_portraits import camera

ROOT=Path(__file__).resolve().parents[2]
SELECTION=[
    ("flows",132,"Six-loop rosette",(35.26438968,45)),
    ("links",26471,"Four linked cycles",(76,-15)),
    ("flows",180844,"Band around a knot",(27,-45)),
    ("flows",21253,"Lobed oscillator",(75,-15)),
    ("projection",301027,"Nested scrolls",(75,-15)),
    ("projection",370279,"Ribbed fans",(75,-15)),
    ("projection",2189166,"Helical chain",(15,90)),
    ("maps",10803912,"Twisted mesh",(10,145)),
    ("delay",1008712,"Folded delay ring",(27,-45)),
    ("braids",114425,"Knotted solenoid",(75,-15)),
    ("pinched",156,"Pinched quasiperiodic graph",(27,-45)),
    ("flows",306,"Six-sided spiral lattice",(35.26438968,45)),
]


def draw(cache,out):
    paths=json.loads((ROOT/"results/shape_search_20261121/comparison-inputs.json").read_text())
    side,margin=310,34
    im=Image.new("RGB",(4*(side+margin)+margin,3*(side+78)+128),"white")
    d=ImageDraw.Draw(im)
    d.text((margin,17),"Twelve spatial attractor portraits",font=font(40,"bold"),fill="#152e39")
    d.text((margin,72),"Full-state samples from the search of 2026-11-21",font=font(21),fill="#49616d")
    for i,(track,number,label,view) in enumerate(SELECTION):
        x=margin+(i%4)*(side+margin)
        y=125+(i//4)*(side+78)
        if track=="links":
            image=link_tile(cache,number,view,side*2,side*2).resize((side,side),Image.Resampling.LANCZOS)
        else:
            with np.load(cache/paths[f"{track}/{number}"]) as f:
                points=f["points"]
            image=tile(points,camera(*view),side)
        im.paste(image,(x,y))
        d.text((x+4,y+side+8),label,font=font(20,"bold"),fill="#152e39")
        d.text((x+5,y+side+36),f"Case {number}",font=font(17),fill="#49616d")
    im.save(out,optimize=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    draw(args.cache,args.out)
