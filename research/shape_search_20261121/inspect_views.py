"""Render equal-unit orthographic thumbnails of saved clouds for visual selection."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
from render_portraits import COLOR,camera

VIEWS = [(27,-45),(35.26438968,45),(75,-15),(5,0),(15,90),(10,145),(46,-124)]
COLORS = np.round(255*COLOR(np.linspace(0,1,256))[:,:3]).astype("uint8")


def tile(points, projection, size):
    pos = points @ projection.T
    pos -= (pos.max(0)+pos.min(0))/2
    r = np.abs(pos).max()*1.035
    h = np.histogram2d(pos[:,0],pos[:,1],bins=size,range=[[-r,r],[-r,r]])[0].T[::-1]
    assert int(h.sum()) == len(points)
    density = np.log1p(h)
    sat = np.quantile(density[h>0],.998)
    return Image.fromarray(COLORS[np.round(255*np.minimum(1,density/sat)**.72).astype("uint8")])


def render(inputs,cache,out,all_angles):
    side,caption = 204,26
    font = ImageFont.truetype(findfont("DejaVu Sans"),13)
    cols = len(VIEWS) if all_angles else 6
    rows = 6 if all_angles else 5
    per_page = rows if all_angles else cols*rows
    for start in range(0,len(inputs),per_page):
        page_items = list(inputs.items())[start:start+per_page]
        sheet = Image.new("RGB",(side*cols,(side+caption)*(len(page_items) if all_angles else (len(page_items)+cols-1)//cols)), "white")
        draw = ImageDraw.Draw(sheet)
        for j,(label,path) in enumerate(page_items):
            with np.load(cache/path) as f:
                x = f["points"]
            if len(x)>100000:
                choose = np.random.default_rng(1713).choice(len(x),100000,replace=False)
                x=x[choose]
            for k,(e,a) in enumerate(VIEWS if all_angles else VIEWS[:1]):
                col,row = (k,j) if all_angles else (j%cols,j//cols)
                left,top = col*side,row*(side+caption)
                sheet.paste(tile(x,camera(e,a),side),(left,top+caption))
                detail = f"   {k}" if all_angles else ""
                draw.text((left+5,top+4),label+detail,font=font,fill="#243d49")
        path = out/f"sheet-{start//per_page:02d}.png"
        sheet.save(path)


if __name__ == "__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs",type=Path,required=True)
    ap.add_argument("--cache",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--views",action="store_true")
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    render(json.loads(args.inputs.read_text(encoding="utf-8")),args.cache,args.out,args.views)
