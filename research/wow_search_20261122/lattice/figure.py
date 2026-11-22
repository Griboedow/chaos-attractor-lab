"""Draw the constraint surface, with measured and exact paths depth-tested against it."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
from PIL import Image, ImageDraw

from geodesic import TAU, sample, start

SIDE=860
SIZE=4.92
BACKGROUND=(9,17,28)


@njit(cache=True)
def surface(t,origin,view,n,k):
    x,y,z=origin+t*view
    return np.cos(x)+np.cos(y)+k*np.cos(x-y+n*z)


@njit(cache=True)
def draw_surface(horizontal,vertical,view,n,k):
    depth=np.full((SIDE,SIDE),-100.)
    pixels=np.empty((SIDE,SIDE,3),np.uint8)
    for iy in range(SIDE):
        for ix in range(SIDE):
            for j in range(3):pixels[iy,ix,j]=BACKGROUND[j]
            origin=((2*ix+1)/SIDE-1)*SIZE*horizontal+(1-(2*iy+1)/SIDE)*SIZE*vertical
            tmin=-9.;tmax=9.
            for axis in range(3):
                a=(-np.pi-origin[axis])/view[axis]
                b=(np.pi-origin[axis])/view[axis]
                tmin=max(tmin,min(a,b));tmax=min(tmax,max(a,b))
            if tmin>=tmax:continue
            step=.065
            phase=tmax
            previous=surface(phase,origin,view,n,k)
            while phase>tmin:
                next_phase=max(phase-step,tmin)
                current=surface(next_phase,origin,view,n,k)
                if current*previous<=0:
                    lo,hi=next_phase,phase
                    for _ in range(15):
                        mid=(lo+hi)/2
                        if surface(mid,origin,view,n,k)*previous>0:hi=mid
                        else:lo=mid
                    d=(lo+hi)/2
                    depth[iy,ix]=d
                    x,y,z=origin+d*view
                    w=x-y+n*z;u=k*np.sin(w)
                    normal=np.array([-np.sin(x)-u,-np.sin(y)+u,-n*u])
                    normal/=np.sqrt(np.sum(normal*normal))
                    if np.dot(view,normal)<0:normal=-normal
                    light=np.array([-.35,-.35,.864]);L=max(0.,np.dot(normal,light))
                    edge=(1-max(0.,np.dot(view,normal)))**3
                    t=max(0.,min(1.,.5+.5*z/np.pi))
                    diffuse=.24+.62*L+.12*max(0.,min(1.,(d+4)/8))
                    for j,base in enumerate((46+40*(1-t),130+50*t,168+41*t)):
                        pixels[iy,ix,j]=min(250,max(0,base*diffuse+85*edge))
                    break
                phase=next_phase;previous=current
    return pixels,depth


def overlay_scene(elevation=36,azimuth=26):
    elev,azim=np.deg2rad([elevation,azimuth])
    horizontal=np.array([-np.sin(azim),np.cos(azim),0.])
    vertical=np.array([-np.sin(elev)*np.cos(azim),-np.sin(elev)*np.sin(azim),np.cos(elev)])
    view=np.cross(horizontal,vertical)
    pixels,depth=draw_surface(horizontal,vertical,view,7,.375)
    im=Image.fromarray(pixels)
    ink=ImageDraw.Draw(im)
    def trace(q,color,width):
        wrap=(q+np.pi)%TAU-np.pi
        xy=np.column_stack(((wrap@horizontal/SIZE+1)*SIDE/2,(1-wrap@vertical/SIZE)*SIDE/2))
        d=wrap@view
        for j in range(1,len(xy)):
            x,y=xy[j].astype(int)
            if 0<=x<SIDE and 0<=y<SIDE and max(abs(xy[j]-xy[j-1]))<12 and d[j]+.063>=depth[y,x]:
                ink.line([tuple(xy[j-1]),tuple(xy[j])],fill=color,width=width)
    s=start(7,.375,.321)
    cloud=sample(s,7,.375,.00625,2,27000)
    trace(cloud[:,:3],(239,178,86),2)
    s=start(7,.375,0)
    trace(s[:3]+np.linspace(0,63.457061038504285,30001)[:,None]*s[3:],(234,252,255),3)
    return im


def render(packet,path):
    data=json.loads(packet.read_text())
    scene=overlay_scene()
    fig=plt.figure(figsize=(13.2,8.25),facecolor='#09111c')
    ax=fig.add_axes([.015,.035,.65,.88]);ax.imshow(scene);ax.axis('off')
    fig.text(.055,.915,'A straight line through a chaotic crystal',color='#f4f7fa',weight='bold',size=25)
    fig.text(.057,.865,'Three cosine waves form a periodic surface with fifteen handles per cell.',color='#a3b9c6',size=13)
    fig.text(.685,.72,'ONE CLOSED RAY',color='#86d6e0',weight='bold',size=11)
    fig.text(.685,.651,'(7, −7, −2)',color='#f4f7fa',size=26)
    fig.text(.685,.575,'Windings in x, y and z before\nthe straight path repeats.',color='#b2c4d1',size=13,linespacing=1.5)
    fig.text(.685,.45,'TRANSVERSE STRETCH PER LAP',color='#86d6e0',weight='bold',size=11)
    fig.text(.685,.378,'4.46 million ×',color='#f4f7fa',size=25)
    fig.text(.685,.288,'A tiny change of direction grows\nalong the unstable mode. The full\nray closes after 63.46 time units.',color='#b2c4d1',size=13,linespacing=1.45)
    fig.text(.062,.08,'━',color='#eafcff',fontsize=19)
    fig.text(.09,.083,'Exact straight orbit',color='#cbdbe2',fontsize=12)
    fig.text(.285,.08,'━',color='#efb256',fontsize=19)
    fig.text(.315,.083,'One geodesic, 337.5 time units',color='#cbdbe2',fontsize=12)
    fig.text(.685,.115,'Opposite faces are identified.\nOnly visible path segments are drawn.',color='#849bac',size=10,linespacing=1.5)
    assert data['selected']['genus']==15
    fig.savefig(path,dpi=155,facecolor=fig.get_facecolor());plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packet',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    render(args.packet,args.out)
