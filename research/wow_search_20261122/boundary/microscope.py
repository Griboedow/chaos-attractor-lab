"""Microscopic forward basin image from stable near-pole coordinates and Arb spot checks."""
import argparse,json,math
import numpy as np
from numba import njit
from search import classify,RAW,OUTPUT,roots
from inverse import refined
from selected import D,K

EPS=1e-30
BOX=(-1.,1.,-.8,.8)

@njit(cache=True,nogil=True)
def micro_image(a,pole,roots,box,nx,ny):
    labels=np.zeros((ny,nx),np.int16);steps=labels.copy();c=(D-1)/D
    residue=((D-1)*pole**D-a)/(D*(D-1)*pole**(D-2))
    lo,hi,bot,top=box
    for y in range(ny):
        for x in range(nx):
            u=complex(lo+(x+.5)*(hi-lo)/nx,top-(y+.5)*(top-bot)/ny)
            # z1 = residue/(epsilon*u) + O(1), relative error O(1e-30) for the shown nonzero u.
            w=residue/(EPS*u)
            # In the external region N(w)=c*w+a*(c*w-1)/(D*w**(D-1)-a).
            # At |w|>=1e5 this absolute correction is smaller than 1e-20.
            prior=max(0,int(np.log(abs(w)/1e5)/np.log(1/c)))
            w*=c**prior
            fate,more=classify(w,D,a,K,roots,400)
            if fate>=D:fate=D+(fate-D-1-prior)%K
            labels[y,x]=fate;steps[y,x]=1+prior+more
    return labels,steps

def main(nx,ny):
    from flint import arb,acb,ctx
    from certificate import configure,nmap
    ctx.prec=656
    c,A,cyc,rb,data=configure();p0=(A/D).root(D-1);p=complex(float(p0.real.mid()),float(p0.imag.mid()))
    labels,steps=micro_image(refined(),p,roots(D,refined()),BOX,nx,ny)
    # Select independent random pixel centers and five centers from every color.
    rng=np.random.default_rng(69022026);target=list(rng.integers(0,nx*ny,224))
    for label in range(D+K):
        loc=np.flatnonzero(labels==label)
        if len(loc):target+=list(rng.choice(loc,min(5,len(loc)),replace=False))
    target=np.unique(target);pass_count=0;step_mismatch=0;reported=[]
    for idx in target:
        y,x=divmod(int(idx),nx)
        # Rational pixel-center coordinates; the y-extent is the exact interval [-4/5,4/5].
        dx=-arb(1)+arb(2*x+1)/nx;dy=arb(4)/5-arb(4*(2*y+1))/(5*ny)
        z=p0+arb('1e-30')*acb(dx,dy)
        for it in range(700):
            caught=[j for j,r in enumerate(rb) if abs(z-r)<arb('2e-7')]
            if caught:out=caught[0];break
            if abs(z)<arb('2e-7'):out=D+((-it)%K);break
            z=nmap(z,A)
        else:out=-1
        # Assignment at 2e-7 is separately rechecked one further step / period below at the certified 1e-7 capture radius.
        if out>=0:
            if out<D:
                zz=nmap(z,A);assert abs(zz-rb[out])<arb('1e-7')
            else:
                zz=z
                for j in range(K):zz=nmap(zz,A)
                assert abs(zz)<arb('1e-7')
        match=out==int(labels[y,x]);pass_count+=int(match);step_mismatch+=int(it!=steps[y,x])
        if not match:reported.append([x,y,int(labels[y,x]),out])
    path=RAW/'case-6902-microscopic.npz';np.savez_compressed(path,labels=labels,steps=steps,box=BOX,native_unit=EPS)
    report=dict(case=6902,box_in_1e_minus_30_units=list(BOX),native_unit='1e-30',shape=[ny,nx],
        normal_form='First step uses the Laurent principal part; external iterations above modulus 1e5 use z -> (5/6)z until the usual forward kernel takes over. Remaining color classification is numerical.',
        pixel_values=np.bincount(labels.ravel()+1,minlength=D+K+1).tolist(),interval_pixel_checks=len(target),
        interval_outcome_agreements=pass_count,step_disagreements=step_mismatch,disagreements=reported,pixel_seed=69022026)
    (OUTPUT/'microscopic-check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report,flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--width',type=int,default=1000);ap.add_argument('--height',type=int,default=800);args=ap.parse_args();main(args.width,args.height)
