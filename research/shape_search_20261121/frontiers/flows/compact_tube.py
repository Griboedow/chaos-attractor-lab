"""Test a compact symmetry-reduced return tube, including independent transverse offsets."""
import json
import numpy as np
from reduced import reduced
from native_section import RAW,RESULT,packed
from return_tube import CENTER,graph
ZL=.98071
ZH=.99943
W=1.e-6

def main():
    g=graph();p=packed()
    z,d=np.meshgrid(np.linspace(ZL,ZH,1537),np.linspace(-W,W,17),indexing='ij')
    x=g(z-CENTER)+d;a=reduced(x.ravel(),z.ravel(),.00625,p).reshape(*z.shape,4)
    assert np.all(a[:,:,3]==1)
    zp=a[:,:,2];dp=a[:,:,1]-g(zp-CENTER)
    q=dict(z_domain=[ZL,ZH],delta_halfwidth=W,grid=list(z.shape),count=int(z.size),
         z_image=[zp.min(),zp.max()],delta_image=[dp.min(),dp.max()],
         return_time=[a[:,:,0].min(),a[:,:,0].max()],z_margin=min(zp.min()-ZL,ZH-zp.max()),delta_margin=W-np.max(abs(dp)))
    np.savez_compressed(RAW/'reduced-small.npz',z=z,d=d,x=x,returns=a)
    print(q,flush=True);(RESULT/'reduced-tube.json').write_text(json.dumps(q,indent=2)+'\n')

if __name__=='__main__':main()
