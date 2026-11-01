"""A fixed three-dimensional ensemble approaching two Adamax candidates."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.spatial import cKDTree
from scipy.stats import qmc
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from adamax_map import spectra
from render_twist import camera,density,colorize


@njit(parallel=True)
def evolve(starts,p,times):
    snapshots=np.empty((len(times),len(starts),3));tails=np.empty((len(starts),16,3))
    for lane in prange(len(starts)):
        x,m,u=starts[lane];b,c,e=p;cursor=0
        for t in range(times[-1]+1):
            if t==times[cursor]:
                snapshots[cursor,lane]=np.array([x,m,u]);cursor=min(cursor+1,len(times)-1)
            if t>times[-1]-16:tails[lane,t-(times[-1]-16)-1]=np.array([x,m,u])
            m=b*m+(1-b)*x;u=max(c*u,abs(x)+e);x-=m/u
    return snapshots,tails


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    source=np.load(args.source);render=json.loads(args.render.read_text())['runs'];ids=[847,4099]
    times=np.array([0,1,4,16,128,2048]);unit=qmc.Sobol(3,scramble=True,seed=2026110125).random_base2(15)
    boxes=[np.array([[-.7,-.2,.4],[.7,.2,.8]]),np.array([[-1.6,-.4,.8],[1.6,.4,1.6]])]
    protocol=dict(ids=ids,count_per_case=len(unit),seed=2026110125,box_bounds=[b.tolist() for b in boxes],
        selection='Boxes chosen from rounded candidate extents, with u bounded away from zero. A genuine 3D initial volume, not a cold-start line.',
        times=times.tolist(),map='Autonomous normalized Adamax limit.',
        expansion_cohort='First 128 Sobol starts from each box, declared before simulation; 100000 autonomous warmup, 200000 observation in four 50000 windows.',
        support='Directed endpoint Euclidean distance to the 2000000 preview states, divided by reference extent norm.',
        source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),render_sha256=hashlib.sha256(args.render.read_bytes()).hexdigest(),
        scope='Finite numerical approach from a declared volume; not a rigorous trapping region or basin measure.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for name in ['approach_adamax.py','adamax_map.py','render_twist.py']:(out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    canvas=Image.new('RGB',(2760,1180),'#071320');draw=ImageDraw.Draw(canvas)
    large=ImageFont.truetype(findfont('DejaVu Sans'),43);small=ImageFont.truetype(findfont('DejaVu Sans'),25)
    draw.text((35,20),'FROM A VOLUME TO FOLDED SHEETS',font=large,fill='#e6edf5')
    draw.text((35,80),'32,768 initial states per row / same points followed through time / no interpolation',font=small,fill='#9eb5c9')
    rows=[];all_parameters=[];all_starts=[]
    for row_index,(ident,box) in enumerate(zip(ids,boxes)):
        k=list(source['ids']).index(ident);p=source['parameters'][k];starts=box[0]+unit*(box[1]-box[0])
        snapshots,tails=evolve(starts,p,times);assert np.isfinite(snapshots).all()
        row=next(r for r in render if r['id']==ident);center=np.array(row['center']);basis=np.array(row['basis']);projection=camera(row['azimuth'],.25)
        q=(snapshots-center)@basis;halfwidth=float(np.max(np.linalg.norm(q,axis=2))*1.045)
        top=200+row_index*460
        draw.text((35,top-48),row['title'],font=small,fill='#e6edf5')
        for j,t in enumerate(times):
            h=density(q[j],projection,430,halfwidth);assert int(h.sum())==len(starts)
            im,_=colorize(h,vmax=np.log1p(70));left=25+j*455;canvas.paste(im,(left,top))
            draw.text((left+12,top+12),f'n = {t:,}',font=small,fill='#cbd9e5')
        preview=np.load(args.preview/f'case-{ident}.npy',mmap_mode='r');tree=cKDTree(preview)
        scale=np.linalg.norm(np.ptp(preview,axis=0));dist=tree.query(snapshots[-1])[0]/scale
        periods=np.zeros(len(starts),dtype=np.int16)
        for lag in range(1,9):
            residual=np.max(np.linalg.norm(tails[:,lag:]-tails[:,:-lag],axis=2),axis=1)
            hit=(residual<1e-10*(1+np.linalg.norm(np.ptp(tails,axis=1),axis=1)))&(periods==0);periods[hit]=lag
        np.savez_compressed(out/f'case-{ident}.npz',parameters=p,starts=starts,times=times,snapshots=snapshots,tails=tails,periods=periods,endpoint_distances=dist)
        rows.append(dict(id=ident,endpoint_distances_q50_q95_q99_max=np.quantile(dist,[.5,.95,.99,1]).tolist(),
            tail_period_up_to8_counts={str(k):int((periods==k).sum()) for k in np.unique(periods)},
            halfwidth=halfwidth,projection=projection.tolist(),center=center.tolist(),basis=basis.tolist(),spatially_clipped=0))
        all_parameters.extend([p.copy() for _ in range(128)]);all_starts.extend(starts[:128]);print('Approach',ident,rows[-1]['endpoint_distances_q50_q95_q99_max'],flush=True)
    draw.text((35,1130),'Each row uses one fixed spatial scale at every time. Final clusters remain discrete iterates, not drawn connecting curves.',font=small,fill='#9eb5c9')
    canvas.save(out/'adamax-approach.png')
    windows,tails,ends,resets,margins=spectra(np.array(all_parameters),np.array(all_starts),100000,200000,50000,4096,False)
    assert np.isfinite(windows).all() and np.isfinite(tails).all()
    np.savez_compressed(out/'volume-checks.npz',parameters=all_parameters,starts=all_starts,windows=windows,tails=tails,ends=ends,reset_counts=resets,minimum_margins=margins)
    for i,row in enumerate(rows):
        w=windows[128*i:128*(i+1)]
        row.update(volume_checks=128,all_four_windows_above005=int(np.all(w[:,:,0]>.005,axis=1).sum()),minimum_leading_window=float(w[:,:,0].min()))
        print('Volume expansion',row['id'],row['all_four_windows_above005'],'/128',row['minimum_leading_window'],flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--render',type=Path,required=True);p.add_argument('--preview',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
