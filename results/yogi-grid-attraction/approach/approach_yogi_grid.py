"""Follow declared three-dimensional ensembles toward two Yogi grid regimes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from scipy.stats import qmc
from scipy.spatial import cKDTree
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from yogi_map import spectra
from render_twist import camera,density,colorize


@njit(parallel=True)
def evolve(p,starts,times):
    snapshots=np.empty((len(times),len(starts),3));tails=np.empty((len(starts),64,3))
    b,c,e=p
    for i in prange(len(starts)):
        x,m,v=starts[i];cursor=0
        for t in range(times[-1]+1):
            if t==times[cursor]:
                snapshots[cursor,i]=(x,m,v);cursor=min(cursor+1,len(times)-1)
            if t>times[-1]-64:tails[i,t-(times[-1]-64)-1]=(x,m,v)
            if t<times[-1]:
                s=1 if v>x*x else (-1 if v<x*x else 0)
                m=b*m+(1-b)*x;v=v-(1-c)*s*x*x;x=x-m/(np.sqrt(v)+e)
    return snapshots,tails


def main(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);set_num_threads(4)
    d=np.load(args.source);preview=json.loads((args.preview/'analysis.json').read_text())['runs']
    ids=[670,125];times=np.array([0,1,16,256,8192]);unit=qmc.Sobol(3,scramble=True,seed=2026110151).random_base2(15)
    boxes=[];parameters=[]
    for ident in ids:
        k=list(d['ids']).index(ident);cloud=d['tails'][k];extent=np.ptp(cloud,axis=0)
        lo=cloud.min(0)-.1*extent;hi=cloud.max(0)+.1*extent;lo[2]=max(lo[2],1e-8)
        boxes.append([lo,hi]);parameters.append(d['parameters'][k])
    protocol=dict(grid_ids=ids,source=str(args.source),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        count_per_case=len(unit),seed=2026110151,times=times.tolist(),boxes=np.array(boxes).tolist(),
        selection='Full3D five-million-tail bounding box expanded10% in each direction, v lower bound clipped to1e-8. Same Sobol unit triples in both cases.',
        implementation='Independent scalar autonomous Yogi recurrence; no interpolation, artificial noise or connecting lines.',
        expansion_cohort='First128 Sobol starts per box, declared before simulation;100000 warmup,200000 observation,4x50000 windows,4096 retained tail.',
        support='Directed nearest-neighbor distance to the full2m preview, divided by its extent norm; all endpoints retained.',
        scope='Finite approach from a declared initial volume, not a trapping-region theorem or basin probability.')
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for n in ['approach_yogi_grid.py','yogi_map.py','render_twist.py']:(out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    canvas=Image.new('RGB',(2850,1430),'#071320');draw=ImageDraw.Draw(canvas)
    large=ImageFont.truetype(findfont('DejaVu Sans'),48);small=ImageFont.truetype(findfont('DejaVu Sans'),26)
    draw.text((40,25),'YOGI / FROM A VOLUME TO PERSISTENT STRUCTURE',font=large,fill='#e8f0f5')
    draw.text((40,100),'32,768 initial states per row / same states followed through time / complete ensembles',font=small,fill='#9ab0c1')
    rows=[];cohort_p=[];cohort_starts=[]
    for row,(ident,p,box) in enumerate(zip(ids,parameters,boxes)):
        starts=box[0]+unit*(box[1]-box[0]);snapshots,tails=evolve(p,starts,times)
        assert np.isfinite(snapshots).all() and np.isfinite(tails).all()
        pr=next(r for r in preview if r['id']==ident);basis=np.array(pr['pca_basis_columns']);center=np.array(pr['center'])
        q=(snapshots-center)@basis;halfwidth=float(np.max(np.linalg.norm(q,axis=2))*1.045);projection=camera(0 if ident==670 else np.pi/4,.25)
        top=230+row*565;draw.text((40,top-50),f'G{ident:04d} / '+('FEATHERED ARC' if ident==670 else 'TWO-DIRECTION CANDIDATE'),font=small,fill='#e8f0f5')
        for j,t in enumerate(times):
            h=density(q[j],projection,535,halfwidth);assert h.sum()==len(starts)
            im,_=colorize(h,vmax=np.log1p(45));left=30+565*j;canvas.paste(im,(left,top))
            draw.text((left+14,top+14),f'n = {t:,}',font=small,fill='#e8f0f5')
        reference=np.load(args.preview/f'case-{ident}.npy');tree=cKDTree(reference);scale=float(np.linalg.norm(np.ptp(reference,axis=0)))
        distance=tree.query(snapshots[-1])[0]/scale
        fixed=np.max(np.linalg.norm(np.diff(tails,axis=1),axis=2),axis=1)<1e-10*(1+np.max(abs(tails),axis=(1,2)))
        np.savez_compressed(out/f'case-{ident}.npz',parameters=p,starts=starts,times=times,snapshots=snapshots,tails=tails,endpoint_distances=distance,fixed_tail=fixed)
        rows.append(dict(grid_id=ident,endpoint_distances_q50_q95_q99_max=np.quantile(distance,[.5,.95,.99,1]).tolist(),
            detected_fixed_tails=int(fixed.sum()),farther_than_002_and_nonfixed=int(np.sum((distance>.02)&~fixed)),
            center=center.tolist(),basis=basis.tolist(),projection=projection.tolist(),halfwidth=halfwidth,clipped_points=0,
            reference=str(args.preview/f'case-{ident}.npy'),reference_sha256=hashlib.sha256((args.preview/f'case-{ident}.npy').read_bytes()).hexdigest(),reference_extent_norm=scale))
        cohort_p.extend([p.copy() for _ in range(128)]);cohort_starts.extend(starts[:128]);print('Approach',ident,rows[-1]['detected_fixed_tails'],rows[-1]['farther_than_002_and_nonfixed'],flush=True)
    draw.text((40,1360),'One fixed camera and spatial scale per row at all times. Competing fixed outcomes are included.',font=small,fill='#9ab0c1')
    canvas.save(out/'yogi-grid-approach.png')
    windows,tails,ends,margins=spectra(np.array(cohort_p),np.array(cohort_starts),100000,200000,50000,4096)
    assert np.isfinite(windows).all() and np.isfinite(tails).all()
    np.savez_compressed(out/'volume-checks.npz',parameters=cohort_p,starts=cohort_starts,windows=windows,tails=tails,ends=ends,minimum_margins=margins)
    for i,r in enumerate(rows):
        w=windows[i*128:(i+1)*128];r.update(volume_checks=128,all_leading_windows_above005=int(np.sum(np.all(w[:,:,0]>.005,axis=1))),all_second_windows_above001=int(np.sum(np.all(w[:,:,1]>.001,axis=1))))
        print('Volume rates',r['grid_id'],r['all_leading_windows_above005'],r['all_second_windows_above001'],flush=True)
    (out/'analysis.json').write_text(json.dumps(dict(runs=rows),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['source','preview','output']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
