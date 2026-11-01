"""Native-state views of a cubical trapping candidate and finite approach runs.

The illustrations do not identify the numerical cloud with the maximal
invariant attracting set. All selected starts, checkpoints, projection counts
and rigid camera coordinates are retained.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from numba import njit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from scipy.spatial import cKDTree
from search_raar_circles import step
from render_twist import camera

INK, BLUE, GOLD = '#203d4b', '#167b9a', '#db8a24'


def sha(path):
    d = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            d.update(b)
    return d.hexdigest()


@njit(cache=False)
def approach(starts, p, times, mask, shape, origin, h, offset):
    saved = np.empty((len(times), len(starts), 3))
    minimum, outside, failed = 1e100, 0, 0
    ny, nz = shape[1], shape[2]
    for k in range(len(starts)):
        x = starts[k].copy(); at = 0
        for t in range(times[-1]+1):
            if t == times[at]:
                saved[at, k] = x; at += 1
            index = np.floor((x-origin)/h).astype(np.int64)+offset
            if np.any(index < 0) or np.any(index >= shape):
                outside += 1
            else:
                key = (index[0]*ny+index[1])*nz+index[2]
                if (mask[key//8] & (1 << (key%8))) == 0:
                    outside += 1
            if t == times[-1]: break
            x, _, axis, flag = step(x, p)
            minimum = min(minimum, axis); failed += flag
    return saved, minimum, outside, failed


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[1]
    protocol = json.loads((args.input/'protocol.json').read_text())
    mask = np.load(args.input/'mask.npy', mmap_mode='r')
    keys = np.load(args.input/'keys.npy', mmap_mode='r')
    assert sha(args.input/'mask.npy') == protocol['mask_sha256']
    assert sha(args.input/'keys.npy') == protocol['keys_sha256']
    shape = np.array(protocol['shape'], dtype=np.int64)
    origin = np.array(protocol['origin']); offset = np.array(protocol['grid_offset'])
    h = protocol['step']; p = np.array(protocol['parameters'])
    reference = repo/'results/raar-circles/render/case-1051-animation.npz'
    ref = np.load(reference); states = ref['states']; assert np.array_equal(p, ref['params'])
    projections = [np.zeros((shape[0], shape[1]), np.uint16),
                   np.zeros((shape[0], shape[2]), np.uint16),
                   np.zeros((shape[1], shape[2]), np.uint16)]
    # Count every occupied source cube exactly, one x slab at a time.
    for ix in range(shape[0]):
        first, last = ix*shape[1]*shape[2], (ix+1)*shape[1]*shape[2]
        bits = np.unpackbits(np.asarray(mask[first//8:(last+7)//8]), bitorder='little')
        slab = bits[first%8:first%8+last-first].reshape(shape[1], shape[2])
        projections[0][ix] = slab.sum(1); projections[1][ix] = slab.sum(0)
        projections[2] += slab
    assert all(int(a.sum()) == len(keys) for a in projections)
    rng = np.random.default_rng(202611011135)
    selected = np.asarray(keys[rng.integers(0, len(keys), 2048)])
    indices = np.column_stack([selected//(shape[1]*shape[2]), (selected//shape[2])%shape[1], selected%shape[2]])
    fractions = rng.random((len(selected), 3))
    starts = origin+h*(indices-offset+fractions)
    times = np.array([0,1,2,4,8,16,32,64,128,256,512,1024,2048], dtype=np.int64)
    checkpoints, axis, outside, failed = approach(starts, p, times, mask, shape, origin, h, offset)
    assert outside == failed == 0
    distances = cKDTree(states).query(checkpoints.reshape(-1,3))[0].reshape(len(times),len(starts))
    print(json.dumps({'stage':'data ready','starts':len(starts),'steps':int(times[-1]),'outside':int(outside),'minimum_sampled_axis':float(axis)}),flush=True)
    view = next(r for r in json.loads((repo/'results/raar-circles/render/analysis.json').read_text())['runs'] if r['id']==1051)
    center = np.array(view['center']); rotation = np.array(view['basis'])@camera(view['azimuth'],view['elevation']).T
    assert np.max(np.abs(rotation.T@rotation-np.eye(3))) < 1e-12
    xy = (states-center)@rotation; projected = (checkpoints-center)@rotation
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':'#b8cbd2','xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42})
    palette = LinearSegmentedColormap.from_list('arch',['#ffffff','#d1e5ec','#69adc1','#17617e','#102f41'])
    stored = {}
    # Literal coordinate projections: pale pixels are occupied columns of Q.
    fig, axes = plt.subplots(1,3,figsize=(13.5,5.3))
    fig.subplots_adjust(left=.065,right=.985,bottom=.14,top=.79,wspace=.25)
    for ax, counts, pair in zip(axes,projections,[(0,1),(0,2),(1,2)]):
        a,b=pair; bounds0=origin+h*(-offset); bounds1=origin+h*(shape-offset)
        extent=[bounds0[a],bounds1[a],bounds0[b],bounds1[b]]
        rgba=np.zeros((*counts.T.shape,4)); rgba[...,:3]=[.75,.86,.89];rgba[...,3]=(counts.T>0)*.8
        ax.imshow(rgba,origin='lower',extent=extent,interpolation='nearest',rasterized=True)
        hist,xe,ye=np.histogram2d(states[:,a],states[:,b],bins=1100)
        overlay=np.ma.masked_where(hist.T==0,np.log1p(hist.T))
        ax.imshow(overlay,origin='lower',extent=[xe[0],xe[-1],ye[0],ye[-1]],interpolation='nearest',cmap=palette,norm=PowerNorm(.65),rasterized=True)
        occupied=np.argwhere(counts>0); low=origin[[a,b]]+h*(occupied.min(0)-offset[[a,b]])
        high=origin[[a,b]]+h*(occupied.max(0)+1-offset[[a,b]])
        pad=.04*(high-low); ax.set_xlim(low[0]-pad[0],high[0]+pad[0]);ax.set_ylim(low[1]-pad[1],high[1]+pad[1])
        ax.set_aspect('equal');ax.set_xlabel('xyz'[a]);ax.set_ylabel('xyz'[b]);ax.spines[['top','right']].set_visible(False)
        stored[f'projection_{a}{b}']=counts
        stored[f'reference_histogram_{a}{b}']=hist.astype(np.uint32)
        stored[f'reference_xedges_{a}{b}']=xe;stored[f'reference_yedges_{a}{b}']=ye
    fig.text(.045,.94,'A FINITE VOLUME AROUND THE SERRATED ARCH',weight='bold',fontsize=18)
    fig.text(.045,.875,'Pale blue: projections of all 21,405,802 cubes in Q.  Dark blue: 1,000,000 actual reference states.',fontsize=10.5,color='#55727f')
    fig.text(.045,.035,'Q is a surrounding region, not the attractor itself. Coordinate projections overlap depth; they do not show holes in three dimensions.',fontsize=9,color='#55727f')
    for ext in ['png','pdf']:fig.savefig(args.out/f'raar-trapping-region.{ext}',dpi=220,facecolor='white')
    plt.close(fig)
    # Same rigid view, same equal coordinate scale in all three panels.
    low=np.minimum(xy[:,:2].min(0),projected[0,:,:2].min(0));high=np.maximum(xy[:,:2].max(0),projected[0,:,:2].max(0))
    pad=.04*(high-low);low-=pad;high+=pad
    hist,xe,ye=np.histogram2d(xy[:,0],xy[:,1],bins=(1200,1200),range=list(zip(low,high)))
    fig,axes=plt.subplots(1,3,figsize=(13.5,5.8));fig.subplots_adjust(left=.025,right=.99,bottom=.14,top=.78,wspace=.025)
    for ax,t in zip(axes,[0,16,256]):
        pos=int(np.where(times==t)[0][0])
        ax.imshow(np.log1p(hist.T),origin='lower',extent=[low[0],high[0],low[1],high[1]],cmap=palette,interpolation='nearest',norm=PowerNorm(.7),alpha=.4,rasterized=True)
        ax.scatter(projected[pos,:,0],projected[pos,:,1],s=2.8,c=GOLD,alpha=.75,linewidths=0,rasterized=True)
        ax.set_xlim(low[0],high[0]);ax.set_ylim(low[1],high[1]);ax.set_aspect('equal');ax.axis('off')
        ax.set_title('VOLUME OF STARTS' if t==0 else f'AFTER {t} ITERATIONS',fontsize=11,weight='bold',pad=10)
    fig.text(.04,.94,'HOW THE CLOUD FORMS FROM A VOLUME OF STARTS',weight='bold',fontsize=18)
    fig.text(.04,.88,'2,048 starts sampled from Q. Same states evolved by the original RAAR map; one rigid camera, no shape warping.',fontsize=10.2,color='#55727f')
    fig.text(.04,.05,'Gold: current numerical iterates. Pale blue: finite reference cloud. This illustration does not identify that cloud with the entire attracting set.',fontsize=9,color='#55727f')
    for ext in ['png','pdf']:fig.savefig(args.out/f'raar-trapping-approach.{ext}',dpi=230,facecolor='white')
    plt.close(fig)
    np.savez_compressed(args.out/'figure-data.npz',selected_keys=selected,fractions=fractions,starts=starts,times=times,checkpoints=checkpoints,reference_distances=distances,center=center,rotation=rotation,rigid_reference_histogram=hist.astype(np.uint32),rigid_xedges=xe,rigid_yedges=ye,**stored)
    result={'protocol_sha256':sha(args.input/'protocol.json'),'reference_sha256':sha(reference),'source_sha256':sha(Path(__file__)),
            'cell_count':len(keys),'projection_counts':[int(v.sum()) for v in projections], 'starts':len(starts),'times':times.tolist(),
            'membership_checks':len(starts)*(int(times[-1])+1),'numerical_outside_count':int(outside),'projection_failures':int(failed),
            'minimum_sampled_axis':float(axis),'reference_distance_quantiles':np.quantile(distances,[.1,.5,.9,.99],axis=1).T.tolist(),
            'scope':'Exact integer counts of occupied grid columns; actual finite numerical iterations and rigid projections. Grid membership of sampled states is a binary64 diagnostic. No identification of the cloud with the attracting invariant set, basin completeness, transitivity or novelty claim.'}
    (args.out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    main(parser.parse_args())
