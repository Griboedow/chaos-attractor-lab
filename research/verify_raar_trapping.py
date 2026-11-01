"""Reconstruct exact inputs, audit exhaustive coverage and replay figure data.

This is not a second implementation of all 21 million interval enclosures.
Run the declared batch commands to recompute those enclosures.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.spatial import cKDTree
import audit_raar_trapping_closure as closure
import restore_raar_trapping as restore
import verify_raar_symbolic as symbolic
from verify_raar_cycles import inclusion, interval
from render_raar_trapping import approach
from audit_raar_trapping_batch import file_hash as sha


def main(args):
    repo=Path(__file__).resolve().parents[1];packet=repo/'results/raar-trapping';runs=packet/'runs'
    args.out.mkdir(parents=True,exist_ok=False)
    read=lambda p:json.loads(p.read_text(encoding='utf-8'))
    manifest=read(packet/'manifest.json')['files']
    actual={p.relative_to(packet).as_posix() for p in packet.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts}
    assert actual=={v['path'] for v in manifest}
    for row in manifest:
        p=packet/row['path'];assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],p
    for row in read(packet/'provenance.json'):
        assert sha(packet/row['path'])==row['sha256']
    info=read(packet/'packet.json')
    assert sha(repo/'results/raar-symbolic-dynamics/manifest.json')==info['symbolic_manifest_sha256']
    restore.main(SimpleNamespace(packet=packet,out=args.out/'restored'))
    closure.main(SimpleNamespace(input=args.out/'restored',batches=[runs/f'raar-trapping-audit-round2-part{i}' for i in range(4)],out=args.out/'accounting'))
    current=read(args.out/'accounting/report.json');saved=read(runs/'raar-trapping-closure-audit01/report.json')
    assert current==saved
    assert current['status']=='closed_accounting_verified' and current['new_unprocessed_cells']==0
    symbolic.main(SimpleNamespace(scratch=args.out,output=args.out/'symbolic-dependency.json'))
    cycle7=repo/'results/raar-cycle-structure/cycle-proof/base-unstable-period-7.json'
    inclusion(read(cycle7))
    structure=read(repo/'results/raar-cycle-structure/structure-audit/report.json')
    assert structure['input_sha256']==sha(cycle7) and interval(structure['determinant'])[1]<0
    a=np.load(runs/'raar-interval-grid-validation05/sample.npz');b=np.load(runs/'raar-interval-grid-validation06/sample.npz')
    assert a.files==b.files and all(np.array_equal(a[k],b[k]) for k in a.files)
    figures=runs/'raar-trapping-render02';data=np.load(figures/'figure-data.npz');render=read(figures/'analysis.json')
    protocol=read(packet/'candidate/protocol.json');mask=np.load(args.out/'restored/mask.npy',mmap_mode='r')
    shape=np.array(protocol['shape'],dtype=np.int64);origin=np.array(protocol['origin']);offset=np.array(protocol['grid_offset']);h=protocol['step']
    assert sha(packet/'candidate/protocol.json')==render['protocol_sha256']
    counts=[np.zeros((shape[0],shape[1]),np.uint16),np.zeros((shape[0],shape[2]),np.uint16),np.zeros((shape[1],shape[2]),np.uint16)]
    for x in range(shape[0]):
        lo,hi=x*shape[1]*shape[2],(x+1)*shape[1]*shape[2]
        bits=np.unpackbits(np.asarray(mask[lo//8:(hi+7)//8]),bitorder='little')
        slab=bits[lo%8:lo%8+hi-lo].reshape(shape[1],shape[2])
        counts[0][x]=slab.sum(1);counts[1][x]=slab.sum(0);counts[2]+=slab
    reference=repo/'results/raar-circles/render/case-1051-animation.npz';assert sha(reference)==render['reference_sha256']
    states=np.load(reference)['states']
    for count,pair in zip(counts,[(0,1),(0,2),(1,2)]):
        suffix=f'{pair[0]}{pair[1]}'
        assert np.array_equal(count,data['projection_'+suffix]) and int(count.sum())==21405802
        hist=np.histogram2d(states[:,pair[0]],states[:,pair[1]],bins=[data['reference_xedges_'+suffix],data['reference_yedges_'+suffix]])[0]
        assert np.array_equal(hist,data['reference_histogram_'+suffix]) and int(hist.sum())==len(states)
    assert data['rotation'].shape==(3,2) and np.max(abs(data['rotation'].T@data['rotation']-np.eye(2)))<1e-12
    xy=(states-data['center'])@data['rotation']
    hist=np.histogram2d(xy[:,0],xy[:,1],bins=[data['rigid_xedges'],data['rigid_yedges']])[0]
    assert np.array_equal(hist,data['rigid_reference_histogram']) and int(hist.sum())==1000000
    keys=data['selected_keys'];indices=np.column_stack([keys//(shape[1]*shape[2]),(keys//shape[2])%shape[1],keys%shape[2]])
    starts=origin+h*(indices-offset+data['fractions']);assert np.array_equal(starts,data['starts'])
    assert np.all(((mask[keys//8]>>(keys%8))&1)==1)
    checkpoints,axis,outside,failed=approach(starts,np.array(protocol['parameters']),data['times'],mask,shape,origin,h,offset)
    assert np.array_equal(checkpoints,data['checkpoints']) and axis==render['minimum_sampled_axis'] and outside==failed==0
    distances=cKDTree(states).query(checkpoints.reshape(-1,3))[0].reshape(len(data['times']),len(starts))
    assert np.array_equal(distances,data['reference_distances'])
    for name in ['region','approach']:
        assert sha(figures/f'raar-trapping-{name}.png')==sha(repo/'docs/assets'/f'raar-trapping-{name}.png')
    result={'status':'passed','files':len(manifest),'manifest_sha256':sha(packet/'manifest.json'),
            'source_cells':current['source_cells'],'accepted_leaves':current['accepted_leaves'],'axis_lower_bound':current['axis_lower_bound'],
            'image_margin_exact':current['image_margin_exact'],'figure_starts':len(starts),'figure_iteration_membership_checks':render['membership_checks'],
            'scope':'Exact input restoration, exhaustive accounting, prior symbolic certificate-output audit, located cycle inclusion, and full finite figure-data replay. Does not independently recompute every interval enclosure; use the full batch replay commands for that.'}
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args())
