"""Unaligned geometric comparison of transcription and actual SciPy states."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from render_twist import camera,density,colorize


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);rows=[];clouds={}
    for ident in [335,448,1890,2006]:
        source=args.gallery/f'case-{ident}.npz';actual=args.upstream/f'case-{ident}.npz'
        with np.load(source) as f:reference=f['states']
        with np.load(actual) as f:upstream=f['states'][-200000:,:3]
        a=cKDTree(reference).query(upstream,workers=4)[0];b=cKDTree(upstream).query(reference,workers=4)[0]
        floor=cKDTree(reference[:100000]).query(reference[100000:],workers=4)[0]
        levels=[.1,.5,.9,.99,1]
        row=dict(id=ident,reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),upstream_sha256=hashlib.sha256(actual.read_bytes()).hexdigest(),points_each=200000,
            actual_to_reference_quantiles=np.quantile(a,levels).tolist(),reference_to_actual_quantiles=np.quantile(b,levels).tolist(),reference_halfcloud_quantiles=np.quantile(floor,levels).tolist(),quantile_levels=levels)
        np.savez_compressed(args.out/f'case-{ident}.npz',actual_to_reference=a,reference_to_actual=b,reference_halfcloud=floor)
        rows.append(row);clouds[ident]=(reference,upstream)
    render=json.loads((args.render/'analysis.json').read_text())['runs'];fig,axes=plt.subplots(2,2,figsize=(12,10),facecolor='#071320');meta=[]
    for row,ident in enumerate([335,448]):
        r=next(v for v in render if v['id']==ident);center=np.array(r['center']);basis=np.array(r['basis']);proj=camera(r['azimuth'],r['elevation'])
        q=[(v-center)@basis for v in clouds[ident]];half=max(np.linalg.norm(v,axis=1).max() for v in q)*1.04
        counts=[]
        for col,(z,label) in enumerate(zip(q,['Mathematical transcription','Unmodified SciPy 1.18.1'])):
            h=density(z,proj,850,half);assert int(h.sum())==len(z);im,_=colorize(h,vmax=np.log1p(35));ax=axes[row,col];ax.imshow(im);ax.set_axis_off();ax.set_title(f'D{ident} / {label}',fontsize=13,color='#e8f0f5');counts.append(int(h.sum()))
        meta.append(dict(id=ident,basis=r['basis'],center=r['center'],projection=proj.tolist(),halfwidth=float(half),counts=counts))
    fig.text(.04,.965,'SAME COORDINATES / INDEPENDENT LIBRARY CLOUDS',fontsize=20,color='#e8f0f5',weight='bold')
    fig.text(.04,.025,'200,000 consecutive late states per panel. One rigid camera, common scale and density colors in each row.\nNo fitting, alignment, rescaling or nonlinear deformation between the two calculations.',fontsize=11,color='#b5c8d8')
    fig.subplots_adjust(top=.91,bottom=.08,left=.02,right=.98,hspace=.15,wspace=.03);fig.savefig(args.out/'rk45-library-clouds.png',dpi=160,facecolor='#071320');plt.close(fig)
    report=dict(runs=rows,render=meta,scope='Unaligned finite-sample geometry comparison in equal actual-state units; small sampling distances are not proof of identical invariant sets or infinite-time shadowing.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for n in ['compare_rk45_clouds.py','render_twist.py']:(args.out/n).write_bytes(Path(__file__).with_name(n).read_bytes())
    print(json.dumps(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--gallery',type=Path,required=True);p.add_argument('--upstream',type=Path,required=True);p.add_argument('--render',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
