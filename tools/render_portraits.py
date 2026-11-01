"""Render a consistent scientific atlas from archived counts and state arrays."""
from pathlib import Path
import argparse,hashlib,json,math,subprocess,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT=Path(__file__).resolve().parents[1]
CATALOG=ROOT/'provenance/portrait-catalog.json'
COLOR=LinearSegmentedColormap.from_list('portrait_ink',['#ffffff','#70a9b7','#164e63','#122f42'])

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def coordinates(spec):
    with np.load(ROOT/spec['source']) as data:
        x=data[spec['key']].copy();mode=spec.get('transform','identity')
        if 'time_key' in spec:
            t=data[spec['time_key']];mask=np.ones(len(t),bool)
            if 'last_duration' in spec:mask &= t>=t[-1]-spec['last_duration']
            if 'after' in spec:mask &= t>spec['after']
            x=x[mask]
        if 'last_count' in spec:x=x[-spec['last_count']:]
        x=x[::spec.get('stride',1)]
        if mode=='motor_cosine':
            matrix=np.array([[1,-1,0],[1,1,-2],[1,1,1]],float)/np.sqrt([[2],[6],[3]])
            x=x[:,:3]@matrix.T
        elif mode=='motor_pair':
            p=data['parameters'];x=np.c_[p[1]*(x[:,0]+x[:,1]),x[:,1:]]
        elif mode=='synchronematic':
            n=x.shape[1]//2;theta=x[:,:n];chi=x[:,n:]
            a=np.exp(1j*(theta+chi)).mean(1);b=np.exp(1j*(theta-chi)).mean(1)
            Q=np.exp(2j*chi).mean(1);C=np.sqrt((abs(a)**2+abs(b)**2)/2)
            x=np.c_[Q.real,Q.imag,C]
        elif mode=='lce':
            p=data['parameters'];q=p[2]/(1+p[2]);x=x[:,1:4]/np.array([q,np.sqrt(p[0])*q,q])
        elif mode=='two_strain':x=x[:,:3]*np.array([1,1000,1000])
        elif mode=='optical':x=np.c_[x[:,0],x[:,3],x[:,2]-x[:,5]]
        elif mode!='identity':raise ValueError(mode)
        if 'columns' in spec:x=x[:,spec['columns']]
    assert x.ndim==2 and x.shape[1]==3 and np.isfinite(x).all()
    return x

def camera(elevation=27,azimuth=-45):
    e,a=np.deg2rad([elevation,azimuth])
    return np.array([[-np.sin(a),np.cos(a),0],[-np.cos(a)*np.sin(e),-np.sin(a)*np.sin(e),np.cos(e)]])

def prepare(row):
    records=[];arrays=[]
    points=[coordinates(p) if p['mode']=='points' else None for p in row['panels']]
    common=None
    if row.get('common_scale'):
        assert all(p is not None for p in points)
        allpoints=np.vstack(points);common=(allpoints.min(0)+allpoints.max(0))/2
        radius=max(float(np.max(abs((p-common)@camera().T))) for p in points)*1.06
    for spec,z in zip(row['panels'],points):
        record={'input':spec,'input_sha256':digest(ROOT/spec['source'])}
        if spec['mode']=='histogram':
            with np.load(ROOT/spec['source']) as data:
                h=data[spec['key']]
                if 'index' in spec:h=h[spec['index']]
            assert h.ndim==2 and np.isfinite(h).all() and np.all(h>=0)
            record.update(original_shape=list(h.shape),original_count=int(h.sum()))
            occupied=np.argwhere(h>0);lo=occupied.min(0);hi=occupied.max(0)+1
            h=h[lo[0]:hi[0],lo[1]:hi[1]]
            margin=max(1,int(np.ceil(max(h.shape)*.04)))
            h=np.pad(h,margin)
            record.update(nonempty_crop=[lo.tolist(),hi.tolist()],empty_margin=margin)
            # Sum whole square pixel blocks; pad only with empty bins. No stretching.
            factor=max(1,math.ceil(max(h.shape)/(1400 if len(row['panels'])==1 else 760)))
            height,width=h.shape;pad=(-height%factor,-width%factor)
            h=np.pad(h,((0,pad[0]),(0,pad[1]))).reshape((height+pad[0])//factor,factor,(width+pad[1])//factor,factor).sum((1,3))
            record.update(square_rebin_factor=factor,zero_padding=list(pad))
        else:
            center=common if common is not None else (z.min(0)+z.max(0))/2
            projection=camera();xy=(z-center)@projection.T
            half=radius if common is not None else float(abs(xy).max()*1.06)
            h=np.histogram2d(xy[:,0],xy[:,1],bins=760,range=[[-half,half]]*2)[0].T[::-1]
            assert int(h.sum())==len(z)
            record.update(original_count=len(z),coordinate_sha256=hashlib.sha256(z.astype('<f8').tobytes()).hexdigest(),center=center.tolist(),camera=projection.tolist(),half_width=half,spatially_clipped=0)
            record['_xy']=xy
        assert int(h.sum())==record['original_count']
        arrays.append(h);records.append(record)
    return arrays,records

def render(row,output):
    arrays,records=prepare(row);fig=plt.figure(figsize=(12,6.75),dpi=160,facecolor='white')
    for k,(h,rec,spec) in enumerate(zip(arrays,records,row['panels'])):
        width=.94/len(arrays);ax=fig.add_axes([.03+k*width,.07,width-.02,.80])
        levels=np.log1p(h);positive=levels[h>0];vmax=float(np.quantile(positive,.9995))
        intensity=np.minimum(levels/vmax,1)**.72
        xy=rec.pop('_xy',None);mark=spec.get('mark','density')
        if mark in ['line','points']:
            assert xy is not None
            if mark=='line':ax.plot(xy[:,0],xy[:,1],color='#164e63',lw=.65,alpha=.78)
            else:ax.scatter(xy[:,0],xy[:,1],color='#164e63',s=1.,alpha=.8,linewidths=0)
            half=rec['half_width'];ax.set(xlim=(-half,half),ylim=(-half,half));ax.set_aspect('equal')
        else:ax.imshow(intensity,cmap=COLOR,vmin=0,vmax=1,interpolation='nearest',aspect='equal',origin='upper')
        ax.set_axis_off();fig.text(.03+k*width+(width-.02)/2,.925,spec['label'],ha='center',fontsize=13,color='#253d48')
        rec.update(rendered_count=int(h.sum()),log_saturation=vmax,occupied_bins=int((h>0).sum()))
    path=output/(row['id']+'.png');fig.savefig(path,dpi=160,facecolor='white');plt.close(fig)
    return {'id':row['id'],'file':path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path),'sha256':digest(path),'panels':records}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--ids',nargs='*');parser.add_argument('--record-dates',action='store_true');args=parser.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    rows=json.loads(CATALOG.read_text(encoding='utf-8'))['rows'];records=[]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12})
    for row in rows:
        if args.ids and row['id'] not in args.ids:continue
        result=render(row,output);records.append(result)
        if args.record_dates:
            subprocess.run([sys.executable,str(ROOT/'tools/render_dates.py'),'record',result['file'],'--source-portrait',row['source_portrait']],check=True)
        print(json.dumps({'rendered':row['id'],'counts':[p['rendered_count'] for p in result['panels']]}),flush=True)
    record={'schema':1,'renderer_sha256':digest(Path(__file__)),'catalog_sha256':digest(CATALOG),'size':[1920,1080],'style':{'background':'white','palette':['#ffffff','#70a9b7','#164e63','#122f42'],'intensity':'min(log1p(count)/q0.9995,1)^0.72','count_meaning':'Sampling counts; independent saturation per panel. Brightness is not comparable across studies.','lines':'When declared, straight segments join consecutive saved continuous-flow states; no new states are inferred.','geometry':'Only empty histogram margins are removed, followed by square-bin aggregation or declared coordinate projections. No smoothing, image synthesis, anisotropic camera scaling or trajectory replication.'},'rows':records}
    (output/'render.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()
