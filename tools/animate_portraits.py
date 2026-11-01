"""Slow rigid-camera rotations of the seven main-page state portraits."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,time
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.font_manager import findfont
from render_portraits import ROOT,COLOR,coordinates,camera,digest,panel_label
from data_cache import data_path

WIDTH,HEIGHT=800,480
PANEL=360
FRAMES,DURATION=600,40
QUARTERS=[0,FRAMES//4,FRAMES//2,3*FRAMES//4]
PALETTE=(COLOR(np.linspace(0,1,64))[:,:3]*255).round().astype(np.uint8)
PALETTE[0]=255

def full_camera(projection):
    p=np.asarray(projection)[:2]
    C=np.vstack([p,np.cross(p[0],p[1])])
    assert np.max(abs(C@C.T-np.eye(3)))<2e-12 and np.linalg.det(C)>.999999999
    return C

def old_camera(az,el):
    return np.array([[np.cos(az),0,np.sin(az)],[np.sin(az)*np.sin(el),np.cos(el),-np.cos(az)*np.sin(el)]])

def load_views(row):
    views=[];records=[];family=row['id']
    if family=='icosahedral':
        source=ROOT/'results/portrait-presentation/icosa-animation.npz'
        with np.load(data_path(source)) as d:
            for ident in [177,183]:
                z=d[f'points_{ident}'];cam_source=ROOT/f'results/icosahedral-flow/icosa-final/portrait-{ident}.npz'
                with np.load(data_path(cam_source)) as f:C=full_camera(f['camera'])
                views.append(z@C.T);records.append({'source':source.relative_to(ROOT).as_posix(),'key':f'points_{ident}','source_sha256':digest(source),'camera_source':cam_source.relative_to(ROOT).as_posix(),'camera_source_sha256':digest(cam_source),'camera':C.tolist(),'samples':len(z),'mark':'density','selection':'Every fourth state of the original dense 16000-unit trajectory; bitwise checked against the frozen sample.'})
    elif family=='anderson':
        metadata_path=ROOT/'results/unexpected-search/visuals/provenance.json';metadata=json.loads(metadata_path.read_text())
        for ident in [54,449]:
            source=ROOT/f'results/unexpected-search/visuals/anderson-{ident}-sample.npz'
            with np.load(data_path(source)) as d:z=d['points']
            meta=next(r for r in metadata if r.get('file')==f'anderson-{ident}.png');C=full_camera(meta['camera']);q=(z-np.array(meta['center']))@C.T
            views.append(q);records.append({'source':source.relative_to(ROOT).as_posix(),'source_sha256':digest(source),'metadata_source':metadata_path.relative_to(ROOT).as_posix(),'metadata_sha256':digest(metadata_path),'samples':len(z),'camera':C.tolist(),'center':meta['center'],'mark':'density','selection':'All states in the archived stride-32 sample of the 3200000-state rendering.'})
    elif family in ['raar','rk45','yogi-tanh']:
        folder={'raar':'raar-circles/render','rk45':'rk45-controller/oscillator/render','yogi-tanh':'yogi-tanh/render'}[family]
        metadata_path=ROOT/'results'/folder/'analysis.json';rows=json.loads(metadata_path.read_text())['runs']
        for spec in row['panels']:
            ident=int(spec['label']);meta=next(r for r in rows if r['id']==ident);source=ROOT/'results'/folder/f'case-{ident}-animation.npz'
            with np.load(data_path(source)) as data:z=data['states'][-250000:]
            basis=np.array(meta['basis']);C=full_camera(old_camera(meta['azimuth'],meta['elevation']))
            assert np.max(abs(basis.T@basis-np.eye(3)))<1e-12 and np.linalg.det(basis)>.999999999
            q=((z-np.array(meta['center']))@basis)@C.T;views.append(q)
            records.append({'source':source.relative_to(ROOT).as_posix(),'source_sha256':digest(source),'metadata_source':metadata_path.relative_to(ROOT).as_posix(),'metadata_sha256':digest(metadata_path),'samples':len(z),'center':meta['center'],'basis':basis.tolist(),'camera':C.tolist(),'mark':'density','selection':'Last 250000 consecutive states of the archived unthinned one-million-state animation input.'})
    elif family in ['motor-four','motor-eight']:
        z=[coordinates(spec) for spec in row['panels']];combined=np.vstack(z);center=(combined.min(0)+combined.max(0))/2;C=full_camera(camera())
        for spec,points in zip(row['panels'],z):
            views.append((points-center)@C.T)
            records.append({'input':spec,'source':spec['source'],'source_sha256':digest(ROOT/spec['source']),'samples':len(points),'center':center.tolist(),'camera':C.tolist(),'mark':spec['mark'],'selection':'Same saved states and physical coordinates as the static presentation.'})
    else:raise ValueError(family)
    radii=[float(np.linalg.norm(q,axis=1).max()*1.045) for q in views]
    if row.get('common_scale'):radii=[max(radii)]*len(radii)
    for record,q,radius in zip(records,views,radii):
        record['radius']=radius;record['view_coordinates_sha256']=hashlib.sha256(q.astype('<f8').tobytes()).hexdigest()
    return views,records

def projected(q,angle):
    return np.column_stack([np.cos(angle)*q[:,0]+np.sin(angle)*q[:,2],q[:,1]])

def histogram(q,angle,radius):
    xy=projected(q,angle);bins=np.floor(xy/(2*radius)*PANEL+PANEL/2).astype(np.int32)
    assert np.all((bins>=0)&(bins<PANEL))
    h=np.bincount(bins[:,1]*PANEL+bins[:,0],minlength=PANEL*PANEL).reshape(PANEL,PANEL)[::-1]
    assert h.sum()==len(q)
    return h

def render(row,output):
    begin=time.time();views,records=load_views(row);font=ImageFont.truetype(findfont('DejaVu Sans'),19)
    for q,record in zip(views,records):
        probes=[histogram(q,a,record['radius']) for a in [0,np.pi/3,2*np.pi/3]]
        record['log_saturation']=max(float(np.quantile(np.log1p(h[h>0]),.9995)) for h in probes)
    frames=[];counts=[];pixel_bounds=[]
    for k,angle in enumerate(np.linspace(0,2*np.pi,FRAMES,endpoint=False)):
        frame=Image.new('P',(WIDTH,HEIGHT),0);frame.putpalette(np.pad(PALETTE,((0,192),(0,0)),constant_values=255).ravel().tolist());draw=ImageDraw.Draw(frame)
        row_counts=[]
        for j,(q,record,spec) in enumerate(zip(views,records,row['panels'])):
            left=20+400*j;top=80;mark=record['mark']
            label=panel_label(spec);box=draw.textbbox((0,0),label,font=font);draw.text((left+PANEL/2-(box[2]-box[0])/2,24),label,font=font,fill=63)
            record['display_label']=label
            if mark=='density':
                h=histogram(q,angle,record['radius']);values=np.ceil(np.minimum(np.log1p(h)/record['log_saturation'],1)**.72*63).astype(np.uint8)
                tile=Image.fromarray(values,mode='P');tile.putpalette(frame.getpalette());frame.paste(tile,(left,top));row_counts.append(int(h.sum()))
            else:
                xy=projected(q,angle)/record['radius']*(PANEL/2);xy[:,0]+=left+PANEL/2;xy[:,1]=top+PANEL/2-xy[:,1]
                assert np.all((xy[:,0]>=left)&(xy[:,0]<left+PANEL)&(xy[:,1]>=top)&(xy[:,1]<top+PANEL))
                if mark=='line':draw.line([tuple(v) for v in xy],fill=48,width=1)
                else:
                    for x,y in xy:draw.ellipse((x-1,y-1,x+1,y+1),fill=43)
                row_counts.append(len(q))
            if k in QUARTERS:
                xy=projected(q,angle);pixel_bounds.append({'frame':k,'panel':j,'minimum':xy.min(0).tolist(),'maximum':xy.max(0).tolist()})
        frames.append(frame);counts.append(row_counts)
        if k%100==0:print(json.dumps({'family':row['id'],'frame':k,'seconds':time.time()-begin}),flush=True)
    path=output/(row['id']+'.webp')
    frames[0].save(path,format='WEBP',save_all=True,append_images=frames[1:],duration=DURATION,loop=0,lossless=True,method=1)
    with Image.open(path) as gif:
        assert gif.n_frames==FRAMES and gif.info['loop']==0
        durations=[]
        for k in range(gif.n_frames):gif.seek(k);gif.load();durations.append(gif.info['duration'])
    assert sum(durations)==FRAMES*DURATION
    result={'id':row['id'],'file':path.relative_to(ROOT).as_posix(),'sha256':digest(path),'bytes':path.stat().st_size,'frames':FRAMES,'duration_ms':durations,'seconds_per_revolution':sum(durations)/1000,'frame_point_counts':counts,'checkpoint_projection_bounds':pixel_bounds,'panels':records,'elapsed_seconds':time.time()-begin}
    assert result['bytes']<50_000_000
    subprocess.run([sys.executable,str(ROOT/'tools/render_dates.py'),'record',result['file'],'--source-portrait',row['source_portrait']],check=True)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--ids',nargs='*');args=parser.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    assert output.is_relative_to(ROOT),'Publication animation output must be inside this repository.'
    catalog=json.loads((ROOT/'provenance/portrait-catalog.json').read_text(encoding='utf-8'));rows=[r for r in catalog['rows'] if r.get('featured') and (not args.ids or r['id'] in args.ids)]
    results=[]
    for row in rows:results.append(render(row,output))
    manifest={'schema':1,'renderer_sha256':digest(Path(__file__)),'static_renderer_sha256':digest(ROOT/'tools/render_portraits.py'),'catalog_sha256':digest(ROOT/'provenance/portrait-catalog.json'),'size':[WIDTH,HEIGHT],'format':'Animated WebP, lossless','palette_size':64,'motion':'Rigid camera rotation about the initial vertical axis. State arrays, labels, point selection, physical units and radial bounds are fixed throughout. No time evolution or geometry morphing.','angle':f'2*pi*frame/{FRAMES}, endpoint excluded','rows':results}
    manifest['data_resolver_sha256']=digest(ROOT/'tools/data_cache.py')
    (output/'animation.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'animations':len(results),'total_bytes':sum(r['bytes'] for r in results)}),flush=True)

if __name__=='__main__':main()
