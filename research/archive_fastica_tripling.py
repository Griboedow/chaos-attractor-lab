"""Archive compact FastICA evidence and append its numerical portraits once."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import render_portraits as static
import animate_portraits as animated

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);a=p.parse_args();src=a.input.resolve()
    assert not src.is_relative_to(ROOT),'Calculation input must be external.'
    packet=ROOT/'results/fastica-tripling';assert not packet.exists(),'Archive only into a new packet.'
    packet.mkdir();(packet/'code').mkdir()
    for name in ['model.py','controls.py','exact_check.py','long_check.py','precision_check.py','library_check.py']:
        shutil.copyfile(src/name,packet/'code'/name)
    for name in ['controls.json','exact-check.json','precision.json','library.json','pilot-01/protocol.json','pilot-01/report.json','long-01/protocol.json','long-01/report.json','long-01/fresh-starts.json','long-01/portrait-040.json','long-01/portrait-150.json']:
        dest=packet/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src/name,dest)
    presentation=src/'presentation';presentation.mkdir(exist_ok=False)
    points=[]
    for key in ['040','150']:
        with np.load(src/f'long-01/portrait-{key}.npz') as d:z=d['states'][:,:3].copy()
        assert z.shape==(500000,3) and np.isfinite(z).all();points.append(z)
        np.savez_compressed(presentation/f'points-{key}.npz',points=z)
    joined=np.vstack(points);center=(joined.min(0)+joined.max(0))/2
    for key in ['040','150']:
        write(packet/f'presentation/view-{key}.json',{'center':center.tolist(),'basis':np.eye(3).tolist(),'projection':static.camera().tolist(),'coordinates':'(A,B,C), discarding D; every consecutive state, equal units','source':f'long-01/portrait-{key}.npz','source_sha256':sha(src/f'long-01/portrait-{key}.npz')})
    registry_path=ROOT/'records/generated-data.json';registry=json.loads(registry_path.read_text());new=[]
    for rel in ['pilot-01/survey.npz','long-01/runs.npz','long-01/portrait-040.npz','long-01/portrait-150.npz','presentation/points-040.npz','presentation/points-150.npz']:
        path=src/rel;checksum=sha(path);cached=Path('D:/Codex/chaos-attractor-cache')/checksum[:2]/(checksum+'.npz');cached.parent.mkdir(parents=True,exist_ok=True)
        if cached.exists():assert sha(cached)==checksum
        else:shutil.copyfile(path,cached)
        with np.load(path) as data:arrays=[{'key':k,'shape':list(data[k].shape),'dtype':str(data[k].dtype),'bytes':data[k].nbytes,'retained':False} for k in data.files]
        row={'path':'results/fastica-tripling/'+rel,'bytes':path.stat().st_size,'sha256':checksum,'group':'fastica-tripling','compact':None,'arrays':arrays,'recovery':'recompute','recompute_command':'See results/fastica-tripling/README.md; run research/fastica_tripling.py into a new external directory. Presentation arrays are states[:, :3] from long-01/portrait-*.npz. Recomputed chaotic trajectories can differ across environments.'}
        assert row['path'] not in {r['path'] for r in registry['files']};new.append(row)
    registry['files'].extend(new);registry['removed_bytes']+=sum(r['bytes'] for r in new);write(registry_path,registry);write(packet/'omitted-data.json',new)
    with (ROOT/'.gitignore').open('a',encoding='utf-8',newline='\n') as f:
        f.write('\n# Recomputed structured FastICA data stay in the external cache.\n'+'\n'.join('/'+r['path'] for r in new)+'\n')
    with (ROOT/'.gitattributes').open('a',encoding='utf-8',newline='\n') as f:f.write('\n# Preserve executed structured FastICA sources and numerical receipts.\nresults/fastica-tripling/** -text whitespace=cr-at-eol\n')
    catalog_path=ROOT/'records/portrait-catalog.json';old_catalog_sha=sha(catalog_path);catalog=json.loads(catalog_path.read_text());old_rows=json.dumps(catalog['rows'],sort_keys=True)
    row={'id':'fastica-tripling','title':'Structured-data FastICA','group':'Numerical algorithms','report':'docs/fastica-tripling.md','source_portrait':'docs/assets/portraits/fastica-tripling.png','featured':True,'common_scale':True,
         'description':'An explicit finite-data construction has a certified attracting positive-entropy set for 0 <= epsilon <= 1/25. The epsilon=0.15 companion is numerical only. Expanding-circle dynamics are known; historical priority of this realization is unresolved.',
         'coordinates':'Projection (A,B,C) of a unit weight in S3; D is omitted. Equal units, common scale, 500000 consecutive states per panel; no sign folding.',
         'panels':[{'label':'epsilon = '+str(e),'source':f'results/fastica-tripling/presentation/points-{key}.npz','key':'points','mode':'points','transform':'identity'} for e,key in [(.04,'040'),(.15,'150')]],
         'animation_panels':[{'source':f'results/fastica-tripling/presentation/points-{key}.npz','key':'points','metadata_source':f'results/fastica-tripling/presentation/view-{key}.json','last_count':500000} for key in ['040','150']]}
    assert row['id'] not in {r['id'] for r in catalog['rows']};catalog['rows'].append(row);assert json.dumps(catalog['rows'][:-1],sort_keys=True)==old_rows;write(catalog_path,catalog)
    output=ROOT/'docs/assets/portraits';r=static.render(row,output)
    subprocess.run([sys.executable,str(ROOT/'tools/render_dates.py'),'record',r['file']],check=True)
    r.update(renderer_sha256=sha(ROOT/'tools/render_portraits.py'),data_resolver_sha256=sha(ROOT/'tools/data_cache.py'))
    extension={'study':row['id'],'previous_catalog_sha256':old_catalog_sha,'current_catalog_sha256':sha(catalog_path),'scope':'Append one study; existing catalog entries and render rows preserved. New rows record their actual renderer hashes separately from the original batch.'}
    manifest_path=output/'render.json';manifest=json.loads(manifest_path.read_text());manifest['rows'].append(r);manifest['catalog_sha256']=sha(catalog_path);manifest.setdefault('extensions',[]).append(extension);write(manifest_path,manifest)
    r=animated.render(row,output)
    r.update(renderer_sha256=sha(ROOT/'tools/animate_portraits.py'),static_renderer_sha256=sha(ROOT/'tools/render_portraits.py'),data_resolver_sha256=sha(ROOT/'tools/data_cache.py'))
    manifest_path=output/'animation.json';manifest=json.loads(manifest_path.read_text());manifest['rows'].append(r);manifest['catalog_sha256']=sha(catalog_path);manifest.setdefault('extensions',[]).append(extension);write(manifest_path,manifest)
    print(json.dumps({'packet':str(packet),'omitted_bytes':sum(r['bytes'] for r in new),'new_portraits':2}),flush=True)

if __name__=='__main__':main()
