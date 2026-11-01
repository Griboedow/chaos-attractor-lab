"""Add the checked nonuniform pair using the repository's common renderers."""
import hashlib,json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tools'))
from render_portraits import render as still
from animate_portraits import render as animation
import subprocess


def main():
    catalog_path=root/'records/portrait-catalog.json'
    catalog=json.loads(catalog_path.read_text(encoding='utf-8'))
    assert not any(r['id']=='product-irregular' for r in catalog['rows'])
    source='results/product-highorder/runs/product-irregular-dense01'
    row=dict(id='product-irregular',title='Nonuniform product posteriors',group='Numerical algorithms',report='docs/product-alphabet/manuscript.md#32-nonuniform-alphabets-and-folded-shells',source_portrait='docs/assets/experimental/product-irregular-dense/product-alphabet-6.png',description='Folded-shell numerical Cases 6 and 7: positive finite growth at 33 starts and 70/100 digits. Toroidal motifs are known; attraction, topology and historical uniqueness remain open.',featured=True,coordinates='Native discrete three-state iterates; orthogonal cameras and equal units within panels. Independent panel scales and density saturation.',panels=[dict(label=str(i),source=f'{source}/case-{i}-histogram.npz',key='histogram',mode='histogram') for i in [6,7]],animation_panels=[dict(source=f'{source}/case-{i}-states.npz',key='states',metadata_source=f'{source}/case-{i}.json',last_count=250000) for i in [6,7]])
    index=next(i for i,r in enumerate(catalog['rows']) if r['id']=='product-alphabet')+1
    catalog['rows'].insert(index,row)
    catalog_path.write_text(json.dumps(catalog,indent=2)+'\n',encoding='utf-8',newline='\n')
    output=root/'docs/assets/portraits';r=still(row,output)
    subprocess.run([sys.executable,str(root/'tools/render_dates.py'),'record',r['file'],'--source-portrait',row['source_portrait']],check=True)
    a=animation(row,output)
    digest=hashlib.sha256(catalog_path.read_bytes()).hexdigest()
    for name,new in [('render.json',r),('animation.json',a)]:
        p=output/name;data=json.loads(p.read_text(encoding='utf-8'));byid={v['id']:v for v in data['rows']};byid[new['id']]=new
        data['catalog_sha256']=digest;data['rows']=[byid[v['id']] for v in catalog['rows'] if v['id'] in byid]
        p.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(portrait=r['file'],animation=a['file'],frames=a['frames'],seconds=a['seconds_per_revolution'])),flush=True)


if __name__=='__main__':main()
