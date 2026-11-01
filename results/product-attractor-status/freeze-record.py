"""Freeze the declared Case 988 status runs without altering earlier packets."""
import hashlib,json,shutil,sys
from pathlib import Path
import numpy as np
root=Path('D:/Codex/chaos-attractor-research-20261101');raw=Path('D:/Codex/chaos-attractor-search-20261031')
packet=root/'results/product-attractor-status';cache=Path('D:/Codex/chaos-attractor-cache')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
def write(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
registry=read(root/'records/generated-data.json');assert not any(r.get('group')=='product-attractor-status' for r in registry['files'])
omissions=[];inventory=[];ignores=[];large=[]
for run in sorted(raw.glob('product-attractor988-*')):
    assert run.is_dir()
    for source in sorted(run.rglob('*')):
        if not source.is_file() or '__pycache__' in source.parts:continue
        relative=source.relative_to(run);dest=packet/'runs'/run.name/relative;size=source.stat().st_size;digest=sha(source)
        entry=dict(run=run.name,file=relative.as_posix(),bytes=size,sha256=digest,storage='retained')
        if source.suffix=='.npz' and size>=262144:
            entry['storage']='optional cache';cachefile=cache/digest[:2]/(digest+'.npz');cachefile.parent.mkdir(parents=True,exist_ok=True)
            if not cachefile.exists():shutil.copyfile(source,cachefile)
            assert sha(cachefile)==digest
            arrays=[];keep={}
            with np.load(source,allow_pickle=False) as data:
                for k in data.files:
                    a=data[k];retained=a.nbytes<131072 and k not in ['tail','attempts','occupied','pending','new']
                    arrays.append(dict(key=k,shape=list(a.shape),dtype=str(a.dtype),bytes=a.nbytes,retained=retained))
                    if retained:keep[k]=a
            compact=None
            if keep:
                cp=root/'results/compact-data'/(digest+'.npz');np.savez_compressed(cp,**keep)
                compact=dict(path=cp.relative_to(root).as_posix(),bytes=cp.stat().st_size,sha256=sha(cp),keys=list(keep))
                with np.load(cp) as saved,np.load(source) as original:
                    for k in keep:assert saved[k].dtype==original[k].dtype and saved[k].shape==original[k].shape and saved[k].tobytes()==original[k].tobytes()
            row=dict(path=dest.relative_to(root).as_posix(),bytes=size,git_blob=None,sha256=digest,group='product-attractor-status',compact=compact,arrays=arrays,
                recovery='recompute',recompute_command=f'See results/product-attractor-status/README.md for run {run.name}; use its retained source in a fresh external directory.',original_run_path=str(source))
            omissions.append(row);ignores.append('/'+row['path'])
        else:
            dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest);assert sha(dest)==digest
            if source.suffix=='.json' and size>=262144:
                reason='Full preimage-node selection and parent graph: essential frozen certificate replay input, not a thinned trajectory.' if source.name=='nodes.json' else 'Exact interval certificate or full selection record, retained as directly linked mathematical evidence.'
                large.append(dict(path=dest.relative_to(root).as_posix(),bytes=size,reason=reason))
        inventory.append(entry)
registry['files'].extend(omissions);registry['retained_large_inputs'].extend(large)
registry['removed_bytes']=sum(r['bytes'] for r in registry['files'])
compacts={r['compact']['path']:r['compact'] for r in registry['files'] if r['compact']}
registry['compact_bytes']=sum(r['bytes'] for r in compacts.values());write(root/'records/generated-data.json',registry)
with (root/'.gitignore').open('a',encoding='utf-8',newline='\n') as f:f.write('\n'+'\n'.join(ignores)+'\n')
write(packet/'inventory.json',dict(schema=1,scope='All scientific files from every product-attractor988-* run, including failed audits, all proposals, two figures and nested replays. Python bytecode caches excluded.',files=inventory))
write(packet/'omissions.json',dict(schema=1,policy='Full arrays remain external optional payloads; compact arrays are explicitly partial and never replace full archives.',files=omissions))
shutil.copyfile(raw/'product-alphabet-dense01/case-988.json',packet/'reference-camera.json')
shutil.copyfile(__file__,packet/'freeze-record.py')
accepted=read(raw/'product-attractor988-chain-audit04/report.json');verification=read(raw/'product-attractor988-fate-verification01/report.json')
assert accepted['status']==verification['status']=='passed'
write(packet/'validation.json',dict(status='passed',exact_audit='runs/product-attractor988-chain-audit04/report.json',exact_audit_sha256=sha(raw/'product-attractor988-chain-audit04/report.json'),finite_verification='runs/product-attractor988-fate-verification01/report.json',finite_verification_sha256=sha(raw/'product-attractor988-fate-verification01/report.json'),
    figure_sha256=sha(root/'docs/assets/diagnostics/product-alphabet-attractor-status.png'),figure_visual_review='Both versions inspected; final render02 has readable labels, complete domains, improved cloud contrast and no annotation/boundary overlap.',
    numerical_scope=verification['scope'],certificate_scope=accepted['scope'],all_raw_hashes_verified=True,all_new_compact_arrays_bit_checked=True))
write(packet/'literature-review.json',dict(date='2026-11-02',scope='Definition and method context only; no new priority, conjugacy or attractor theorem inferred from these sources.',sources=[
    dict(title='Milnor (1985), On the concept of attractor',url='https://link.springer.com/article/10.1007/BF01212280',read_scope='Publisher abstract and bibliographic page. Full Project Euclid text was inaccessible behind a security check.',use='The probabilistic/asymptotic view of attraction is established. No existence theorem transferred to this singular map.'),
    dict(title='Dellnitz and Junge (1998), An adaptive subdivision technique for the approximation of attractors and invariant measures',url='https://portal.fis.tum.de/en/publications/an-adaptive-subdivision-technique-for-the-approximation-of-attrac-2/',read_scope='Institutional abstract and metadata; no theorem or figure audit.',use='Set-oriented subdivision is an established method, not a claimed innovation here.')]))
write(packet/'accounting.json',dict(runs=len({r['run'] for r in inventory}),raw_files=len(inventory),retained_raw_files=sum(r['storage']=='retained' for r in inventory),optional_arrays=len(omissions),optional_bytes=sum(r['bytes'] for r in omissions),all_raw_hashes_verified=True))
write(raw/'product-attractor-status-pages.json',['README.md','docs/product-alphabet/manuscript.md','docs/reproduction.md','results/README.md','results/product-attractor-status/README.md'])
print(json.dumps(read(packet/'accounting.json')))
