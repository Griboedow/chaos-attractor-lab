"""Compile case evidence, stage counts and checksums of the external checkpoints."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def read(path): return json.loads(path.read_text(encoding='utf8'))


def write(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8')


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main(args):
    stages=[];checkpoints=[];pair=np.zeros(4,dtype=int)
    for geom,prefix,select in [('single','','inspection'),('linked','linked-','linked-inspection'),('smooth','smooth-','smooth-select')]:
        counts=np.zeros(5,dtype=int)
        for path in sorted((args.cache/(prefix+'survey')).glob('batch-*.npz')):
            j=int(path.stem.split('-')[1])
            with np.load(path) as f:n=int(f['settings'][0]);fast=len(f['flags'])
            deep=args.cache/(prefix+'deep')/f'deep-{j:05d}.npz'
            two=args.cache/(prefix+'two-growth')/deep.name
            with np.load(deep) as f:submitted=len(f['all_ids']);survived=len(f['ids'])
            with np.load(two) as f:hyper=len(f['ids'])
            counts += [n,fast,submitted,survived,hyper]
            checkpoints.append(dict(geometry=geom,batch=j,screen_sha256=sha(path),deep_sha256=sha(deep),two_rate_sha256=sha(two),counts=[n,fast,submitted,survived,hyper]))
            if geom=='smooth':
                with np.load(path) as a,np.load(args.cache/'linked-survey'/path.name) as b:
                    h=np.min(b['metrics'][:4],axis=0)>.015;s=np.min(a['metrics'][:4],axis=0)>.015
                    pair+=np.bincount(2*h.astype(int)+s.astype(int),minlength=4)
        with np.load(args.cache/(select+'.npz')) as f:eligible=int(f['total_eligible'][0]);inspection=len(f['ids'])
        cpu=read(args.out/(prefix+'two-rate-audit.json'))
        stages.append(dict(geometry=geom,screened=int(counts[0]),short_positive=int(counts[1]),spatial_rescreened=int(counts[2]),deep_positive=int(counts[3]),shape_eligible=eligible,shape_inspected=inspection,gpu_two_rate=int(counts[4]),cpu_two_rate=cpu['cpu_passes']))
    write(args.out/'checkpoints.json',dict(cache_prefix='projection/',batch_size=4096,columns=['screened','short_positive','spatial_rescreened','deep_positive','gpu_two_rate'],checkpoints=checkpoints))
    manifest=[];rows=[]
    proof_audit=args.out/'snapback'/'audit.json'
    entropy={v['case']:v for v in read(proof_audit)['cases']} if proof_audit.exists() else {}
    for path in sorted(args.out.glob('case-*.json'),key=lambda f:int(f.stem.split('-')[1])):
        rec=read(path);case=rec['case']
        with np.load(rec['cloud']) as f:z=f['points'].astype(float)
        eig=np.linalg.eigvalsh(np.cov(z.T));eig/=eig.sum()
        near=args.out/'neighborhood'/f'near-{case}.json'
        hood=read(near)['rows'] if near.exists() else []
        prec=[int(p.stem.split('-')[-1]) for p in (args.out/'precision').glob(f'mpfr-{case}-*.json')]
        two=rec['cohort'][0].get('qr_windows') is not None
        both=sum(all(r[1]>.015 for r in v['qr_windows']) and v['period'] is None for v in rec['cohort']) if two else None
        cyc=args.out/'cycles'/f'cycles-{case}.json'
        sink=args.out/'sinks'/f'cycles-{case}.json'
        ncyc=len(read(cyc)['cycles']) if cyc.exists() else 0
        nsink=len(read(sink)['cycles']) if sink.exists() else 0
        record=dict(case=case,geometry=rec['family'],label=rec['label'],passes=rec['passes'],lambda1=rec['rates'][0],lambda2=rec['rates'][1],lambda3=rec['rates'][2],min_window1=rec['minimum_leading_window'],min_window2=rec['minimum_second_window'],positive_starts=rec['cohort_positive'],two_positive_starts=both,pca_smallest=float(eig[0]),pca_middle=float(eig[1]),pca_largest=float(eig[2]),mpfr_bits='/'.join(map(str,sorted(prec))),positive_neighbors=sum(v['min_leading']>.015 for v in hood),tested_neighbors=len(hood),certified_cycles=ncyc,certified_sinks=nsink)
        proof=entropy.get(case)
        record['entropy_lower_denominator']=proof['return_steps'] if proof else None
        rows.append(record)
        if rec['passes']:
            asset=Path(rec['cloud'])
            assert sha(asset)==rec['cloud_sha256']
            manifest.append(dict(case=case,family=rec['family'],label=rec['label'],study='base survey',entropy_certificate='snapback/'+proof['packet'] if proof else None,entropy_lower_denominator=record['entropy_lower_denominator'],cloud=str(asset),path=str(asset),cloud_relative='projection/final/'+asset.name,record=path.name,rates=rec['rates'],positive_starts=rec['cohort_positive'],total_starts=len(rec['cohort']),two_positive_starts=both,source='research/shape_search_20261121/projection/reproduce.py',points_saved=rec['points_saved'],burn=rec['burn'],observed=rec['windows']*rec['window_length'],cloud_sha256=rec['cloud_sha256'],mpfr_bits=sorted(prec),certified_cycles=ncyc,certified_sinks=nsink))
    with (args.out/'cases.csv').open('w',newline='',encoding='utf8') as stream:
        out=csv.DictWriter(stream,fieldnames=list(rows[0]));out.writeheader();out.writerows(rows)
    base_retained=len(manifest)
    extra=args.out/'three'/'manifest.json'
    if extra.exists():
        for item in read(extra):
            item['record']='three/'+item['record']
            if item.get('entropy_certificate'):item['entropy_certificate']='three/'+item['entropy_certificate']
            manifest.append(item)
    write(args.out/'manifest.json',manifest)
    write(args.out/'summary.json',dict(stages=stages,matched_linked_short_screen=dict(neither=int(pair[0]),only_blended=int(pair[1]),only_hard=int(pair[2]),both=int(pair[3])),long_tested=len(rows),retained=base_retained,additional_three_positive=len(manifest)-base_retained,aggregate_retained=len(manifest),retained_with_second_min_above_0015=sum(q['passes'] and q['min_window2']>.015 for q in rows),certified_cycles=sum(q['certified_cycles'] for q in rows),certified_sinks=sum(q['certified_sinks'] for q in rows)))
    print(len(rows),len(manifest),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
