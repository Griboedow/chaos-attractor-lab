"""Freeze compact numerical records and a portable cloud registry for the extension."""
import argparse
import csv
import gzip
import json
from pathlib import Path


def read(p):
    return json.loads(p.read_text())


def write(p,v):
    data=(json.dumps(v,indent=2,allow_nan=False)+'\n').encode()
    p.write_bytes(gzip.compress(data,mtime=0) if p.suffix=='.gz' else data)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();b=a.cache;r=a.out;r.mkdir(exist_ok=True,parents=True)
    d=read(b/'events-a/data.json');write(r/'bifurcations.json',d['events'])
    write(r/'orbit-controls.json',{p.stem:read(p) for p in (b/'sections-a').glob('*.json')})
    if (b/'basins-a/basins.json').exists():write(r/'basins.json',read(b/'basins-a/basins.json'))
    cases=[dict(v,kind='periodic_orbit',point_dtype='float64',label=k.replace('_',' '),recipe='rotated_nodes')
           for k,v in d['representatives'].items()]
    co=b/'coexistence-a/report.json'
    if co.exists():
        info=read(co);write(r/'chaotic-replays.json',info)
        if 'cloud' in info:cases.append(dict(info['cloud'],case=180844,family='ring',identifier='coexisting_chaotic_cloud',
                                      e=info['e'],kind='finite_time_chaotic_cloud',recipe='rk4',label='Long-lived positive-growth trajectory at the stable-loop coefficient'))
    if (b/'long-fates-a/report.json').exists():write(r/'long-fates.json',read(b/'long-fates-a/report.json'))
    for source,file,tag in [('daughters-a/results.json','postflip.json','full_dop853'),('large-period-a/report.json','period-32.json','single_relative_sector'),('period64-a/report.json','period-64.json','single_relative_sector')]:
        if (b/source).exists():
            info=read(b/source);write(r/file,info)
            for row in (info if isinstance(info,list) else info['postflip']):
                if 'cloud_relative' in row:
                    name='postflip_'+str(row['order_relative_to_2T']) if isinstance(info,list) else ('period_64' if '64' in file else 'period_32')
                    cases.append(dict(row,case=180844,family='ring',kind='periodic_orbit',label=name.replace('_',' '),
                                      identifier=name,recipe=tag,point_dtype='float64'))
    write(r/'manifest.json',dict(campaign_cache_root=b.parent.parent.as_posix(),cases=cases))
    fields=['branch','direction','i','arclength','parameter_fraction','e','period','multiplier','neutral','relative_node_spread','de_ds','residual']
    with (r/'continuation.csv').open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=fields);writer.writeheader()
        for folder in ['asym-down','asym-up','sym-down']:
            for line in (b/folder/'branch.jsonl').read_text().splitlines():
                rec=json.loads(line);row={k:rec[k] for k in fields if k!='direction'};row.update(direction=folder,multiplier=rec['multiplier'][0])
                writer.writerow(row)
    print('Portable cloud records',len(cases))


if __name__=='__main__':
    main()
