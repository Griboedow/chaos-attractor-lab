"""Cross-link the covering, cycle and presentation evidence without resampling."""
import argparse,hashlib,json,sys
from fractions import Fraction as Q
from pathlib import Path
import numpy as np

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bounds(s):return list(map(Q,s.strip('[]').split(',')))
def disjoint(a,b):return any(x[1]<y[0] or y[1]<x[0] for x,y in zip(a,b))

def main(args):
    root=Path(__file__).resolve().parents[1];args.out.mkdir(parents=True,exist_ok=False);reports={};files={}
    names=['soft-cycles01','soft-returns01','soft-covering-proposals01','soft-covering-certificate01','soft-covering-audit01','soft-covering-ad01','soft-covering-ad-audit01','soft-cycle-certificates01','soft-cycle-ad01','soft-symbolic-render01','soft-symbolic-render02']
    for name in names:
        folder=args.base/name;r=json.loads((folder/'report.json').read_text());reports[name]=r;files[name+'/report.json']=sha(folder/'report.json')
        for source,digest in r.get('sources',{}).items():assert sha(folder/source)==digest
    assert len(reports['soft-cycles01']['attempts'])==640 and len(reports['soft-cycles01']['cycles'])==65 and len(reports['soft-covering-proposals01']['proposals'])==39
    first=reports['soft-covering-certificate01'];second=reports['soft-covering-ad01'];assert first['status']==second['status']=='passed' and first['input_sha256']==second['input_sha256'];assert second['primary_report_sha256']==files['soft-covering-certificate01/report.json']
    for cert,audit,total in [('soft-covering-certificate01','soft-covering-audit01',4737),('soft-covering-ad01','soft-covering-ad-audit01',9918)]:
        r=reports[audit];assert r['status']=='passed' and r['complete_graph'] and r['window0_in_trap'];assert r['total_leaves']==total and r['source_report_sha256']==files[cert+'/report.json']
    first=reports['soft-cycle-certificates01'];second=reports['soft-cycle-ad01'];assert first['status']==second['status']=='passed' and second['primary_report_sha256']==files['soft-cycle-certificates01/report.json'];assert first['trap_sha256']==second['trap_sha256']==sha(args.trap)
    all_orbits=[];rows=[]
    for original,control in zip(first['cases'],second['cases']):
        assert original['id']==control['id'] and original['period']==control['period'] and original['unstable_dimensions']==control['unstable_dimensions'];ident=original['id'];assert sha(args.base/'soft-cycles01'/f'cycle-{ident}.npz')==original['source_sha256']
        for row in [original,control]:
            box=list(map(bounds,row['box']));K=list(map(bounds,row['krawczyk']));assert all(a[0]<b[0]<=b[1]<a[1] for a,b in zip(box,K));assert bounds(row['contraction_bound'])[1]<1;roots=[];unstable=0
            for bracket in row['multiplier_brackets']:
                x=bounds(bracket['interval']);a=bounds(bracket['at_left']);b=bounds(bracket['at_right']);assert a[1]<0<b[0] or b[1]<0<a[0];expanding=x[1]<-1 or x[0]>1;assert expanding or -1<x[0]<=x[1]<1;assert expanding==bracket['expanding'];unstable+=expanding;roots.append(x)
            assert all(roots[j][1]<roots[j+1][0] for j in range(2)) and unstable==row['unstable_dimensions'];assert {v['divisor'] for v in row['proper_divisor_exclusions']}=={d for d in range(1,row['period']) if row['period']%d==0}
            for v in row['proper_divisor_exclusions']:assert any(a>0 or b<0 for a,b in map(bounds,v['difference']))
            assert len(row['trapping_membership'])==row['period']
        larger=list(map(bounds,control['box']));smaller=list(map(bounds,original['box']));assert all(a[0]<b[0]<=b[1]<a[1] for a,b in zip(larger,smaller))
        orbit=[list(map(bounds,v['box'])) for v in control['trapping_membership']]
        for previous in all_orbits:assert all(disjoint(a,b) for a in orbit for b in previous)
        all_orbits.append(orbit);rows.append(dict(id=ident,period=original['period'],unstable_dimensions=original['unstable_dimensions']))
    render=reports['soft-symbolic-render02'];assert render['states_count']==render['histogram_count']==2000000 and render['spatial_clipping']==0
    for name,digest in render['figures'].items():assert sha(root/'docs/soft-symbolic-dynamics/figures'/name)==digest
    with np.load(args.base/'soft-symbolic-render02/display-coordinates.npz') as f:
        for row in first['cases']:
            expected=np.array([[sum(map(float,v.strip('[]').split(',')))/2 for v in q['box']] for q in row['trapping_membership']]);assert np.array_equal(f[f"orbit_{row['id']}"],expected)
    result=dict(status='passed',report_hashes=files,all_source_snapshots_match=True,declared_search_attempts=640,numerical_cycles=65,window_proposals=39,complete_covering_leaves=[4737,9918],same_periodic_roots_independently_certified=rows,all_three_orbits_pairwise_disjoint=True,final_figure_hashes=render['figures'],scope='Cross-packet identities and exact saved-inequality checks, distinct periodic orbits and display provenance. Complete interval recomputation is supplied by the two separate engines and their audits, not repeated here.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');(args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps({k:v for k,v in result.items() if k!='report_hashes'}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--trap',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
