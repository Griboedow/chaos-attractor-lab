"""Exact affine-arithmetic audit and frozen-source replay of finite continuation."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
from fractions import Fraction as Q
from audit_product_orientation import endpoint,interval,point,add,sub,mul,isum,contains,mm,mv


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    read=lambda p:json.loads(p.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name in [Path(__file__).name,'audit_product_orientation.py']:
        (args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    prior=read(args.certificate);atlas=read(args.certificate.parent/'atlas.json');rows=[]
    K=list(map(interval,prior['root']['krawczyk']));eps=Q(prior['epsilon']);L=Q(prior['graph']['L'])
    for folder in args.runs:
        report=read(folder/'report.json');trail=read(folder/'trajectory.json')
        assert report['status']=='completed' and sha(folder/'trajectory.json')==report['trajectory_sha256']
        for row in trail:
            E=[[point(v) for v in r] for r in row['E']];I=[[point(v) for v in r] for r in row['inverse']]
            assert mm(E,I)==[[point(i==j) for j in range(3)] for i in range(3)]
            center=[point(endpoint(v)) for v in row['center_binary']];delta=list(map(interval,row['chart_error']))
            physical=mv(E,delta);box=list(map(interval,row['box']))
            assert all(contains(b,add(c,v)) for b,c,v in zip(box,center,physical))
            # Upper distance is independently checked without floating square roots.
            differences=[sub(b,point(Q(v))) for b,v in zip(box,row['nearest_state'])]
            squared_upper=sum(max(abs(a),abs(b))**2 for a,b in differences)
            distance=interval(row['distance']);assert distance[0]>=0 and distance[1]**2>=squared_upper
            k=row['iteration']//4;assert row['iteration']==4*k
            if k==0:
                original=mv(I,[sub(add(K[i],mul(E[i][0],point(eps))),center[i]) for i in range(3)])
                original=[original[0],add(original[1],(-L*eps,L*eps)),add(original[2],(-L*eps,L*eps))]
                assert all(contains(a,b) for a,b in zip(delta,original))
            if k+1<len(trail):
                step=row['step'];assert len(step['radii'])==4 and all(interval(v)[0]>0 for v in step['radii'])
                In=[[point(v) for v in r] for r in trail[k+1]['inverse']]
                J=[[interval(v) for v in r] for r in step['derivative']]
                T=[[interval(v) for v in r] for r in step['transition']]
                actual=mm(mm(In,J),E)
                assert all(contains(T[i][j],actual[i][j]) for i in range(3) for j in range(3))
                propagated=mv(T,delta);nxt=list(map(interval,trail[k+1]['chart_error']))
                assert all(contains(nxt[i],add(interval(step['center_error'][i]),propagated[i])) for i in range(3))
        escaped=next(r for r in trail if r['iteration']==80);gaps=[]
        for a in atlas:
            b=list(map(interval,escaped['box']));c=list(map(interval,a['box']))
            gaps.append(max(max(b[i][0]-c[i][1],c[i][0]-b[i][1]) for i in range(3)))
        assert min(gaps)>Q('0.00011729734758')
        last=trail[-1];assert last['iteration']==200 and interval(last['distance'])[1]<Q('0.000011199')
        dest=args.out/'replays'/folder.name
        command=[sys.executable,str(folder/'continue_product_unstable_graph.py'),'--input',str(args.input),
                 '--certificate',str(args.certificate),'--source',str(args.source),'--mpmath-path',str(args.mpmath_path),
                 '--engine',report['engine'],'--digits',str(report['digits']),'--steps','200','--out',str(dest)]
        subprocess.run(command,check=True,capture_output=True,text=True)
        assert sha(dest/'trajectory.json')==sha(folder/'trajectory.json')
        rows.append(dict(run=folder.name,report_sha256=sha(folder/'report.json'),steps=200,
                         minimum_escape_gap=float(min(gaps)),endpoint_distance_upper=float(interval(last['distance'])[1]),
                         exact_chart_and_distance_audit=True,exact_trajectory_replay=True))
    result=dict(status='passed',runs=rows,scope='Exact rational audit of moving charts, initial graph-point enclosure, propagated errors, atlas escape and distance upper bounds; frozen-source replay of both interval expressions. The prior graph theorem is imported, and the reference state is numerical data, not a certified invariant-set point.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['input','certificate','source','mpmath-path','out']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--runs',type=Path,nargs=2,required=True);main(p.parse_args())
