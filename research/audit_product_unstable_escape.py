"""Exact rational audit and frozen-source replay of unstable-graph escapes."""
import argparse, hashlib, json, subprocess, sys
from fractions import Fraction as Q
from pathlib import Path
from audit_product_orientation import endpoint, interval, point, add, sub, mul, isum, contains, mm, mv, det


def read(p): return json.loads(p.read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(v): return max(abs(x) for x in v)


def audit(folder, windows, cycles):
    r=read(folder/'report.json');assert r['status']=='passed'
    assert r['windows_sha256']==sha(windows) and r['cycles_sha256']==sha(cycles)
    for name,h in r['sources'].items():assert sha(folder/name)==h
    root=r['root']; x=list(map(point,root['center_decimal'])); radius=Q(root['radius'])
    B=list(map(interval,root['box'])); K=list(map(interval,root['krawczyk']))
    C=[[point(v) for v in row] for row in root['preconditioner']]
    H=[[interval(v) for v in row] for row in root['derivative']]
    R=[[interval(v) for v in row] for row in root['remainder']]
    CH=mm(C,H);res=list(map(interval,root['residual']))
    correction=mv(C,res)
    for i in range(3):
        assert contains(B[i],add(x[i],(-radius,radius)))
        for j in range(3):assert contains(R[i][j],sub(point(i==j),CH[i][j]))
        bound=add(sub(x[i],correction[i]),isum(mul(R[i][j],(-radius,radius)) for j in range(3)))
        assert contains(K[i],bound) and B[i][0]<K[i][0]<=K[i][1]<B[i][1]
    assert max(sum(norm(v) for v in row) for row in R)<1
    assert det(C)[0]!=0 and det(C)[0]==det(C)[1]
    g=r['graph'];E=[[point(v) for v in row] for row in g['E']];I=[[point(v) for v in row] for row in g['inverse']]
    assert mm(E,I)==[[point(i==j) for j in range(3)] for i in range(3)]
    graph_radius,L=Q(g['radius']),Q(g['L']);box=list(map(interval,g['physical_box']))
    offsets=mv(E,[(-graph_radius,graph_radius),(-graph_radius*L,graph_radius*L),(-graph_radius*L,graph_radius*L)])
    for b,k,o in zip(box,K,offsets):assert contains(b,add(k,o))
    J=[[interval(v) for v in row] for row in g['jacobian']]
    assert J[0][0][1]<0 or J[0][0][0]>0
    a=interval(g['a'])[0];b=interval(g['b'])[1];c=interval(g['c'])[1];d=interval(g['d'])[1]
    assert a<=min(abs(v) for v in J[0][0]) and b>=sum(norm(v) for v in J[0][1:])
    assert c>=max(norm(J[i][0]) for i in [1,2]) and d>=max(sum(norm(v) for v in J[i][1:]) for i in [1,2])
    expansion=a-b*L
    assert expansion>1 and c+d*L<L*expansion and (1-d)*expansion>(c+d*L)*b and c+d*L<L
    trail=read(folder/'trajectory.json');atlas=read(folder/'atlas.json')
    assert sha(folder/'trajectory.json')==r['trajectory_sha256'] and sha(folder/'atlas.json')==r['atlas_sha256']
    assert len(atlas)==92 and {(v['window'],v['phase']) for v in atlas}=={(i,j) for i in range(23) for j in range(4)}
    raw=read(windows)
    unique=r['selfloop_uniqueness']
    tree={};volume=Q(0)
    for cell in unique['cells']:
        node=tree;mid=[Q(0)]*3;rad=[Q(1)]*3;address=cell['address'];assert len(address)%2==0
        for k in range(0,len(address),2):
            axis=int(address[k]);side=int(address[k+1]);assert axis in [0,1,2] and side in [0,1] and '#' not in node
            rad[axis]/=2;mid[axis]+=(2*side-1)*rad[axis];node=node.setdefault(address[k:k+2],{})
        assert not node;node['#']=True;volume+=Q(1,2**(len(address)//2))
        spread=mv([[point(v) for v in row] for row in raw['exact_windows'][0]['matrix']],[(m-h,m+h) for m,h in zip(mid,rad)])
        for j in range(3):assert contains(interval(cell['outer_box'][j]),add(point(Q(raw['centers'][0][j])),spread[j]))
        H0=[[interval(v) for v in row] for row in cell['derivative']]
        R0=[[interval(v) for v in row] for row in cell['remainder']]
        CH0=mm(C,H0)
        for i in range(3):
            for j in range(3):assert contains(R0[i][j],sub(point(i==j),CH0[i][j]))
        assert max(sum(norm(v) for v in row) for row in R0)<1
    def complete(node):
        if '#' in node:assert set(node)=={'#'};return
        keys=sorted(node);assert len(keys)==2 and keys[0][0]==keys[1][0] and keys[0][1]=='0' and keys[1][1]=='1'
        for child in node.values():complete(child)
    complete(tree);assert volume==1
    local=mv([[point(v) for v in row] for row in raw['exact_windows'][0]['inverse']],
             [sub(K[i],point(Q(raw['centers'][0][i]))) for i in range(3)])
    for stored,value in zip(unique['root_location'],local):assert contains(interval(stored),value) and -1<value[0]<=value[1]<1
    for i,w in enumerate(raw['exact_windows']):
        spread=mv([[point(v) for v in row] for row in w['matrix']],[(-Q(1),Q(1))]*3)
        enclosure=list(map(interval,next(v for v in atlas if v['window']==i and v['phase']==0)['box']))
        for j in range(3):assert contains(enclosure[j],add(point(Q(raw['centers'][i][j])),spread[j]))
    eps=Q(r['epsilon']);assert 0<eps<graph_radius
    initial=mv(E,[point(eps),(-eps*L,eps*L),(-eps*L,eps*L)])
    initial=[add(k,o) for k,o in zip(K,initial)]
    nominal0=[point(endpoint(t)) for t in trail[0]['center_binary']]
    error=mv(I,[sub(v,n) for v,n in zip(initial,nominal0)])
    for stored,needed in zip(trail[0]['chart_error'],error):assert contains(interval(stored),needed)
    for k,t in enumerate(trail):
        assert t['iteration']==4*k
        current=[point(endpoint(v)) for v in t['center_binary']];delta=list(map(interval,t['chart_error']))
        physical=mv(E,delta)
        for j in range(3):assert contains(interval(t['box'][j]),add(current[j],physical[j]))
        if k+1<len(trail):
            step=t['step'];assert len(step['radii'])==4 and all(interval(v)[0]>0 for v in step['radii'])
            prop=mv([[interval(v) for v in row] for row in step['chart_derivative']],delta)
            for j in range(3):assert contains(interval(trail[k+1]['chart_error'][j]),add(interval(step['center_error'][j]),prop[j]))
    final=r['separated'];assert final['iteration']==trail[-1]['iteration'] and final['box']==trail[-1]['box']
    gaps=[]
    for check,ab in zip(final['separations'],atlas):
        j=check['axis'];endbox=interval(final['box'][j]);targetbox=interval(ab['box'][j]);gap=sub(endbox,targetbox) if check['sign']==1 else sub(targetbox,endbox)
        # The record stores the endpoint separation as a degenerate outward interval.
        assert interval(check['gap'])[0]<=gap[0]<=interval(check['gap'])[1] and gap[0]>Q('1e-4')
        gaps.append(gap[0])
    return dict(run=folder.name,report_sha256=sha(folder/'report.json'),atlas_boxes=92,uniqueness_partition_cells=len(unique['cells']),four_step_blocks=len(trail)-1,
                minimum_separation=float(min(gaps)),graph_expansion=float(expansion),graph_contraction=float(d+(c+d*L)*b/expansion))


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);rows=[];replays=[]
    for name in [Path(__file__).name,'audit_product_orientation.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    for folder in args.runs:
        rows.append(audit(folder,args.windows,args.cycles));r=read(folder/'report.json')
        dest=args.replays/folder.name
        cmd=[sys.executable,str(folder/'certify_product_unstable_escape.py'),'--windows',str(args.windows),'--cycles',str(args.cycles),'--mpmath-path',str(args.mpmath_path),'--out',str(dest),'--engine',r['engine'],'--digits',str(r['digits']),'--root-radius',r['root']['radius'],'--graph-radius',r['graph']['radius'],'--slope',r['graph']['L'],'--epsilon',r['epsilon']]
        subprocess.run(cmd,check=True,capture_output=True,text=True)
        replay=read(dest/'report.json')
        for key in ['root','graph','selfloop_uniqueness','separated','atlas_sha256','trajectory_sha256']:assert replay[key]==r[key]
        replays.append(dict(run=folder.name,exact_payloads=True))
    for name in [Path(__file__).name,'audit_product_orientation.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    report=dict(status='passed',certificates=rows,replays=replays,
                scope='Exact rational checks of Krawczyk arithmetic, the unstable-graph inequalities, chart propagation and separation from all 92 atlas boxes; full frozen-source replays. Transcendental enclosures share mpmath interval arithmetic and are checked through two posterior formulas.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['windows','cycles','out','replays','mpmath-path']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--runs',type=Path,nargs=2,required=True);main(p.parse_args())
