"""Targeted comparison with Bischi--Tramontana (2010), Figures 3 and 8.

Finite clouds at five published parameter values. Local alignment does not
establish conjugacy, topology, or historical novelty.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numba import njit
from compare_attractor_shapes import align,normalize,fingerprint,distance
from compare_affine_shapes import whiten


@njit(cache=True)
def simulate(a,warm,count):
    A=np.array([[-1.,a,-.5],[-.5,-1.,a],[a,-.5,-1.]])
    x=np.array([.2,.3,.4]);w=np.array([1.,.371,-.217]);w/=np.linalg.norm(w)
    tail=np.empty((100000,3));sums=np.zeros(10);minimum=1e100;maximum=-1e100
    for it in range(warm+count):
        linear=2+A@x;y=x*linear;J=np.diag(linear)+np.diag(x)@A
        assert np.isfinite(y).all() and y.min()>0 and abs(y).max()<10
        w=J@w;length=np.linalg.norm(w);w/=length;x=y
        if it>=warm:sums[(it-warm)//(count//10)]+=np.log(length)/(count//10)
        if it>=warm+count-100000:tail[it-(warm+count-100000)]=x
        minimum=min(minimum,x.min());maximum=max(maximum,x.max())
    return tail,sums,minimum,maximum


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for name in [Path(__file__).name,'compare_attractor_shapes.py','compare_affine_shapes.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    with np.load(args.source) as f:z=f['states']
    parameters=[.65,.6575,.665,.66824,.675];refs={};numerics=[]
    for a in parameters:
        cloud,rates,low,high=simulate(a,100000,1000000);refs[a]=cloud
        np.savez_compressed(args.out/f'lv-{a}.npz',states=cloud,rates=rates,params=[2.,-1.,a,-.5],initial=[.2,.3,.4])
        numerics.append(dict(a=a,leading_blocks=rates.tolist(),minimum_coordinate=low,maximum_coordinate=high))
    rng=np.random.default_rng(2026110158);comparisons={};aligned={}
    for name,transform in [('similarity',normalize),('affine',whiten)]:
        q=transform(z);qa=transform(z[:len(z)//2]);qb=transform(z[len(z)//2:]);ia=rng.choice(len(qa),1500,False);ib=rng.choice(len(qb),1500,False)
        self_rms,_=align(qa[ia],qb[ib]);rows=[]
        for a,cloud in refs.items():
            ref=transform(cloud);ai=rng.choice(len(q),1500,False);bi=rng.choice(len(ref),1500,False)
            rms,fit=align(q[ai],ref[bi]);fp=float(distance(fingerprint(q),fingerprint(ref)))
            rows.append(dict(a=a,local_alignment_rms=rms,fingerprint_distance=fp))
            key=f'{name}-{a}';aligned[key+'-candidate']=fit;aligned[key+'-reference']=normalize(ref[bi]);aligned[key+'-candidate-indices']=ai;aligned[key+'-reference-indices']=bi
        comparisons[name]=dict(self_time_halves_rms=self_rms,references=rows)
    np.savez_compressed(args.out/'alignments.npz',**aligned)
    # Algebraic Jacobian cross-check at declared interior probes.
    A=np.array([[-1.,.675,-.5],[-.5,-1.,.675],[.675,-.5,-1.]])
    errors=[]
    for x in refs[.675][::3125]:
        f=lambda q:q*(2+A@q);J=np.diag(2+A@x)+np.diag(x)@A;h=1e-7
        fd=np.column_stack([(f(x+h*np.eye(3)[j])-f(x-h*np.eye(3)[j]))/(2*h) for j in range(3)]);errors.append(float(np.max(abs(J-fd))))
    assert max(errors)<1e-7
    report=dict(status='completed',source_sha256=sha(args.source),paper='https://doi.org/10.1016/j.cnsns.2009.10.021',
                author_pdf='https://www.mdef.it/public/fileadmin/whoiswho/bischi/Selected_publications/Research_papers/Bi-Tra-Commun_Nonlin_Sci-2010.pdf',
                model='x_i_next=x_i*(2+(A*x)_i); A=[[-1,a,-0.5],[-0.5,-1,a],[a,-0.5,-1]]',parameters=parameters,initial=[.2,.3,.4],warm=100000,observed=1000000,retained_last=100000,
                seed=2026110158,numerics=numerics,comparisons=comparisons,maximum_jacobian_fd_error=max(errors),seconds=time.time()-begin,
                protocol='Equation (1), Figure 3 parameters, Figure 8(b,c) at a=0.675. Initial point declared here; not specified by the inspected figure caption. 1500-point local ICP, 48 PCA starts, 20 iterations, reflections allowed.',scope=__doc__)
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(dict(comparisons=comparisons,seconds=report['seconds'])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
