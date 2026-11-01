"""Source accounting, reverse-order replay, and independent interval entry."""
import argparse,hashlib,json,sys,time
from fractions import Fraction
from pathlib import Path
import numpy as np
from audit_soft_trapping import audit
from propose_raar_trapping_bitset import mark


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);begin=time.time()
    if args.mpmath_path:sys.path.append(str(args.mpmath_path))
    import mpmath as mp
    mp.mp.dps=100;mp.iv.dps=80
    record=json.loads((args.audit/'report.json').read_text(encoding='utf-8'))
    assert record['status']=='all_cubes_passed' and record['source_sha256']==sha(args.source)
    for name,digest in record['sources'].items():
        assert sha(args.audit/name)==digest and sha(Path(__file__).with_name(name))==digest,name
    assert record['processed']==record['cubes'] and record['leaves']==record['cubes']+7*record['splits']
    offset=0
    for row in record['history']:
        assert row['offset']==offset and row['processed']==row['count'] and row['status']==0
        assert row['leaves']==row['count']+7*row['splits'];offset+=row['count']
    assert offset==record['cubes']
    with np.load(args.source) as f:beta,tau=f['params'];points=f['points'];origin=f['origin'];shape=f['shape'];h=float(f['step']);occupied=f['occupied']
    with np.load(args.reference) as f:
        assert np.array_equal(f['params'],[beta,tau]) and np.array_equal(f['points'],points)
        z=f['states'];indices=np.floor((z-origin)/h).astype(np.int64);keys=(indices[:,0]*shape[1]+indices[:,1])*shape[2]+indices[:,2]
    assert h==2.**-8 and record['margin']==2.**-35 and np.all(np.diff(occupied)>0)
    assert occupied[0]>=0 and occupied[-1]<np.prod(shape)
    for x in origin:assert (Fraction(float(x))/Fraction(h)).denominator==1
    expected=set(occupied.tolist());assert set(keys.tolist())<=expected
    mask=np.zeros((int(np.prod(shape))+7)//8,np.uint8)
    for key in occupied:mark(mask,key)
    decoded=np.flatnonzero(np.unpackbits(mask,bitorder='little'));assert np.array_equal(decoded,occupied)
    leaves=splits=processed=0;axis=1e100;depth=0
    reverse=occupied[::-1].copy()
    for offset in range(0,len(reverse),750):
        block=reverse[offset:offset+750]
        flag,pos,key,missing,nleaf,nsplit,r,d,error=audit(block,mask,shape,origin,h,beta,tau,points,2.**-35,8,2000000)
        assert flag==0 and pos==len(block),(flag,key,missing)
        assert nleaf==len(block)+7*nsplit
        processed+=pos;leaves+=nleaf;splits+=nsplit;axis=min(axis,r);depth=max(depth,d)
        print(json.dumps(dict(processed=processed,seconds=time.time()-begin)),flush=True)
    assert (processed,leaves,splits,axis,depth)==(record['processed'],record['leaves'],record['splits'],record['minimum_axis_lower_bound'],record['maximum_depth'])
    (args.out/'replay.json').write_text(json.dumps(dict(status='passed',processed=processed,leaves=leaves,splits=splits,axis=axis,depth=depth),indent=2)+'\n',encoding='utf-8',newline='\n')
    # Natural interval wrapping needs a much smaller box than the separate
    # binary64 mean-value entry certificate, which retains radius 2^-35.
    # This independent library control uses radius 2^-160.
    def scalar(x):return mp.iv.mpf([mp.mpf(float(x)),mp.mpf(float(x))])
    b=scalar(beta);temperature=scalar(tau);sites=[[scalar(a) for a in q] for q in points]
    start=np.array([.317,-.211,.419]);delta=mp.iv.mpf([-mp.mpf(2)**-160,mp.mpf(2)**-160]);x=[scalar(a)+delta for a in start]
    for k in range(11):
        logits=[(2*sum(q[j]*x[j] for j in range(3))-sum(a*a for a in q))/temperature for q in sites]
        weights=[mp.iv.exp(a) for a in logits];total=sum(weights);weights=[a/total for a in weights]
        m=[sum(weights[t]*sites[t][j] for t in range(len(sites))) for j in range(3)]
        rb=[2*m[j]-x[j] for j in range(3)];r=mp.iv.sqrt(rb[0]**2+rb[1]**2);a=[rb[0]/r,rb[1]/r,scalar(0)]
        x=[b/2*(2*a[j]-rb[j]+x[j])+(1-b)*m[j] for j in range(3)]
    lo=[Fraction.from_float(float(np.nextafter(float(mp.mpf(v.a)),-np.inf))) for v in x]
    hi=[Fraction.from_float(float(np.nextafter(float(mp.mpf(v.b)),np.inf))) for v in x]
    for i in range(3):assert mp.mpf(float(lo[i]))<=mp.mpf(x[i].a)<=mp.mpf(x[i].b)<=mp.mpf(float(hi[i]))
    lower=[int((lo[i]-Fraction(float(origin[i])))//Fraction(h)) for i in range(3)]
    upper=[int((hi[i]-Fraction(float(origin[i])))//Fraction(h)) for i in range(3)]
    assert all(0<=lower[i]<=upper[i]<shape[i] for i in range(3));tested=0
    for ix in range(lower[0],upper[0]+1):
        for iy in range(lower[1],upper[1]+1):
            for iz in range(lower[2],upper[2]+1):assert int((ix*shape[1]+iy)*shape[2]+iz) in expected;tested+=1
    result=dict(status='passed',source_sha256=sha(args.source),audit_report_sha256=sha(args.audit/'report.json'),reference_sha256=sha(args.reference),cubes=processed,leaves=leaves,splits=splits,minimum_axis_lower_bound=axis,maximum_depth=depth,reverse_order_replay='all aggregate results identical with reversed source order and 750-cube batches',mask_membership='independently decoded with numpy.unpackbits; identical sorted addresses',terminal_reference_states=len(z),reference_coverage='every terminal reference state lies in the proposed union',independent_entry=dict(library='mpmath 1.3.0 interval context',digits=80,initial=start.tolist(),initial_coordinate_radius=2.**-160,iterations=11,final_box=[[float(lo[i]),float(hi[i])] for i in range(3)],intersected_cubes=tested,status='contained'),seconds=time.time()-begin,scope='Trapping source/accounting/replay and independent finite entry. No proof of chaos, identity of the rendered cloud with the maximal invariant set, or historical uniqueness.')
    (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    (args.out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--audit',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--mpmath-path',type=Path);main(p.parse_args())
