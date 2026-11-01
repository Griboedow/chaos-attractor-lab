"""Independent full-map Krawczyk checks of every classified inverse-image branch."""
import argparse,hashlib,json,time
from pathlib import Path
import mpmath as mp
from certify_raar_cycles import exactfloat,iterate,eye,matmul,matvec
from certify_raar_multiplier_events import lift


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);mp.mp.dps=200;mp.iv.dps=200;iv=mp.iv;begin=time.time();source=json.loads(args.input.read_text())['cases'];sturm=json.loads(args.sturm.read_text())['cases'];records=[]
    for name in ['check_raar_preimages_interval.py','certify_raar_cycles.py','certify_raar_multiplier_events.py']:(args.out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
    (args.out/'input.json').write_bytes(args.input.read_bytes());(args.out/'sturm-report.json').write_bytes(args.sturm.read_bytes())
    for case in source:
        p=[exactfloat(v,mp.mp) for v in case['params']];target=[exactfloat(v,mp.mp) for v in case['target']];pi=[exactfloat(v,iv) for v in case['params']];ti=[exactfloat(v,iv) for v in case['target']];expected=next(r for r in sturm if r['label']==case['label']);valid=[r for r in expected['roots'] if r['valid_preimage']];assert len(valid)==len(case['preimages']);rows=[];matches=[]
        for seed in case['preimages']:
            x=[mp.mpf(v) for v in seed]
            for _ in range(20):
                y,J,_,_=iterate(x,p,1,mp.mp);res=mp.matrix([y[i]-target[i] for i in range(3)])
                if max(abs(v) for v in res)<mp.mpf('1e-160'):break
                delta=mp.lu_solve(mp.matrix(J),res);x=[x[i]-delta[i] for i in range(3)]
            else:raise ArithmeticError('Full-map Newton failed')
            C=mp.inverse(mp.matrix(J));ci=[[lift(C[i,j],iv) for j in range(3)] for i in range(3)];center=[lift(v,iv) for v in x];radius=iv.mpf('1e-100');box=[iv.mpf([(v-radius).a,(v+radius).b]) for v in center]
            _,Ji,_,axes=iterate(box,pi,1,iv);yc,_,_,_=iterate(center,pi,1,iv);res=[yc[i]-ti[i] for i in range(3)];I=eye(iv);CJ=matmul(ci,Ji);R=[[I[i][j]-CJ[i][j] for j in range(3)] for i in range(3)];cf=matvec(ci,res);tail=matvec(R,[box[i]-center[i] for i in range(3)]);K=[center[i]-cf[i]+tail[i] for i in range(3)];contraction=max(sum(abs(v).b for v in row) for row in R)
            assert all(b.a<k.a and k.b<b.b for b,k in zip(box,K)) and contraction<1 and all(v.a>0 for v in axes)
            found=[]
            for root in valid:
                intervals=[iv.mpf([v.strip() for v in s.strip('[]').split(',')]) for s in root['preimage_box']]
                if all(a.a<b.a and b.b<a.b for a,b in zip(intervals,box)):found.append(root['root'])
            assert len(found)==1;matches.extend(found);rows.append(dict(sturm_root=found[0],center_decimal=[mp.nstr(v,210) for v in x],box=[str(v) for v in box],krawczyk=[str(v) for v in K],contraction=str(contraction),axis_distances=[str(v) for v in axes],newton_residual=mp.nstr(max(abs(v) for v in res),210)))
        assert sorted(matches)==sorted(r['root'] for r in valid);record=dict(label=case['label'],number_of_preimages=len(rows),preimages=rows);records.append(record);(args.out/(case['label']+'.json')).write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(dict(label=case['label'],verified=len(rows))),flush=True)
    report=dict(status='passed',interval_digits=200,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),sturm_report_sha256=hashlib.sha256(args.sturm.read_bytes()).hexdigest(),total_preimages=sum(r['number_of_preimages'] for r in records),cases=records,seconds=time.time()-begin,scope='Independent nested-reflector map and derivative validate every proposed preimage and place it in the independently derived Sturm box. Completeness comes from the Sturm classification; no attractor-membership or historical-priority assertion.')
    (args.out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--sturm',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
