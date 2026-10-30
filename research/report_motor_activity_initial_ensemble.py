"""Check every declared initial condition against the independently found saddle."""
import argparse, datetime, hashlib, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment, minimize_scalar
import motor_rows
from report_motor_activity_gpu import row_summary

ROOT=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def rows(y):
    return np.stack((y[:8],y[8:16],np.r_[y[16:],-sum(y[16:])]),axis=1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new analysis folder')
    protocol=json.loads((args.input/'protocol.json').read_text());data=np.load(args.input/'outcomes.npz')
    baseline=np.load(ROOT/protocol['source_path']);reference=np.load(args.reference)
    checks=[]
    def check(label,passed,**details):
        checks.append(dict(check=label,passed=bool(passed),**details))
    check('all16complete',len(data['status'])==16 and np.all(data['status']==0))
    check('full10000duration',np.all(abs(data['elapsed_original']-10000.)<1e-8))
    check('original_source_hash',sha(ROOT/protocol['source_path'])==protocol['source_sha256'])
    for name,expected in protocol['source_snapshot_sha256'].items():
        check('executed_source_hash',sha(args.input/name)==expected,source=name)
    p=reference['parameters'];period=float(reference['time'][-1]);initial=reference['full_states'][0]
    check('reference_parameters',np.array_equal(data['parameters'],np.tile(p,(16,1))))
    ref=solve_ivp(lambda t,y:motor_rows.rhs(t,y,p,8),(0.,period),initial,method='DOP853',rtol=5e-13,atol=5e-15,max_step=.01,dense_output=True)
    if not ref.success:
        raise RuntimeError(ref.message)
    outcomes=[]
    for i,case in enumerate(protocol['cases']):
        regenerated,_=motor_rows.source_start(8,case['seed'])
        check('seed_reproduction',np.array_equal(regenerated,data['initial_states'][i]),case=i)
        if i<3:
            check('original_initial_state',np.array_equal(data['initial_states'][i],baseline['initial_states'][45+i]),case=i)
        target=rows(data['final_states'][i])
        def fit(t):
            state=ref.sol(float(t)%period);r=rows(state)
            cost=np.sum((target[:,None,:]-r[None,:,:])**2,axis=2)
            a,b=linear_sum_assignment(cost)
            return float(np.mean(cost[a,b])),state,b
        grid=np.linspace(0.,period,1025)[:-1];values=np.array([fit(t)[0] for t in grid])
        minima=[j for j in range(len(grid)) if values[j]<=values[(j-1)%len(grid)] and values[j]<=values[(j+1)%len(grid)]]
        candidates=[]
        for j in minima:
            opt=minimize_scalar(lambda t:fit(t)[0],bounds=(grid[j]-period/1024,grid[j]+period/1024),method='bounded',options=dict(xatol=1e-13))
            candidates.append((fit(opt.x)[0],float(opt.x)))
        _,phase=min(candidates)
        for _ in range(16):
            value,state,assignment=fit(phase);delta=rows(state)[assignment]-target
            velocity=rows(motor_rows.rhs(0.,state,p,8))[assignment]
            correction=float(np.sum(delta*velocity)/np.sum(velocity*velocity))
            if fit(phase-correction)[0]>value:
                break
            phase-=correction
            if abs(correction)<1e-14:
                break
        value,state,assignment=fit(phase)
        follow=solve_ivp(lambda t,y:motor_rows.rhs(t,y,p,8),(0.,period),data['final_states'][i],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.025)
        if not follow.success:
            raise RuntimeError(follow.message)
        outcomes.append(dict(case=i,seed=case['seed'],original_lane=case['original_lane'],
            mean_ftle=float(data['lle'][i]),four2250_tau_windows=data['lle_windows'][i].tolist(),
            screen_positive=bool(data['lle'][i]>.005 and np.all(data['lle_windows'][i]>.001)),
            final_rows=row_summary(data['final_states'][i]),saddle_phase=phase%period,
            target_row_to_reference_row=assignment.tolist(),saddle_full_row_rms_distance=float(np.sqrt(value)),
            tested_reference_period=period,one_period_full_state_closure_error=float(np.max(abs(follow.y[:,-1]-data['final_states'][i])))))
        print('case',i,'saddle_distance',outcomes[-1]['saddle_full_row_rms_distance'],'closure',outcomes[-1]['one_period_full_state_closure_error'],flush=True)
    if not all(c['passed'] for c in checks):
        raise ArithmeticError([c for c in checks if not c['passed']])
    args.output.mkdir(parents=True)
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed',checks=checks,outcomes=outcomes,
        input_hashes={name:sha(args.input/name) for name in ['protocol.json','outcomes.json','outcomes.npz']},
        reference_path=args.reference.as_posix(),reference_sha256=sha(args.reference),reference_period=period,
        source_hashes={Path(__file__).name:sha(Path(__file__)),Path(motor_rows.__file__).name:sha(Path(motor_rows.__file__))},
        limitations=['Finite declared ensemble; no basin-volume inference.',
            'Near-saddle geometry and one-period closure are floating-point diagnostics, not interval inclusion.',
            'Row permutations are mathematical symmetries of this equal-row model, not claims of physical equivalence of labeled axonemes.',
            'A different numerical trajectory from an identical source start can enter a roundoff-synchronized saddle. No asymptotic basin assignment follows.'])
    (args.output/'analysis.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    for source in [Path(__file__),Path(motor_rows.__file__),ROOT/'research/report_motor_activity_gpu.py']:
        (args.output/source.name).write_bytes(source.read_bytes())
    print('passed',len(checks),'checks')


if __name__=='__main__':
    main()
