"""Add explicit, cumulative continuation of saved BDF tangent histories."""
from pathlib import Path
root=Path(__file__).resolve().parent
source=(root/'tangent.py').read_text(encoding='utf-8')
changes=[
("dt=float(c['dt']);origin=float(c['time'])", "dt=float(c['dt']);origin=float(c['time'])\n        saved_tangent={name:c[name].copy() for name in ['dk','dk_previous','dnodes','dnodes_previous','log_growth']} if args.resume_tangent else None"),
("dkprev=dk.copy();dn=np.zeros((n+1,3,count));dnprev=dn.copy()", "dkprev=dk.copy();dn=np.zeros((n+1,3,count));dnprev=dn.copy()\n    initial_observed=0.0\n    if args.resume_tangent:\n        assert args.burn==0\n        dk=saved_tangent['dk'];dkprev=saved_tangent['dk_previous']\n        dn=saved_tangent['dnodes'];dnprev=saved_tangent['dnodes_previous']\n        initial_observed=seed_report['observed_duration']"),
("totals=np.zeros(count);observations=[];success=True;failure=None", "totals=saved_tangent['log_growth'].copy() if args.resume_tangent else np.zeros(count)\n    observations=[];success=True;failure=None"),
("observed_duration=max(0,(step-burn)*dt)", "observed_duration=initial_observed+max(0,(step-burn)*dt)"),
("totals/((step-burn)*dt)", "totals/(initial_observed+(step-burn)*dt)"),
("observed_duration=max(0,(completed-burn)*dt)", "initial_observed_duration=initial_observed,new_observed_duration=max(0,(completed-burn)*dt),\n        observed_duration=initial_observed+max(0,(completed-burn)*dt)"),
("totals/((completed-burn)*dt)", "totals/(initial_observed+(completed-burn)*dt)"),
("initial_tangents='Three independent Gaussian curvature directions, copied to current and previous history; zero initial nodal variations; algebraic startup excluded by the burn interval.'", "initial_tangents='Restored current and previous tangent histories and cumulative log growth' if args.resume_tangent else 'Three independent Gaussian curvature directions; algebraic startup excluded by burn'"),
("ap.add_argument('--random-seed'", "ap.add_argument('--resume-tangent',action='store_true')\n    ap.add_argument('--random-seed'")]
for before,after in changes:
    assert source.count(before)==1,(before,source.count(before))
    source=source.replace(before,after)
(root/'tangent_extended.py').write_text(source,encoding='utf-8')
