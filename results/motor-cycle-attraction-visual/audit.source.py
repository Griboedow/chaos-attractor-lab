import datetime,hashlib,json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
base=root/'results/motor-cycle-attraction-visual'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
checks={};summary=[]
for name,amplitude in [('.',.6),('local-0p001',.001),('local-0p03',.03)]:
    p=base/name;protocol=json.loads((p/'protocol.json').read_text());analysis=json.loads((p/'analysis.json').read_text())
    checks[name+':amplitude']=protocol['full24D_perturbation_norm']==amplitude
    checks[name+':source_hash']=sha(p/'source.py')==protocol['source_sha256']
    checks[name+':field_hash']=sha(p/'motor_rows.py')==protocol['field_source_sha256']
    checks[name+':certificate_hash']=sha(root/protocol['certificate'])==protocol['certificate_sha256']
    checks[name+':protocol_hash']=sha(p/'protocol.json')==analysis['protocol_sha256']
    checks[name+':reference_hash']=sha(p/'reference.npz')==analysis['reference_sha256']
    checks[name+':all_starts']=len(analysis['runs'])==len(protocol['initial_states'])==6
    for i,r in enumerate(analysis['runs']):
        file=p/f'trajectory-{i:02}.npz';data=np.load(file)
        checks[f'{name}:{i}:hash']=sha(file)==r['trajectory_sha256']
        checks[f'{name}:{i}:initial']=np.array_equal(data['initial_state'],np.array(protocol['initial_states'][i])) and np.array_equal(data['main_states'][0],data['initial_state']) and np.array_equal(data['strict_states'][0],data['initial_state'])
        checks[f'{name}:{i}:complete']=data['strict_states'].shape==(3001,23) and abs(data['time'][-1]-6*protocol['period'])<1e-12
        error=float(np.max(abs(data['main_states']-data['strict_states'])))
        checks[f'{name}:{i}:control_record']=error==r['solver_max_state_difference'] and (error<1e-7)==r['integration_control_passed']
        checks[f'{name}:{i}:distance_record']=data['full_state_orbit_distance'][-1]==r['final_distance_to_orbit']
    summary.append(dict(folder=name,amplitude=amplitude,controlled=sum(v['integration_control_passed'] for v in analysis['runs']),
                        approaching_reference=sum(v['approached_reference_cycle'] for v in analysis['runs'])))
checks['displayed_six_pass']=summary[-1]['controlled']==summary[-1]['approaching_reference']==6
checks['small_six_pass']=summary[1]['controlled']==summary[1]['approaching_reference']==6
render=root/'docs/assets/motor-cycle-attraction.json';r=json.loads(render.read_text())
checks['renderer_source_hash']=sha(root/'research/render_motor_cycle_attraction.py')==r['source_sha256']
for path,h in r['inputs'].items():checks['render_input:'+path]=sha(root/path)==h
for path,v in r['outputs'].items():checks['render_output:'+path]=sha(render.parent/path)==v['sha256']
checks={k:bool(v) for k,v in checks.items()}
record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),passed=all(checks.values()),
            checks=checks,attempts=summary,scope='Data/provenance consistency and complete experiment accounting; not a validated basin proof.')
(base/'audit.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
(base/'audit.source.py').write_bytes(Path(__file__).read_bytes())
files=[p for p in base.rglob('*') if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts]
files += [root/'research/illustrate_motor_cycle_attraction.py',root/'research/render_motor_cycle_attraction.py',render]
files += [render.parent/p for p in r['outputs']]
manifest=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed' if all(checks.values()) else 'failed',
              checks=len(checks),file_count=len(files),total_bytes=sum(p.stat().st_size for p in files),
              files=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)])
(base/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(passed=record['passed'],checks=len(checks),attempts=summary,files=len(files),bytes=manifest['total_bytes'])))
assert record['passed']
