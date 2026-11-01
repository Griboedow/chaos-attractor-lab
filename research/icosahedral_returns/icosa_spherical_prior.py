"""Compare finite edge orientation with an explicitly credited integrable field."""
from pathlib import Path
from return_paths import INPUT,STRUCTURE,CYCLES,OUTPUT
import json
import numpy as np

base=OUTPUT;out=base/'icosa-spherical-prior';out.mkdir(exist_ok=False)
d=np.load(INPUT/'icosa-v1/input.npz');group=d['group'];foci=np.load(STRUCTURE/'icosa-itinerary-final/data.npz')['foci'];v=foci/np.linalg.norm(foci,axis=1)[:,None];A=np.load(CYCLES/'icosa-directed-cycles/data.npz')['counts']>0
p=(1+np.sqrt(5))/2;S=np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,1.]])
def invariant(x):
    a,b,c=x@S.T;return (p*p*a*a-b*b)*(p*p*b*b-c*c)*(p*p*c*c-a*a)
def reference(x):
    # Pull the author's polar vector field back through the explicit x/y swap.
    a,b,c=S@x;P=p*p*a*a-b*b;Q=p*p*b*b-c*c;R=p*p*c*c-a*a
    grad=np.array([2*a*(p*p*Q*R-P*Q),2*b*(-Q*R+p*p*P*R),2*c*(-P*R+p*p*P*Q)])
    return S.T@((2.5-np.sqrt(5))*np.cross(2*np.array([a,b,c]),grad))
x=np.array([.31,-.82,.4]);invariance=float(max(abs(invariant(G@x)-invariant(x)) for G in group));assert invariance<1e-14
rows=[]
for i,j in np.argwhere(A):
    mid=v[i]+v[j];mid/=np.linalg.norm(mid);tangent=v[j]-v[i];tangent/=np.linalg.norm(tangent);F=reference(mid);F/=np.linalg.norm(F)
    rows.append({'from':int(i),'to':int(j),'reference_direction_dot_observed_edge':float(F@tangent),'midpoint_invariant':float(invariant(mid))})
z=np.load(INPUT/'icosa-long/sample.npz')['tail'][list(np.load(INPUT/'icosa-long/sample.npz')['ids']).index(177)][::8];r=np.linalg.norm(z,axis=1);levels=np.array([invariant(a/np.linalg.norm(a)) for a in z])
report={'source':'https://doi.org/10.13140/RG.2.2.20034.95689','source_title':'Giedrius Alkauskas, The simplest icosahedral flow on the sphere (26 June 2026 preprint)','read_scope':'Author-uploaded full-text extraction, Sections 1 and 2; PDF link retrieval returned 404. No figure-level comparison is claimed.','coordinate_pullback':S.tolist(),'invariant_group_error':invariance,'edge_count':len(rows),'minimum_direction_dot':min(a['reference_direction_dot_observed_edge'] for a in rows),'maximum_direction_dot':max(a['reference_direction_dot_observed_edge'] for a in rows),'max_abs_midpoint_level':max(abs(a['midpoint_invariant']) for a in rows),'our_sampled_radius_range':[float(r.min()),float(r.max())],'our_normalized_sextic_range':[float(levels.min()),float(levels.max())],'edges':rows,'limits':'The integrable reference preserves radius and sextic level. A directed graph agreement would be prior geometric structure, not an equivalence of full flows or a historical novelty certificate.'}
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='edges'}),flush=True)
