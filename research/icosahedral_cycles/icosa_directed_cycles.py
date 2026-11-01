"""Check threshold sensitivity, oriented face cycles, and numerical orbit itineraries."""
from pathlib import Path
import sys
from cycle_paths import ROOT,RESEARCH,INPUT,STRUCTURE,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json
import numpy as np
from scipy.spatial import cKDTree

base=OUTPUT;out=base/'icosa-directed-cycles';out.mkdir(exist_ok=False)
foci=np.load(STRUCTURE/'icosa-itinerary-final/data.npz')['foci'];tree=cKDTree(foci)
geo=np.load(STRUCTURE/'icosa-graph-geometry/geometry.npz');edge=geo['edge_mask']
z=np.load(INPUT/'icosa-long/177-dense.npy',mmap_mode='r');dist,idx=tree.query(z)
def events(dist,idx,r):
    sample=np.flatnonzero(dist<r);ids=idx[sample];keep=np.r_[True,ids[1:]!=ids[:-1]]
    return np.column_stack((sample[keep],ids[keep]))
def counts(ev):
    c=np.zeros((30,30),int);np.add.at(c,(ev[:-1,1],ev[1:,1]),1);return c
evs={str(r):events(dist,idx,r) for r in [.12,.16,.2,.25,.3]};C=counts(evs['0.2']);A=C>0
faces=[]
for label in ['triangles','pentagons']:
    for face in geo[label]:
        sub=A[np.ix_(face,face)];directed=bool(np.all(sub.sum(0)==1) and np.all(sub.sum(1)==1))
        faces.append({'vertices':face.tolist(),'size':len(face),'directed_cycle':directed})
thresholds=[]
for key,ev in evs.items():
    c=counts(ev);thresholds.append({'radius':float(key),'events':len(ev),'directed_edges':int((c>0).sum()),'same_sequence_as_02':bool(np.array_equal(ev[:,1],evs['0.2'][:,1])),'reverse_edges_vs_02':int(((c>0)&A.T).sum()),'reverse_transitions_vs_02':int(c[A.T].sum())})
# Follow each radius-.12 change through the radius-.2 visits at its exact sample times.
narrow=evs['0.12'];wide=evs['0.2'];classification={};examples=[]
for (t0,a),(t1,b) in zip(narrow[:-1],narrow[1:]):
    if A[a,b]:continue
    left=np.searchsorted(wide[:,0],t0,side='right')-1;right=np.searchsorted(wide[:,0],t1,side='right')-1
    path=wide[left:right+1,1];key=str(len(path)-1);classification[key]=classification.get(key,0)+1
    assert path[0]==a and path[-1]==b and np.all(A[path[:-1],path[1:]])
    if len(examples)<5:examples.append({'narrow':[int(a),int(b)],'wide_path':path.tolist()})
cycles=[]
for index in [0,20]:
    p=np.load(base/f'icosa-cycle-refinement/cycle-{index}.npz')['points'];cd,ci=tree.query(p);ev=events(cd,ci,.2);seq=ev[:,1]
    if seq[0]==seq[-1]:seq=seq[:-1]
    cycles.append({'index':index,'visits_cyclic':seq.tolist(),'distinct_centers':len(np.unique(seq)),'matches_face':any(set(seq.tolist())==set(f['vertices']) and len(seq)==f['size'] for f in faces),'all_edges_follow_02_orientation':bool(np.all(A[seq,np.roll(seq,-1)])),'min_distance_to_focus':float(cd.min()),'max_distance_to_nearest_focus':float(cd.max())})
report={'thresholds':thresholds,'one_direction_per_geometric_edge':bool(np.array_equal(A.astype(int)+A.T.astype(int),edge.astype(int))),'outdegree':A.sum(1).tolist(),'indegree':A.sum(0).tolist(),'faces':faces,'missed_visit_path_lengths':classification,'missed_visit_examples':examples,'cycles':cycles,'limits':'Finite orbit and numerical cycles only; no invariant partition, heteroclinic connection, or certified periodic orbit is asserted.'}
np.savez_compressed(out/'data.npz',counts=C,**{f'events_{k}':v for k,v in evs.items()})
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
