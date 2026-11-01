"""Finite visits between the thirty outer saddle foci of case 177."""
from pathlib import Path
import sys
from structure_paths import ROOT,RESEARCH,INPUT,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json
import numpy as np
from scipy.spatial import cKDTree

base=OUTPUT;inv=json.loads((base/'icosa-equilibria/inventory.json').read_text());case=next(c for c in inv if c['id']==177);out=base/'icosa-itinerary-final';out.mkdir(exist_ok=False)
focus=np.array([r['point'] for r in case['roots'] if abs(r['radius']-1.0408942)<1e-6]);assert len(focus)==30
z=np.load(INPUT/'icosa-long/177-dense.npy',mmap_mode='r');dist,idx=cKDTree(focus).query(z);radius=.12
events=[];last=-1
for step,label in zip(np.flatnonzero(dist<radius),idx[dist<radius]):
    if label!=last:events.append((int(step),int(label)));last=label
events=np.array(events);counts=np.zeros((30,30),dtype=int)
for a,b in zip(events[:-1,1],events[1:,1]):counts[a,b]+=1
directions=focus/np.linalg.norm(focus,axis=1)[:,None];angles=np.arccos(np.clip(directions@directions.T,-1,1))
edges=[{'from':int(i),'to':int(j),'count':int(counts[i,j]),'angle_degrees':float(np.degrees(angles[i,j]))} for i,j in zip(*np.where(counts>0))]
hist=np.bincount(events[:,1],minlength=30)
np.savez_compressed(out/'data.npz',foci=focus,events=events,transition_counts=counts,nearest_distance=dist[::16])
report={'id':177,'radius_threshold':radius,'time_per_saved_state':.0125,'events':len(events),'visited_foci':int((hist>0).sum()),'visits_per_focus':hist.tolist(),'directed_edges':len(edges),'edges':edges,'fraction_near_focus':float((dist<radius).mean()),'limits':'Finite neighborhood-entry itinerary. Consecutive visits to the same focus are merged, even after excursions outside all neighborhoods. It is not a Markov partition, a certified connection graph or a topological invariant.'}
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='edges'}),flush=True)
