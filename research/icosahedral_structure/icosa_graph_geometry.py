"""Identify the finite visit graph against the geometry of its saddle directions."""
from pathlib import Path
import sys
from structure_paths import ROOT,RESEARCH,INPUT,OUTPUT
sys.path.insert(0,str(RESEARCH/'unexpected_search'))
import json
import numpy as np
from scipy.spatial import ConvexHull

base=OUTPUT;d=np.load(base/'icosa-itinerary-final/data.npz');out=base/'icosa-graph-geometry';out.mkdir(exist_ok=False)
points=d['foci'];directions=points/np.linalg.norm(points,axis=1)[:,None];dots=directions@directions.T
edge=(np.abs(dots-np.cos(np.pi/5))<1e-8);np.fill_diagonal(edge,False);assert edge.sum()==120 and np.all(edge.sum(1)==4)
hull=ConvexHull(directions);faces=set()
for eq in hull.equations:
    face=tuple(np.flatnonzero(np.abs(directions@eq[:3]+eq[3])<1e-8).tolist());faces.add(face)
sizes=[len(f) for f in faces];assert sizes.count(3)==20 and sizes.count(5)==12 and len(faces)==32
counts=d['transition_counts'];assert np.all(counts[~edge]==0)
missing=np.argwhere(edge&(counts==0));angles=np.degrees(np.arccos(np.clip(dots[counts>0],-1,1)))
report={'vertices':30,'undirected_edges':60,'vertex_degree':4,'triangular_faces':20,'pentagonal_faces':12,'geometry':'Icosahedron edge-midpoint directions: an icosidodecahedron. This is a numerical coordinate/convex-hull identification, not a formal polyhedral proof.','observed_transitions':int(counts.sum()),'observed_directed_edges':int((counts>0).sum()),'possible_directed_edges':120,'missing_directed_edges':missing.tolist(),'min_transition_angle_degrees':float(angles.min()),'max_transition_angle_degrees':float(angles.max()),'non_edge_transitions':int(counts[~edge].sum()),'limits':'A finite graph describing this trajectory and this neighborhood threshold. It is not an invariant transition rule for all orbits or a chaos proof.'}
(out/'report.json').write_text(json.dumps(report,indent=2));np.savez_compressed(out/'geometry.npz',directions=directions,edge_mask=edge,triangles=np.array([f for f in sorted(faces) if len(f)==3]),pentagons=np.array([f for f in sorted(faces) if len(f)==5]))
print(json.dumps(report),flush=True)
