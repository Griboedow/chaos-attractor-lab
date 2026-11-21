"""Robust sampled covering relations between monotone laps of the fifth-return chart."""
import json
import numpy as np
from scipy.linalg import eigvals
from scipy.sparse.csgraph import connected_components
from native_section import RAW,RESULT
from return_tube import graph,CENTER

def main():
    file=RAW/'tube-n4097-l17-h0.00625.npz'
    a=np.load(file);z=a['z'][:,0];u=a['returns'][:,:,2]
    g=graph();d=a['returns'][:,:,1]-g(u-CENTER)
    slope=np.diff(u[:,u.shape[1]//2])/np.diff(z)
    turns=np.where(slope[1:]*slope[:-1]<0)[0]+1
    endpoints=np.r_[0,turns,len(z)-1]
    # Gaps between lap intervals enforce separated symbolic itineraries.
    left=endpoints[:-1]+2;right=endpoints[1:]-2
    ok=right>left+1;left=left[ok];right=right[ok];n=len(left)
    M=np.zeros((n,n),int);margin=4e-5
    witness=[]
    for i,(l,r) in enumerate(zip(left,right)):
        for j,(sl,sr) in enumerate(zip(left,right)):
            mpos=min(z[sl]-u[l].max(),u[r].min()-z[sr])
            mneg=min(z[sl]-u[r].max(),u[l].min()-z[sr])
            m=max(mpos,mneg)
            if m>margin:
                M[i,j]=1
                witness.append([int(i),int(j),float(m),1 if mpos>mneg else -1])
    ncomp,labels=connected_components(M,connection='strong')
    comps=[]
    for k in range(ncomp):
        which=np.flatnonzero(labels==k)
        A=M[np.ix_(which,which)];rho=max(abs(eigvals(A)))
        if rho>1+1e-8:
            comps.append(dict(indices=which.tolist(),adjacency=A.tolist(),rho=float(rho),
              conditional_flow_entropy=float(np.log(rho)/a['returns'][:,:,0].max())))
    cert=dict(source=file.name,rectangle_gap_native=float(min(z[left[1:]]-z[right[:-1]])),edge_margin_threshold=margin,
           rectangles=np.c_[z[left],z[right]].tolist(),face_grid_indices=np.c_[left,right].tolist(),edges=witness,components=comps,
           largest_transverse_image=float(np.max(abs(d))),edge_count=len(witness),
           interpretation='Conditional on the complete vertical faces satisfying the sampled inequalities, transverse image remaining interior, and a continuous return over each rectangle. No face or flow enclosure inferred from a finite mesh.')
    (RESULT/'coverings.json').write_text(json.dumps(cert,indent=2)+'\n')
    print('laps',n,'edges',len(witness),'components',[(len(c['indices']),c['rho'],c['conditional_flow_entropy']) for c in comps])

if __name__=='__main__':main()
