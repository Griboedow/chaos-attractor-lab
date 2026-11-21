"""Numerical cubic covering charts and an independent geometric braid check."""
import itertools
import numpy as np

PERMS=np.asarray(list(itertools.permutations(range(3))))
OMEGA=np.exp(2j*np.pi/3)


def forcing(t,row):
    t=np.asarray(t)
    c=row["coeffs"];h=row["h"]
    real=np.full(t.shape,float(c[0]))
    for j in range(1,h+1):
        if c[j]:real+=c[j]*np.cos(j*t)
        if c[h+j]:real+=c[h+j]*np.sin(j*t)
    return real/row["divisor"]+1j*row["imag_amplitude"]*np.sin(row["k"]*t)


def all_roots(constant):
    a=-np.asarray(constant,dtype=complex)/2
    radical=np.sqrt(a*a-1)
    v=np.where(abs(a+radical)>abs(a-radical),a+radical,a-radical)**(1/3)
    roots=np.stack([v*q+1/(v*q) for q in (1.,OMEGA,OMEGA**2)],axis=-1)
    assert np.max(abs(roots**3-3*roots+np.asarray(constant)[...,None]))<2e-10
    return roots


def charts(row,n=4096):
    t=np.arange(n+1)*2*np.pi/n
    raw=all_roots(forcing(t,row))
    paths=np.zeros_like(raw)
    paths[0]=sorted(raw[0],key=lambda z:(z.real,z.imag))
    for j in range(1,len(raw)):
        d=abs(raw[j,:,None]-paths[j-1,None,:])
        idx=d.argmin(0)
        if len(set(idx))!=3:raise ValueError(("Grid too coarse",row["case"]))
        paths[j]=raw[j,idx]
    assert np.max(abs(paths[-1]-paths[0]))<1e-10
    return paths


def apply_chart(theta,u,row,strand,core):
    # Continue from the unperturbed sheet along a radial segment in the residual disk.
    n=len(core)-1
    phase=(theta%(2*np.pi))*n/(2*np.pi)
    index=phase.astype(int)%n;weight=phase-np.floor(phase)
    z=(1-weight)*core[index,strand]+weight*core[index+1,strand]
    b=forcing(theta,row)
    for f in np.linspace(0,1,9):
        for _ in range(5):z-=(z*z*z-3*z+b-f*u)/(3*z*z-3)
    assert np.max(abs(z**3-3*z+b-u))<2e-11
    return z


def xyz(theta,z,row):
    rad=row["embedding_radius"]+z.real
    assert np.min(rad)>0
    return np.column_stack((rad*np.cos(theta),rad*np.sin(theta),z.imag))


def geometric_word(loops,angle):
    # Project along a generic transverse direction; exactly one adjacent pair can
    # change order in a sampled segment. Interpolate its crossing depth.
    z=loops*np.exp(1j*angle)
    words=[]
    initial=np.argsort(z[0].real)
    ordering=initial.copy()
    for i in range(1,len(z)):
        after=np.argsort(z[i].real)
        if np.array_equal(after,ordering):continue
        ix=np.flatnonzero(after!=ordering)
        if len(ix)!=2 or ix[1]!=ix[0]+1:raise ValueError("Refine crossing sampling or change projection")
        j=ix[0];l,r=ordering[j:j+2]
        before=z[i-1,l]-z[i-1,r];now=z[i,l]-z[i,r]
        t=before.real/(before.real-now.real)
        depth=(1-t)*before.imag+t*now.imag
        assert abs(depth)>1.e-6 and 0<=t<=1
        words.append(int((j+1)*(1 if depth>0 else -1)))
        ordering=after
    if not np.array_equal(ordering,initial):raise ValueError("Curves fail to close in braid frame")
    return words


def invariants(word):
    # Independently propagate integer unipotent matrices and geometric crossing signs.
    I=np.eye(3,dtype=np.int64)
    meridians=[I.copy() for _ in range(3)]
    meridians[0][0,1]=1;meridians[1][1,2]=1
    longitudes=[I.copy() for _ in range(3)]
    labels=list(range(3));counts={(0,1):0,(0,2):0,(1,2):0}
    def inv(A):
        N=A-I
        return I-N+N@N
    for letter in word:
        i=abs(letter)-1;j=i+1
        counts[tuple(sorted(labels[i:j+1]))]+=int(np.sign(letter))
        X,Y=meridians[i],meridians[j]
        P,Q=longitudes[i],longitudes[j]
        if letter>0:
            meridians[i:j+1]=[X@Y@inv(X),X]
            longitudes[i:j+1]=[X@Q,P]
        else:
            meridians[i:j+1]=[Y,inv(Y)@X@Y]
            longitudes[i:j+1]=[Q,inv(Y)@P]
        labels[i:j+1]=labels[i:j+1][::-1]
    assert labels==[0,1,2] and all(v%2==0 for v in counts.values())
    mu=None if any(counts.values()) else int(inv(longitudes[2])[0,2])
    return [c//2 for c in counts.values()],mu
