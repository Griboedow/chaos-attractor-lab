"""Ross et al. 2025, full five-state continuous chemostat author model.

Pinned author model_full.m: bde791b4c82cd3bebe8d2b34553d8e50c7b47d01.
States N1,N2,R1,R2,R3; parameter order follows code assignments, not its stale
header. min growth limitation is retained literally. Time tau=D*t.
"""
import numpy as np
from numba import njit

SOURCE='https://doi.org/10.1038/s41467-025-63986-y'
AUTHOR_COMMIT='bde791b4c82cd3bebe8d2b34553d8e50c7b47d01'
REFERENCE=np.array([1.,10.,10.,1.,10.,20.,1.,10.,10.,1.,10.,20.,.2,1.,1.,11.10])
REFERENCE_START=np.array([.1,.11,1.,1.,11.10])
PARAMETERS=['mu1_max','k13','k12','q13','q12','q11','mu2_max','k23','k21','q23','q21','q22','D','R1_in','R2_in','R3_in']


@njit(cache=True)
def field_jac(x,p):
    a,b,u,v,g=x
    m12=p[0]*v/(p[2]+v);m13=p[0]*g/(p[1]+g)
    m21=p[6]*u/(p[8]+u);m23=p[6]*g/(p[7]+g)
    d12=p[0]*p[2]/(p[2]+v)**2;d13=p[0]*p[1]/(p[1]+g)**2
    d21=p[6]*p[8]/(p[8]+u)**2;d23=p[6]*p[7]/(p[7]+g)**2
    m1=min(m12,m13);m2=min(m21,m23);D=p[12]
    g1=np.zeros(5);g2=np.zeros(5);gg1=np.zeros(5);gg2=np.zeros(5)
    if m12<m13:g1[3]=d12
    else:g1[4]=d13
    if m21<m23:g2[2]=d21
    else:g2[4]=d23
    gg1[4]=d13;gg2[4]=d23
    f=np.array([(m1-D)*a,(m2-D)*b,
        D*(p[13]-u)+p[5]*(m13-m1)*a-p[10]*m2*b,
        D*(p[14]-v)+p[11]*(m23-m2)*b-p[4]*m1*a,
        D*(p[15]-g)-p[3]*m13*a-p[9]*m23*b])
    j=np.empty((5,5))
    j[0]=a*g1;j[0,0]+=m1-D
    j[1]=b*g2;j[1,1]+=m2-D
    j[2]=p[5]*a*(gg1-g1)-p[10]*b*g2
    j[2,0]+=p[5]*(m13-m1);j[2,1]-=p[10]*m2;j[2,2]-=D
    j[3]=p[11]*b*(gg2-g2)-p[4]*a*g1
    j[3,0]-=p[4]*m1;j[3,1]+=p[11]*(m23-m2);j[3,3]-=D
    j[4]=-p[3]*a*gg1-p[9]*b*gg2
    j[4,0]-=p[3]*m13;j[4,1]-=p[9]*m23;j[4,4]-=D
    return f/D,j/D


def independent_field(x,p):
    """Four physical fluxes and a constant stoichiometric matrix."""
    n1,n2,r1,r2,g=x
    amino1=p[0]*r2/(p[2]+r2);glucose1=p[0]*g/(p[1]+g)
    amino2=p[6]*r1/(p[8]+r1);glucose2=p[6]*g/(p[7]+g)
    # Real-branch selection permits complex-step checks away from switching.
    growth1=amino1 if amino1.real<glucose1.real else glucose1
    growth2=amino2 if amino2.real<glucose2.real else glucose2
    flux=np.array([growth1*n1,growth2*n2,glucose1*n1,glucose2*n2])
    stoich=np.array([[1,0,0,0],[0,1,0,0],[-p[5],-p[10],p[5],0],
                     [-p[4],-p[11],0,p[11]],[0,0,-p[3],-p[9]]])
    return stoich@flux/p[12]+np.r_[0.,0.,p[13:16]]-x


def balance_weights(p):
    """Positive exact linear budget: B'=B_in-B in tau units."""
    return np.array([p[3]+p[4]*p[9]/p[11],p[9]+p[10]*p[3]/p[5],
                     p[3]/p[5],p[9]/p[11],1.])


def balance_in(p):
    return balance_weights(p)[2:]@p[13:16]


def branch_gaps(x,p):
    return np.array([p[0]*(x[3]/(p[2]+x[3])-x[4]/(p[1]+x[4])),
                     p[6]*(x[2]/(p[8]+x[2])-x[4]/(p[7]+x[4]))])
