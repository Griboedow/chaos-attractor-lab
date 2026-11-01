"""Fourth iterate of the sphere/posterior map for periodic-saddle proposals."""
import numpy as np
from numba import njit
from search_soft_sphere import step


@njit(cache=True)
def iterate(x,beta,tau,points,count):
    J=np.eye(3);minimum=np.inf
    for k in range(4*count):
        x,A,r,_,_,flag=step(x,beta,points,tau);minimum=min(minimum,r)
        if flag or not np.isfinite(x).all():return x,J,minimum,1
        J=A@J
    return x,J,minimum,0
