"""Outward interval extension for the finite Gaussian mean / circle map.

The exponential uses a degree-12 Taylor polynomial on [-1/8,1/8],
outward arithmetic, an absolute 2^-71 remainder, and repeated squaring.
Indeed exp(1/8)*(1/8)^13/13! <= (8/7)*(1/8)^13/13! < 2^-71.
No libm exponential or assumption about its ulp error enters this extension.
Basic binary64 arithmetic and sqrt use the same nextafter expansion as the
existing circle certificates; no fastmath or parallel reductions are used.
Parameters and sites represent exact binary64 real numbers.
"""
import numpy as np
from numba import njit
from raar_interval_scalar import add,sub,mul,square,reciprocal,root,up,down


@njit(cache=False,inline='always')
def point(x):return (x,x)


@njit(cache=False)
def exp_point(x):
    assert np.isfinite(x) and abs(x)<=700.
    scaled=x;count=0
    while abs(scaled)>.125:scaled*=.5;count+=1
    result=point(1.)
    for k in range(12,0,-1):
        result=add(point(1.),mul(mul(point(scaled),reciprocal(point(float(k)))),result))
    result=add(result,(-2.**-71,2.**-71))
    for k in range(count):
        result=square(result);result=(max(0.,result[0]),result[1])
    return result


@njit(cache=False)
def exponential(x):return (exp_point(x[0])[0],exp_point(x[1])[1])


@njit(cache=False)
def evaluate(box,beta,tau,points):
    n=len(points);mid=(box[:,0]+box[:,1])*.5
    anchor=np.argmax(2*points@mid-np.sum(points**2,axis=1))
    weights=np.empty((n,2));total=point(0.)
    for k in range(n):
        logit=point(0.)
        if k!=anchor:
            for j in range(3):
                delta=sub(point(points[k,j]),point(points[anchor,j]))
                linear=mul(mul(point(2.),delta),box[j])
                # Squared norms are enclosed from the exact stored coordinates.
                constant=sub(square(point(points[anchor,j])),square(point(points[k,j])))
                logit=add(logit,add(linear,constant))
            logit=mul(logit,reciprocal(point(tau)))
            weights[k]=exponential(logit)
        else:weights[k]=point(1.)
        total=add(total,weights[k])
    inverse=reciprocal(total)
    mean=np.zeros((3,2))
    for k in range(n):
        weights[k]=mul(weights[k],inverse)
        for j in range(3):mean[j]=add(mean[j],mul(weights[k],point(points[k,j])))
    D=np.zeros((3,3,2));factor=mul(point(2.),reciprocal(point(tau)))
    for k in range(n):
        delta=np.empty((3,2))
        for i in range(3):delta[i]=sub(point(points[k,i]),mean[i])
        for i in range(3):
            for j in range(3):
                product=square(delta[i]) if i==j else mul(delta[i],delta[j])
                D[i,j]=add(D[i,j],mul(weights[k],product))
    for i in range(3):
        for j in range(3):D[i,j]=mul(factor,D[i,j])
    v=np.empty((3,2))
    for i in range(3):v[i]=sub(mul(point(2.),mean[i]),box[i])
    radius2=add(square(v[0]),square(v[1]));image=np.zeros((3,2));J=np.zeros((3,3,2))
    if radius2[0]<=0:return image,J,0.,1
    radius=root(radius2);inverse=reciprocal(radius)
    u=np.zeros((3,2));u[0]=mul(v[0],inverse);u[1]=mul(v[1],inverse)
    tangent=np.zeros((3,2));tangent[0]=(-u[1,1],-u[1,0]);tangent[1]=u[0]
    coefficient=sub(point(1.),mul(point(2.),point(beta)))
    for i in range(3):
        image[i]=add(add(mul(point(beta),box[i]),mul(point(beta),u[i])),mul(coefficient,mean[i]))
        for j in range(3):
            reflected=point(0.)
            for k in range(3):
                value=sub(mul(point(2.),D[k,j]),point(1. if k==j else 0.))
                reflected=add(reflected,mul(tangent[k],value))
            normalized=mul(mul(mul(point(beta),inverse),tangent[i]),reflected)
            J[i,j]=add(add(point(beta if i==j else 0.),mul(coefficient,D[i,j])),normalized)
    return image,J,radius[0],0


@njit(cache=False)
def remainder(box,mid,image_center,linear,beta,tau,points,margin):
    radii=np.empty(3);errors=np.zeros(3)
    for j in range(3):radii[j]=up(max(up(mid[j]-box[j,0]),up(box[j,1]-mid[j])))
    _,bounds,axis,flag=evaluate(box,beta,tau,points)
    if flag:return radii,errors,axis,flag
    center_box=np.empty((3,2))
    for j in range(3):center_box[j]=point(mid[j])
    center,_,_,flag=evaluate(center_box,beta,tau,points)
    if flag:return radii,errors,axis,flag
    for i in range(3):
        error=sub(center[i],point(image_center[i]));errors[i]=up(up(max(abs(error[0]),abs(error[1])))+margin)
        for j in range(3):
            diff=sub(bounds[i,j],point(linear[i,j]))
            errors[i]=up(errors[i]+up(max(abs(diff[0]),abs(diff[1]))*radii[j]))
    return radii,errors,axis,0
