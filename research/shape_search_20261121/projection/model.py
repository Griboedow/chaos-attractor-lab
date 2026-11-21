"""Independent binary64 map and Jacobian using the two principal tangent directions."""
import math
import numpy as np
from numba import njit


@njit(cache=True)
def torus_frame(vx,vy,vz,tube,second,dp):
    if second:
        vx,vy,vz=vx-1,vz,vy
    rho=math.sqrt(vx*vx+vy*vy)
    length=math.sqrt((rho-1)**2+vz*vz)
    if min(rho,length)<1e-13:
        raise ValueError('Torus nearest point is ill-conditioned')
    cs=vx/rho
    sn=vy/rho
    nr=(rho-1)/length
    nz=vz/length
    a=(1+tube*nr)/rho
    b=tube/length
    tangent1=np.array([-sn,cs,0.])
    tangent2=np.array([nz*cs,nz*sn,-nr])
    projection=np.array([(1+tube*nr)*cs,(1+tube*nr)*sn,tube*nz])
    if second:
        projection=np.array([1+projection[0],projection[2],projection[1]])
        tangent1=np.array([-sn,0.,cs])
        tangent2=np.array([nz*cs,-nr,nz*sn])
    for i in range(3):
        for j in range(3):
            dp[i,j]=a*tangent1[i]*tangent1[j]+b*tangent2[i]*tangent2[j]
    return projection,rho,length


@njit(cache=True)
def evaluate(x, sites, beta, inv_tau, tube, out, jac, weights, cov, dp, geometry=0, mix_temp=0.):
    n=len(sites)
    maximum=-np.inf
    for k in range(n):
        s=0.
        for j in range(3):
            delta=x[j]-sites[k,j]
            s-=delta*delta
        weights[k]=s*inv_tau
        maximum=max(maximum,weights[k])
    den=0.
    m=np.zeros(3)
    for k in range(n):
        w=math.exp(weights[k]-maximum)
        weights[k]=w
        den+=w
        for j in range(3):
            m[j]+=w*sites[k,j]
    for j in range(3):
        m[j]/=den
    for i in range(3):
        for j in range(i,3):
            v=0.
            for k in range(n):
                v+=weights[k]*(sites[k,i]-m[i])*(sites[k,j]-m[j])
            cov[i,j]=cov[j,i]=2*inv_tau*v/den
    vx=2*m[0]-x[0]
    vy=2*m[1]-x[1]
    vz=2*m[2]-x[2]
    second=False
    margin=np.inf
    dist_a=dist_b=0.
    if geometry:
        # Compare distances to the surfaces, then differentiate the chosen branch.
        dist_a=abs(math.hypot(math.hypot(vx,vy)-1,vz)-tube)
        dist_b=abs(math.hypot(math.hypot(vx-1,vz)-1,vy)-tube)
        second=dist_b<dist_a
        margin=abs(dist_a-dist_b)
    if geometry==2:
        dp2=np.empty((3,3))
        pa,ra,la=torus_frame(vx,vy,vz,tube,False,dp)
        pb,rb,lb=torus_frame(vx,vy,vz,tube,True,dp2)
        weight=1/(1+math.exp(min(700.,max(-700.,(dist_a**2-dist_b**2)/mix_temp))))
        projection=weight*pa+(1-weight)*pb
        coefficient=2*weight*(1-weight)/mix_temp
        for i in range(3):
            for j in range(3):
                dp[i,j]=weight*dp[i,j]+(1-weight)*dp2[i,j]+coefficient*(pa[i]-pb[i])*(pa[j]-pb[j])
        rho=min(ra,rb)
        length=min(la,lb)
    else:
        projection,rho,length=torus_frame(vx,vy,vz,tube,second,dp)
    for i in range(3):
        out[i]=beta*(x[i]+projection[i])+(1-2*beta)*m[i]
    for i in range(3):
        for j in range(3):
            value=0.
            for k in range(3):
                value+=dp[i,k]*cov[k,j]
            jac[i,j]=beta*((i==j)-dp[i,j]+2*value)+(1-2*beta)*cov[i,j]
    return rho,length,margin,second


@njit(cache=True)
def replay(raw, initial, burn, window, windows, saved, full):
    beta,inv_tau,tube,count=raw[0]
    sites=np.ascontiguousarray(raw[2:2+int(count),:3])
    x=initial.copy()
    out=np.zeros(3)
    jac=np.empty((3,3))
    cov=np.empty((3,3))
    dp=np.empty((3,3))
    w=np.zeros(len(sites))
    q=np.eye(3)
    q[0,0]=.271
    q[1,0]=.534
    q[2,0]=-.801
    q[:,0]/=np.linalg.norm(q[:,0])
    nextq=np.zeros((3,3))
    rank=3 if full else 1
    lengths=np.empty(rank)
    rates=np.zeros((windows,rank))
    cloud=np.empty((saved,3))
    minimum=np.ones(3)*np.inf
    branch_counts=np.zeros(2,dtype=np.int64)
    maxima=0.
    for step in range(burn+window*windows):
        rho,length,margin,second=evaluate(x,sites,beta,inv_tau,tube,out,jac,w,cov,dp,int(raw[1,3]),raw[2,3])
        minimum[0]=min(minimum[0],rho)
        minimum[1]=min(minimum[1],length)
        minimum[2]=min(minimum[2],margin)
        for j in range(rank):
            for i in range(3):
                v=0.
                for h in range(3):
                    v+=jac[i,h]*q[h,j]
                nextq[i,j]=v
            # Twice-modified Gram–Schmidt for a full exponent estimate.
            for rep in range(2):
                for h in range(j):
                    u=0.
                    for i in range(3):
                        u+=nextq[i,h]*nextq[i,j]
                    for i in range(3):
                        nextq[i,j]-=u*nextq[i,h]
            sq=0.
            for i in range(3):
                sq+=nextq[i,j]**2
            lengths[j]=math.sqrt(sq)
            if lengths[j]<1e-290:
                raise ValueError('Unresolved tangent rank')
            for i in range(3):
                nextq[i,j]/=lengths[j]
        for i in range(3):
            for j in range(rank):
                q[i,j]=nextq[i,j]
            x[i]=out[i]
        if step>=burn:
            t=step-burn
            branch_counts[int(second)]+=1
            maxima=max(maxima,math.sqrt(x[0]**2+x[1]**2+x[2]**2))
            for j in range(rank):
                rates[t//window,j]+=math.log(lengths[j])/window
            if t>=window*windows-saved:
                cloud[t-(window*windows-saved)]=x
    return cloud,rates,minimum,maxima,branch_counts


def independent(x,raw):
    beta,inv_tau,tube,count=raw[0]
    b=raw[2:2+int(count),:3]
    ln=-np.sum((x-b)**2,axis=1)*inv_tau
    w=np.exp(ln-ln.max())
    m=(w@b)/w.sum()
    v=2*m-x
    center=np.array([v[0],v[1],0.])/np.linalg.norm(v[:2])
    proj=center+tube*(v-center)/np.linalg.norm(v-center)
    if raw[1,3]>0:
        c=np.array([1.,0.,0.])
        delta=v-c
        center2=c+delta*np.array([1.,0.,1.])/np.linalg.norm(delta[[0,2]])
        proj2=center2+tube*(v-center2)/np.linalg.norm(v-center2)
        if raw[1,3]==2:
            logs=-np.array([np.sum((proj-v)**2),np.sum((proj2-v)**2)])/raw[2,3]
            weights=np.exp(logs-logs.max())
            proj=(weights[0]*proj+weights[1]*proj2)/weights.sum()
        elif np.linalg.norm(proj2-v)<np.linalg.norm(proj-v):
            proj=proj2
    # Nested reflector and averaging form.
    return beta*(2*proj-v+x)/2+(1-beta)*m


def derivative_audit(raw,states):
    beta,inv_tau,tube,count=raw[0]
    n=int(count)
    b=raw[2:2+n,:3].copy()
    err=0.
    fd=0.
    for x in states:
        a=np.zeros(3)
        jac=np.zeros((3,3))
        evaluate(x,b,beta,inv_tau,tube,a,jac,np.zeros(n),np.zeros((3,3)),np.zeros((3,3)),int(raw[1,3]),raw[2,3])
        err=max(err,float(np.max(np.abs(a-independent(x,raw)))))
        hs=2e-6*max(1.,np.linalg.norm(x))
        diff=np.column_stack([(independent(x+v*hs,raw)-independent(x-v*hs,raw))/(2*hs) for v in np.eye(3)])
        fd=max(fd,float(np.linalg.norm(jac-diff)/max(1.,np.linalg.norm(jac))))
    return err,fd
