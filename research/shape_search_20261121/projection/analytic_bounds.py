"""Independent interval enclosure from Gaussian covariance and torus tangent frames."""
from certify_cycles import exactfloat,eye,matmul,relaxation


def step(x,rec,ctx):
    raw=rec['parameters']
    inv=exactfloat(raw['inv_tau'],ctx);beta=relaxation(raw,ctx)
    tube=exactfloat(raw['tube_radius'],ctx)
    points=[[exactfloat(v,ctx) for v in row] for row in raw['sites']]
    logits=[(2*sum(v*b for v,b in zip(x,p))-sum(b*b for b in p))*inv for p in points]
    weights=[ctx.exp(v) for v in logits]
    total=sum(weights);w=[v/total for v in weights]
    mean=[sum(t*p[i] for t,p in zip(w,points)) for i in range(3)]
    cov=[[2*inv*sum(t*(p[i]-mean[i])*(p[j]-mean[j]) for t,p in zip(w,points)) for j in range(3)] for i in range(3)]
    v=[2*a-b for a,b in zip(mean,x)]
    def frame(second):
        vx,vy,vz=(v[0]-1,v[2],v[1]) if second else v
        axis=ctx.sqrt(vx**2+vy**2)
        if not axis>0: raise ArithmeticError('Unresolved projection axis')
        length=ctx.sqrt((axis-1)**2+vz**2)
        if not length>0: raise ArithmeticError('Unresolved centerline')
        cs=vx/axis;sn=vy/axis;nr=(axis-1)/length;nz=vz/length
        gain1=(1+tube*nr)/axis;gain2=tube/length
        a=[-sn,cs,0];b=[nz*cs,nz*sn,-nr]
        p=[(1+tube*nr)*cs,(1+tube*nr)*sn,tube*nz]
        if second: a=[a[0],a[2],a[1]];b=[b[0],b[2],b[1]];p=[1+p[0],p[2],p[1]]
        dp=[[gain1*a[i]*a[j]+gain2*b[i]*b[j] for j in range(3)] for i in range(3)]
        return p,dp,(length-tube)**2
    proj,dp,da=frame(False)
    if raw['second_ring']:
        q,dq,db=frame(True)
        if 'mix_temperature' in raw:
            temp=exactfloat(raw['mix_temperature'],ctx)
            weight=1/(1+ctx.exp((da-db)/temp))
            delta=[a-b for a,b in zip(proj,q)]
            dp=[[weight*dp[i][j]+(1-weight)*dq[i][j]+2*weight*(1-weight)*delta[i]*delta[j]/temp for j in range(3)] for i in range(3)]
            proj=[weight*a+(1-weight)*b for a,b in zip(proj,q)]
        elif db<da:proj,dp=q,dq
        elif not da<db:raise ArithmeticError('Unresolved branch switch')
    ident=eye(ctx)
    derivative=matmul(dp,[[2*cov[i][j]-ident[i][j] for j in range(3)] for i in range(3)])
    jac=[[beta*(ident[i][j]+derivative[i][j])+(1-2*beta)*cov[i][j] for j in range(3)] for i in range(3)]
    return [beta*(a+b)+(1-2*beta)*c for a,b,c in zip(x,proj,mean)],jac


def iterate(x,period,rec,ctx):
    matrix=eye(ctx);orbit=[]
    for _ in range(period):
        orbit.append(x)
        x,j=step(x,rec,ctx)
        matrix=matmul(j,matrix)
    return x,matrix,orbit
