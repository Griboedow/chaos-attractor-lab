"""Two expressions for the posterior map through second derivatives."""
from product_binary_interval import model, parameters


def binary_second(x, par, ctx):
    beta,tau,c,a=par;y,J,r=model(x,par,ctx);m=[];d=[];e=[]
    for j in range(3):
        k=2*a[j]/tau;w=ctx.exp(2*k*(x[j]-c[j]));t=(w-1)/(w+1)
        m.append(c[j]+a[j]*t);d.append(4*a[j]*k*w/(w+1)**2);e.append(-2*k*t*d[-1])
    u=[(2*m[j]-x[j])/r for j in range(3)];v1=[2*q-1 for q in d]
    H=[[[beta*(-((i==j)*u[k]+(i==k)*u[j]+(j==k)*u[i]-3*u[i]*u[j]*u[k])/r**2*v1[j]*v1[k]
                 +((i==j)-u[i]*u[j])/r*2*e[j]*(j==k))+(1-2*beta)*e[i]*(i==j)*(i==k)
          for k in range(3)] for j in range(3)] for i in range(3)]
    return y,J,H,r


class Jet:
    """Value, gradient and Hessian with arithmetic in the supplied context."""
    def __init__(self,v,d=None,h=None):
        self.v=v;self.d=[v*0]*3 if d is None else d
        self.h=[[v*0]*3 for _ in range(3)] if h is None else h
    def __add__(self,b):
        if not isinstance(b,Jet):b=Jet(self.v*0+b)
        return Jet(self.v+b.v,[x+y for x,y in zip(self.d,b.d)],[[self.h[i][j]+b.h[i][j] for j in range(3)] for i in range(3)])
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,[-v for v in self.d],[[-v for v in row] for row in self.h])
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        if not isinstance(b,Jet):b=Jet(self.v*0+b)
        return Jet(self.v*b.v,[self.d[i]*b.v+self.v*b.d[i] for i in range(3)],
                   [[self.h[i][j]*b.v+self.d[i]*b.d[j]+self.d[j]*b.d[i]+self.v*b.h[i][j] for j in range(3)] for i in range(3)])
    __rmul__=__mul__
    def compose(self,value,first,second):
        return Jet(value,[first*d for d in self.d],[[second*self.d[i]*self.d[j]+first*self.h[i][j] for j in range(3)] for i in range(3)])
    def reciprocal(self):return self.compose(1/self.v,-1/self.v**2,2/self.v**3)
    def __truediv__(self,b):return self*(b.reciprocal() if isinstance(b,Jet) else 1/b)
    def exp(self,ctx):
        v=ctx.exp(self.v);return self.compose(v,v,v)
    def sqrt(self,ctx):
        v=ctx.sqrt(self.v);return self.compose(v,1/(2*v),-1/(4*v**3))


def full_second(x, raw, ctx):
    beta,tau,_,_=parameters(raw['params'],__import__('numpy').array(raw['points']),ctx)
    exact=lambda z:ctx.mpf(float(z).as_integer_ratio()[0])/float(z).as_integer_ratio()[1]
    points=[[exact(v) for v in row] for row in raw['points']]
    q=[Jet(v,[ctx.mpf(i==j) for j in range(3)]) for i,v in enumerate(x)]
    weights=[(-sum((q[j]-p[j])*(q[j]-p[j]) for j in range(3))/tau).exp(ctx) for p in points]
    total=sum(weights);m=[sum(w*p[j] for w,p in zip(weights,points))/total for j in range(3)]
    reflected=[2*m[j]-q[j] for j in range(3)];radius=sum(v*v for v in reflected).sqrt(ctx)
    reflected_sphere=[2*v/radius-v for v in reflected]
    y=[beta/2*(q[j]+reflected_sphere[j])+(1-beta)*m[j] for j in range(3)]
    return [v.v for v in y],[v.d for v in y],[v.h for v in y],radius.v


def determinant_gradient(J,H):
    from certify_product_cycles import det
    return [sum(det([[H[i][j][k] if row==i else J[row][j] for j in range(3)] for row in range(3)]) for i in range(3)) for k in range(3)]


def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
