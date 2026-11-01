"""Independent arbitrary-precision explicit cubic and directional derivative."""
from pathlib import Path
import hashlib,json,time
import mpmath as mp

def raw(w,e):
    A,B,C,D=w
    h=mp.matrix([A**3-3*A*B**2+18*e*A*A*C+e**2*(36*A*C*C-18*A*D*D)+18*e**3*C**3,
      B**3-3*B*A*A+18*e*B*B*D+e**2*(36*B*D*D-18*B*C*C)+18*e**3*D**3,
      6*e*A**3+e**2*(36*A*A*C-18*B*B*C)+54*e**3*A*C*C+e**4*(-18*C**3-108*C*D*D),
      6*e*B**3+e**2*(36*B*B*D-18*A*A*D)+54*e**3*B*D*D+e**4*(-18*D**3-108*D*C*C)])
    H=mp.matrix(4)
    H[0,0]=3*A*A-3*B*B+36*e*A*C+e**2*(36*C*C-18*D*D)
    H[1,1]=3*B*B-3*A*A+36*e*B*D+e**2*(36*D*D-18*C*C)
    H[2,2]=e**2*(36*A*A-18*B*B)+108*e**3*A*C-e**4*(54*C*C+108*D*D)
    H[3,3]=e**2*(36*B*B-18*A*A)+108*e**3*B*D-e**4*(54*D*D+108*C*C)
    for j,k,v in [(0,1,-6*A*B),(0,2,18*e*A*A+72*e**2*A*C+54*e**3*C*C),
                  (0,3,-36*e**2*A*D),(1,2,-36*e**2*B*C),
                  (1,3,18*e*B*B+72*e**2*B*D+54*e**3*D*D),(2,3,-216*e**4*C*D)]:H[j,k]=H[k,j]=v
    return h,H

def run(digits,epsilon):
    mp.mp.dps=digits;e=mp.mpf(epsilon);w=mp.matrix(list(map(mp.mpf,['.83','.29','-.13','.11'])));w/=mp.norm(w)
    v=mp.matrix([1,1,1,1]);v-=w*mp.fdot(w,v);v/=mp.norm(v);sums=[mp.mpf(0)]*5;mn=mp.inf;res=mp.mpf(0)
    # Independently differentiate the explicit vector at the initial state.
    h,H=raw(w,e);derivative=max(abs(H[j,k]-mp.diff(lambda x:raw(mp.matrix([x if l==k else w[l] for l in range(4)]),e)[0][j],w[k])) for j in range(4) for k in range(4))
    for it in range(22000):
        h,H=raw(w,e);hn=mp.norm(h);mn=min(mn,hn);f=h/hn;v=H*v;v=(v-f*mp.fdot(f,v))/hn;gain=mp.norm(v);v/=gain;w=f;res=max(res,abs(mp.fdot(w,w)-1))
        if it>=2000:sums[(it-2000)//4000]+=mp.log(gain)/4000
    return {'digits':digits,'epsilon_decimal_exact':epsilon,'warmup':2000,'observed':20000,'window':4000,'rates':list(map(str,sums)),'minimum_update_norm':str(mn),'sphere_residual':str(res),'initial_derivative_check':str(derivative),'input_semantics':'Exact decimal epsilon and initial state; explicit cubic and derivatives evaluated in mpmath.'}

if __name__=='__main__':
    rows=[];out=Path(__file__).with_name('precision.json')
    for digits,epsilon in [(70,'.04'),(100,'.04'),(70,'.15'),(100,'.15')]:
        t=time.time();row=run(digits,epsilon);row.update(wall_seconds=time.time()-t,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),mpmath_version=mp.__version__);rows.append(row);out.write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
