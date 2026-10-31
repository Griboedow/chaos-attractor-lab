#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#ifndef N
#define N 4
#endif
#define NP (1+4*(N-1))
#define ND (2*N)

// Fixed-size Cholesky, local bead resistance is positive definite.
inline void factor(const double *m,double *L) {
    for(int i=0;i<N*N;i++)L[i]=0.;
    for(int i=0;i<N;i++)for(int j=0;j<=i;j++) {
        double v=m[i*N+j];for(int k=0;k<j;k++)v-=L[i*N+k]*L[j*N+k];
        L[i*N+j]=(i==j)?sqrt(v):v/L[j*N+j];
    }
}
inline void solve(const double *L,const double *v,double *r) {
    double y[N];
    for(int i=0;i<N;i++){double a=v[i];for(int j=0;j<i;j++)a-=L[i*N+j]*y[j];y[i]=a/L[i*N+i];}
    for(int i=N-1;i>=0;i--){double a=y[i];for(int j=i+1;j<N;j++)a-=L[j*N+i]*r[j];r[i]=a/L[i*N+i];}
}
inline void rhs(const double *z,const double *p,double *out) {
    double l[N],drag[N],spring[N],rest[N],tail[N];
    l[0]=1.;drag[0]=1.;spring[0]=1.;rest[0]=0.;
    for(int i=1;i<N;i++) {
        spring[i]=p[i];drag[i]=p[(N-1)+i];l[i]=p[2*(N-1)+i];
        rest[i]=rest[i-1]+p[3*(N-1)+i];
    }
    tail[N-1]=drag[N-1];for(int i=N-2;i>=0;i--)tail[i]=tail[i+1]+drag[i];
    double m[N*N],dm[N*N],L[N*N],g[N],dg[N],f[N],df[N];
    for(int i=0;i<N;i++) {
        g[i]=p[0]*l[i]*sin(z[i]-z[N-1]);
        dg[i]=p[0]*l[i]*cos(z[i]-z[N-1])*(z[N+i]-z[ND-1]);
        for(int j=0;j<N;j++) {
            double b=tail[max(i,j)]*l[i]*l[j],delta=z[i]-z[j];
            m[N*i+j]=b*cos(delta);dm[N*i+j]=-b*sin(delta)*(z[N+i]-z[N+j]);
            double k=0.;
            if(i==j)k=spring[i]+((i<N-1)?spring[i+1]:0.);
            if(i==j+1)k=-spring[i];if(j==i+1)k=-spring[j];
            g[i]-=k*(z[j]-rest[j]);dg[i]-=k*z[N+j];
        }
    }
    factor(m,L);solve(L,g,f);
    for(int i=0;i<N;i++)for(int j=0;j<N;j++)dg[i]-=dm[N*i+j]*f[j];
    solve(L,dg,df);for(int i=0;i<N;i++){out[i]=f[i];out[i+N]=df[i];}
}
__kernel void evaluate(__global const double *states,__global const double *params,__global double *values) {
    int id=get_global_id(0);double p[NP],z[ND],f[ND];
    for(int i=0;i<NP;i++)p[i]=params[id*NP+i];for(int i=0;i<ND;i++)z[i]=states[id*ND+i];
    rhs(z,p,f);for(int i=0;i<ND;i++)values[id*ND+i]=f[i];
}
__kernel void advance(__global double *states,__global const double *params,__global double *sums,
                      __global int *status,const double h,const int steps,const int qr) {
    int id=get_global_id(0);if(status[id])return;
    double p[NP],z[ND],a[ND],b[ND],c[ND],d[ND],v[ND],logs=0.;
    for(int i=0;i<NP;i++)p[i]=params[id*NP+i];for(int i=0;i<ND;i++)z[i]=states[id*ND+i];
    for(int step=0;step<steps;step++) {
        rhs(z,p,a);for(int i=0;i<ND;i++)v[i]=z[i]+.5*h*a[i];rhs(v,p,b);
        for(int i=0;i<ND;i++)v[i]=z[i]+.5*h*b[i];rhs(v,p,c);
        for(int i=0;i<ND;i++)v[i]=z[i]+h*c[i];rhs(v,p,d);
        for(int i=0;i<ND;i++)z[i]+=h*(a[i]+2.*b[i]+2.*c[i]+d[i])/6.;
        if((step+1)%qr==0) {
            double norm=0.,angle=0.;for(int i=0;i<N;i++){norm+=z[N+i]*z[N+i];angle+=fabs(z[i]);}
            norm=sqrt(norm);
            if(!isfinite(norm)||norm<=0.||!isfinite(angle)||angle>1.e6){status[id]=1;return;}
            logs+=log(norm);for(int i=N;i<ND;i++)z[i]/=norm;
        }
    }
    for(int i=0;i<ND;i++)states[id*ND+i]=z[i];sums[id]+=logs;
}
