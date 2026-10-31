#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// p: force, k2,k3, drag2,drag3, length2,length3, rest_turn2,rest_turn3
inline void solve(const double *m, const double *v, double *r) {
    double a=m[0],b=m[1],c=m[2],d=m[4],e=m[5],f=m[8];
    double aa=d*f-e*e,bb=c*e-b*f,cc=b*e-c*d,dd=a*f-c*c,ee=b*c-a*e,ff=a*d-b*b;
    double det=a*aa+b*bb+c*cc;
    r[0]=(aa*v[0]+bb*v[1]+cc*v[2])/det;
    r[1]=(bb*v[0]+dd*v[1]+ee*v[2])/det;
    r[2]=(cc*v[0]+ee*v[1]+ff*v[2])/det;
}

inline void rhs(const double *z,const double *p,double *out) {
    double l[3]={1.,p[5],p[6]}, tail[3]={1.+p[3]+p[4],p[3]+p[4],p[4]};
    double rest[3]={0.,p[7],p[7]+p[8]};
    double k[9]={1.+p[1],-p[1],0.,-p[1],p[1]+p[2],-p[2],0.,-p[2],p[2]};
    double m[9],dm[9],g[3],dg[3],f[3],df[3];
    for(int i=0;i<3;i++) {
        g[i]=p[0]*l[i]*sin(z[i]-z[2]);
        dg[i]=p[0]*l[i]*cos(z[i]-z[2])*(z[3+i]-z[5]);
        for(int j=0;j<3;j++) {
            double b=tail[max(i,j)]*l[i]*l[j],delta=z[i]-z[j];
            m[3*i+j]=b*cos(delta);
            dm[3*i+j]=-b*sin(delta)*(z[3+i]-z[3+j]);
            g[i]-=k[3*i+j]*(z[j]-rest[j]);
            dg[i]-=k[3*i+j]*z[3+j];
        }
    }
    solve(m,g,f);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)dg[i]-=dm[3*i+j]*f[j];
    solve(m,dg,df);
    for(int i=0;i<3;i++){out[i]=f[i];out[i+3]=df[i];}
}

__kernel void evaluate(__global const double *states,__global const double *params,__global double *values) {
    int id=get_global_id(0); double p[9],z[6],f[6];
    for(int i=0;i<9;i++)p[i]=params[id*9+i];
    for(int i=0;i<6;i++)z[i]=states[id*6+i];
    rhs(z,p,f);for(int i=0;i<6;i++)values[id*6+i]=f[i];
}

__kernel void advance(__global double *states,__global const double *params,__global double *sums,
                      __global int *status,const double h,const int steps,const int qr) {
    int id=get_global_id(0);if(status[id])return;
    double p[9],z[6],a[6],b[6],c[6],d[6],v[6],logs=0.;
    for(int i=0;i<9;i++)p[i]=params[id*9+i];
    for(int i=0;i<6;i++)z[i]=states[id*6+i];
    for(int step=0;step<steps;step++) {
        rhs(z,p,a);
        for(int i=0;i<6;i++)v[i]=z[i]+.5*h*a[i];rhs(v,p,b);
        for(int i=0;i<6;i++)v[i]=z[i]+.5*h*b[i];rhs(v,p,c);
        for(int i=0;i<6;i++)v[i]=z[i]+h*c[i];rhs(v,p,d);
        for(int i=0;i<6;i++)z[i]+=h*(a[i]+2.*b[i]+2.*c[i]+d[i])/6.;
        if((step+1)%qr==0) {
            double norm=sqrt(z[3]*z[3]+z[4]*z[4]+z[5]*z[5]);
            if(!isfinite(norm)||norm<=0.||!isfinite(z[0]+z[1]+z[2])||fabs(z[0])+fabs(z[1])+fabs(z[2])>1.e6){status[id]=1;return;}
            logs+=log(norm);for(int i=3;i<6;i++)z[i]/=norm;
        }
    }
    for(int i=0;i<6;i++)states[id*6+i]=z[i];
    sums[id]+=logs;
}
