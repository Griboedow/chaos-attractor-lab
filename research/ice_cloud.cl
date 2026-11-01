/* Literal logarithmic three-state ice-cloud field, piecewise analytic tangent, adaptive DP5(4). */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double16 rhs(double16 z,double16 p){
    double u=z[0],v=z[1],s=z[2];
    double a=p[0],d=p[1],e=p[2],c=p[3],f=p[4],g=p[5],h=p[6],i=p[7],k=p[8],sref=p[9];
    double en=exp(k*(s-sref)),N=a*en*exp(-u),M=d*en*exp(-v),H=h*en;
    double G=e*exp(2.*(u-v)/3.),C=c*exp(2.*(v-u)/3.),F=f*exp(2.*(v-u)/3.);
    double I=i*exp((2.*u+v)/3.),r=s-1.,rm=fmin(r,0.);
    double16 ans=(double16)(0.);
    ans[0]=N+G*rm-C;ans[1]=M+G*r-F;ans[2]=g*s-H-I*r;
    ans[3]=(-N+2.*(G*rm+C)/3.)*z[3]-2.*(G*rm+C)/3.*z[4]+(k*N+G*(s<1.))*z[5];
    ans[4]=2.*(G*r+F)/3.*z[3]+(-M-2.*(G*r+F)/3.)*z[4]+(k*M+G)*z[5];
    ans[5]=-2.*I*r/3.*z[3]-I*r/3.*z[4]+(g-k*H-I)*z[5];
    return ans;
}
__kernel void evaluate(__global const double16 *s,__global const double16 *p,__global double16 *f){
    int i=get_global_id(0);f[i]=rhs(s[i],p[i]);
}
__kernel void advance(__global double16 *states,__global const double16 *parameters,
    __global double2 *clocks,__global double4 *logs,__global double16 *ranges,
    __global int *status,__global int2 *counts,int attempts,double warm,double duration,
    double rtol,double atol,double maxstep,int max_attempts){
    int id=get_global_id(0);if(status[id])return;
    double16 s=states[id],p=parameters[id],range=ranges[id];double4 l=logs[id];double2 clock=clocks[id];
    int2 cnt=counts[id];double t=clock.x,h=clock.y;
    for(int j=0;j<attempts&&t<warm+duration;j++){
        if(cnt.x+cnt.y>=max_attempts){status[id]=3;break;}
        int w=t<warm?-1:min(3,(int)floor((t-warm)/(duration/4.)));
        double target=w<0?warm:warm+(w+1)*(duration/4.);
        while(target<=t&&w<3){w++;target=warm+(w+1)*(duration/4.);}
        h=fmin(h,fmin(maxstep,target-t));
        if(h<1e-14||t+h==t){status[id]=2;break;}
        double16 k1=rhs(s,p),k2=rhs(s+h*(1./5.*k1),p);
        double16 k3=rhs(s+h*(3./40.*k1+9./40.*k2),p);
        double16 k4=rhs(s+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double16 k5=rhs(s+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double16 k6=rhs(s+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double16 next=s+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double16 k7=rhs(next,p);
        double16 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        double norm0=0.,norm1=0.;for(int k=3;k<6;k++){norm0+=s[k]*s[k];norm1+=next[k]*next[k];}
        double tangent_scale=atol+rtol*sqrt(fmax(norm0,norm1));double err=0.;
        for(int k=0;k<6;k++){
            double scale=k<3?atol+rtol*fmax(fabs(s[k]),fabs(next[k])):tangent_scale;
            double scaled=fabs(error[k])/scale;
            if(!isfinite(next[k])||!isfinite(scaled)){err=INFINITY;break;}err=fmax(err,scaled);
        }
        if(err<=1.){
            s=next;t+=h;cnt.x++;if(target-t<1e-12*fmax(1.,target))t=target;
            double norm=sqrt(norm1);if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(w>=0)l[w]+=log(norm);for(int k=3;k<6;k++)s[k]/=norm;
            if(fabs(s[0])>650.||fabs(s[1])>650.||s[2]<0.||s[2]>5.){status[id]=1;break;}
            if(status[id])break;
            if(w>=0)for(int k=0;k<3;k++){range[2*k]=fmin(range[2*k],s[k]);range[2*k+1]=fmax(range[2*k+1],s[k]);}
        }else cnt.y++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    states[id]=s;clocks[id]=(double2)(t,h);logs[id]=l;ranges[id]=range;counts[id]=cnt;
}
