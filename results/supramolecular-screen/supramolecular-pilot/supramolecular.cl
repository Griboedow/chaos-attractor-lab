/* Sharko2022 simplified polymer moment closure, Eqs5.6/5.12.
   x=a1/kd2,y=d/kd2,z=m0/kd3,tau=kd3*t; k+=1. No clipping. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double8 rhs(double8 q,double4 p){
    double x=q.s0,y=q.s1,z=q.s2,alpha=p.s0,beta=p.s1,total=p.s2,C=p.s3;
    double mass=total-x-y,gap=beta*mass-2*z;
    double slope=1+4*beta*C*z*z*z/(gap*gap*gap);
    double dz=-6-2*beta*C*(3*z*z/(gap*gap)+4*z*z*z/(gap*gap*gap))-2*z;
    return (double8)(alpha*y-2*x*z,-alpha*y+2*z,beta*mass-6*z-2*beta*C*z*z*z/(gap*gap)-z*z,
        -2*z*q.s3+alpha*q.s4-2*x*q.s5,-alpha*q.s4+2*q.s5,
        -beta*slope*(q.s3+q.s4)+dz*q.s5,-2*z-alpha+dz,0.);
}
__kernel void evaluate(__global const double8 *s,__global const double4 *p,__global double8 *f){
    int i=get_global_id(0);f[i]=rhs(s[i],p[i]);
}
__kernel void advance(__global double8 *state,__global const double4 *parameters,
    __global double2 *clock,__global double4 *logs,__global double8 *stats,
    __global int *status,__global int2 *counts,int work,
    double warm,double duration,double rtol,double atol,double maxstep,int budget){
    int id=get_global_id(0);if(status[id])return;
    double8 q=state[id],st=stats[id];double4 p=parameters[id],l=logs[id];
    double2 c=clock[id];int2 cnt=counts[id];double t=c.s0,h=c.s1;
    for(int it=0;it<work&&t<warm+duration;it++){
        if(cnt.s0+cnt.s1>=budget){status[id]=3;break;}
        int win=t<warm?-1:min(3,(int)floor((t-warm)/(duration/4.)));
        double target=win<0?warm:warm+(win+1)*duration/4.;
        while(target<=t&&win<3){win++;target=warm+(win+1)*duration/4.;}
        h=fmin(h,fmin(maxstep,target-t));
        if(h<1e-14||t+h==t){status[id]=2;break;}
        double8 k1=rhs(q,p),k2=rhs(q+h/5.*k1,p);
        double8 k3=rhs(q+h*(3./40.*k1+9./40.*k2),p);
        double8 k4=rhs(q+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double8 k5=rhs(q+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double8 k6=rhs(q+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double8 next=q+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double8 k7=rhs(next,p);
        double8 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        double err=0.;for(int j=0;j<3;j++)err=fmax(err,fabs(error[j])/(atol+rtol*fmax(fabs(q[j]),fabs(next[j]))));
        double qnorm=sqrt(q.s3*q.s3+q.s4*q.s4+q.s5*q.s5),nnorm=sqrt(next.s3*next.s3+next.s4*next.s4+next.s5*next.s5);
        for(int j=3;j<6;j++)err=fmax(err,fabs(error[j])/(atol+rtol*fmax(qnorm,nnorm)));
        if(!all(isfinite(next))||!isfinite(err))err=INFINITY;
        if(err<=1.){
            double divergence=next.s6;q=next;q.s6=0.;q.s7=0.;t+=h;cnt.s0++;
            if(target-t<1e-13*fmax(1.,target))t=target;
            double norm=sqrt(q.s3*q.s3+q.s4*q.s4+q.s5*q.s5);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(win>=0)l[win]+=log(norm);
            q.s3/=norm;q.s4/=norm;q.s5/=norm;
            double gap=p.s1*(p.s2-q.s0-q.s1)-2*q.s2;
            if(fmin(q.s0,fmin(q.s1,q.s2)) < 0. || gap<=0.){status[id]=5;break;}
            if(q.s0>p.s2||q.s1>p.s2||q.s2>p.s1*p.s2/2){status[id]=1;break;}
            if(win>=0){for(int j=0;j<3;j++){st[2*j]=fmin(st[2*j],q[j]);st[2*j+1]=fmax(st[2*j+1],q[j]);}st.s6+=divergence;
                double length=p.s1*(p.s2-q.s0-q.s1)/q.s2;st.s7=st.s7==0.?length:fmin(st.s7,length);}
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    state[id]=q;clock[id]=(double2)(t,h);logs[id]=l;stats[id]=st;counts[id]=cnt;
}
