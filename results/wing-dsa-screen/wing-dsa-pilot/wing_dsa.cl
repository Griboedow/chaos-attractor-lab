/* Literal normalized 2022 stretch-activation field and analytic tangent. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double8 rhs(double8 s,double4 p){
    double a=p.x*(1.+p.y),b=p.y*p.x*p.x,c=p.z*p.x*(1.-p.y);
    return p.w*(double8)(s.s1,-s.s0-s.s1*fabs(s.s1)+s.s2,s.s3,c*s.s1-b*s.s2-a*s.s3,
        s.s5,-s.s4-2.*fabs(s.s1)*s.s5+s.s6,s.s7,c*s.s5-b*s.s6-a*s.s7);
}
__kernel void evaluate(__global const double8 *s,__global const double4 *p,__global double8 *v){
    int id=get_global_id(0);v[id]=rhs(s[id],p[id]);
}
__kernel void advance(__global double8 *states,__global const double4 *parameters,
    __global double2 *clocks,__global double4 *logs,__global double8 *stats,
    __global int *status,__global int2 *counts,int attempts,double warm,double duration,
    double rtol,double atol,double maxstep,int max_attempts){
    int id=get_global_id(0);if(status[id])return;
    double8 s=states[id],st=stats[id];double4 p=parameters[id],l=logs[id];double2 clock=clocks[id];
    int2 cnt=counts[id];double t=clock.x,h=clock.y;
    for(int j=0;j<attempts&&t<warm+duration;j++){
        if(cnt.x+cnt.y>=max_attempts){status[id]=3;break;}
        int window=t<warm?-1:min(3,(int)floor((t-warm)/(duration/4.)));
        double target=window<0?warm:warm+(window+1)*(duration/4.);
        while(target<=t&&window<3){window++;target=warm+(window+1)*(duration/4.);}
        h=fmin(h,fmin(maxstep/p.w,target-t));
        if(h<1e-14||t+h==t){status[id]=2;break;}
        double8 k1=rhs(s,p),k2=rhs(s+h*(1./5.*k1),p);
        double8 k3=rhs(s+h*(3./40.*k1+9./40.*k2),p);
        double8 k4=rhs(s+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double8 k5=rhs(s+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double8 k6=rhs(s+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double8 next=s+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double8 k7=rhs(next,p);
        double8 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        double4 physical=fabs(error.lo)/(atol+rtol*fmax(fabs(s.lo),fabs(next.lo)));
        double4 tangent=fabs(error.hi)/(atol+rtol*fmax(length(s.hi),length(next.hi)));
        double err=0.;for(int z=0;z<4;z++){err=fmax(err,physical[z]);err=fmax(err,tangent[z]);}
        if(!all(isfinite(next))||!all(isfinite(physical))||!all(isfinite(tangent)))err=INFINITY;
        if(err<=1.){
            s=next;t+=h;cnt.x++;if(target-t<1e-12*fmax(1.,target))t=target;
            double norm=length(s.hi);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(window>=0)l[window]+=log(norm);s.hi/=norm;
            for(int z=0;z<4;z++)if(fabs(s[z])>1e7){status[id]=1;break;}
            if(status[id])break;
            if(window>=0)for(int z=0;z<4;z++){st[2*z]=fmin(st[2*z],s[z]);st[2*z+1]=fmax(st[2*z+1],s[z]);}
        }else cnt.y++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    clocks[id]=(double2)(t,h);states[id]=s;logs[id]=l;stats[id]=st;counts[id]=cnt;
}
