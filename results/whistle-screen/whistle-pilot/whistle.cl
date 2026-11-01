/* Literal six-state whistle field, piecewise analytic tangent, adaptive DP5(4). */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double16 rhs(double16 s,double16 p){
    double x=s[0],y=s[1],P=s[2],g=s[3],q=s[4],v=s[5];
    double A=p[0],R=p[1],B=p[2],C=p[3],D=p[4],G=p[5],H=p[6],J=p[7],U=p[8],W=p[9],F=p[10];
    double e=exp(G*fabs(P)),decay=J*e*(1.+H*fabs(x));
    double16 f=(double16)(0.);
    f[0]=A*(1.-P)-R*x-B*x*fabs(x);f[1]=A*P-R*y-B*y*fabs(y);
    f[2]=C*(x-y)-D*g;f[3]=v-decay*g;f[4]=v;f[5]=U*(1.-q*q)*v-W*q+F*P;
    f[6]=(-R-2.*B*fabs(x))*s[6]-A*s[8];
    f[7]=(-R-2.*B*fabs(y))*s[7]+A*s[8];
    f[8]=C*(s[6]-s[7])-D*s[9];
    f[9]=-g*J*e*H*sign(x)*s[6]-g*decay*G*sign(P)*s[8]-decay*s[9]+s[11];
    f[10]=s[11];f[11]=F*s[8]+(-2.*U*q*v-W)*s[10]+U*(1.-q*q)*s[11];
    return f;
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
        double norm0=0.,norm1=0.;for(int k=6;k<12;k++){norm0+=s[k]*s[k];norm1+=next[k]*next[k];}
        double tangent_scale=atol+rtol*sqrt(fmax(norm0,norm1));double err=0.;
        for(int k=0;k<12;k++){
            double scale=k<6?atol+rtol*fmax(fabs(s[k]),fabs(next[k])):tangent_scale;
            double scaled=fabs(error[k])/scale;
            if(!isfinite(next[k])||!isfinite(scaled)){err=INFINITY;break;}err=fmax(err,scaled);
        }
        if(err<=1.){
            s=next;t+=h;cnt.x++;if(target-t<1e-12*fmax(1.,target))t=target;
            double norm=sqrt(norm1);if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(w>=0)l[w]+=log(norm);for(int k=6;k<12;k++)s[k]/=norm;
            for(int k=0;k<6;k++)if(fabs(s[k])>1e6){status[id]=1;break;}
            if(status[id])break;
            if(w>=0)for(int k=0;k<6;k++){range[2*k]=fmin(range[2*k],s[k]);range[2*k+1]=fmax(range[2*k+1],s[k]);}
        }else cnt.y++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    states[id]=s;clocks[id]=(double2)(t,h);logs[id]=l;ranges[id]=range;counts[id]=cnt;
}
