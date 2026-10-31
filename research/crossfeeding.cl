/* Ross et al. 2025 full author chemostat. Literal min branches, tau=D*t.
   q[0:5] physical state; q[5:10] tangent; q[10] integrated trace. */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
inline double16 rhs(double16 q,double16 p){
    double a=q.s0,b=q.s1,u=q.s2,v=q.s3,g=q.s4,D=p.sc;
    double m12=p.s0*v/(p.s2+v),m13=p.s0*g/(p.s1+g);
    double m21=p.s6*u/(p.s8+u),m23=p.s6*g/(p.s7+g);
    double d12=p.s0*p.s2/pow(p.s2+v,2),d13=p.s0*p.s1/pow(p.s1+g,2);
    double d21=p.s6*p.s8/pow(p.s8+u,2),d23=p.s6*p.s7/pow(p.s7+g,2);
    double m1=fmin(m12,m13),m2=fmin(m21,m23);
    double dm13=d13*q.s9,dm23=d23*q.s9;
    double dm1=m12<m13?d12*q.s8:dm13,dm2=m21<m23?d21*q.s7:dm23;
    double16 f=(double16)(0.);
    f.s0=(m1-D)*a;f.s1=(m2-D)*b;
    f.s2=D*(p.sd-u)+p.s5*(m13-m1)*a-p.sa*m2*b;
    f.s3=D*(p.se-v)+p.sb*(m23-m2)*b-p.s4*m1*a;
    f.s4=D*(p.sf-g)-p.s3*m13*a-p.s9*m23*b;
    f.s5=(m1-D)*q.s5+a*dm1;f.s6=(m2-D)*q.s6+b*dm2;
    f.s7=-D*q.s7+p.s5*((m13-m1)*q.s5+a*(dm13-dm1))-p.sa*(m2*q.s6+b*dm2);
    f.s8=-D*q.s8+p.sb*((m23-m2)*q.s6+b*(dm23-dm2))-p.s4*(m1*q.s5+a*dm1);
    f.s9=-D*q.s9-p.s3*(m13*q.s5+a*dm13)-p.s9*(m23*q.s6+b*dm23);
    f.sa=m1+m2-5.*D-p.sa*b*(m21<m23?d21:0.)-p.s4*a*(m12<m13?d12:0.)-p.s3*a*d13-p.s9*b*d23;
    return f/D;
}
__kernel void evaluate(__global const double16 *s,__global const double16 *p,__global double16 *f){
    int i=get_global_id(0);f[i]=rhs(s[i],p[i]);
}
__kernel void advance(__global double16 *state,__global const double16 *parameters,
    __global double2 *clock,__global double4 *logs,__global double16 *stats,
    __global int *status,__global int2 *counts,int work,
    double warm,double duration,double rtol,double atol,double maxstep,int budget){
    int id=get_global_id(0);if(status[id])return;
    double16 q=state[id],st=stats[id],p=parameters[id];double4 l=logs[id];
    double2 c=clock[id];int2 cnt=counts[id];double t=c.s0,h=c.s1;
    double w0=p.s3+p.s4*p.s9/p.sb,w1=p.s9+p.sa*p.s3/p.s5;
    double w2=p.s3/p.s5,w3=p.s9/p.sb,supply=w2*p.sd+w3*p.se+p.sf;
    for(int it=0;it<work&&t<warm+duration;it++){
        if(cnt.s0+cnt.s1>=budget){status[id]=3;break;}
        int win=t<warm?-1:min(3,(int)floor((t-warm)/(duration/4.)));
        double target=win<0?warm:warm+(win+1)*duration/4.;
        while(target<=t&&win<3){win++;target=warm+(win+1)*duration/4.;}
        h=fmin(h,fmin(maxstep,target-t));
        if(h<1e-14||t+h==t){status[id]=2;break;}
        double16 k1=rhs(q,p),k2=rhs(q+h/5.*k1,p);
        double16 k3=rhs(q+h*(3./40.*k1+9./40.*k2),p);
        double16 k4=rhs(q+h*(44./45.*k1-56./15.*k2+32./9.*k3),p);
        double16 k5=rhs(q+h*(19372./6561.*k1-25360./2187.*k2+64448./6561.*k3-212./729.*k4),p);
        double16 k6=rhs(q+h*(9017./3168.*k1-355./33.*k2+46732./5247.*k3+49./176.*k4-5103./18656.*k5),p);
        double16 next=q+h*(35./384.*k1+500./1113.*k3+125./192.*k4-2187./6784.*k5+11./84.*k6);
        double16 k7=rhs(next,p);
        double16 error=h*(71./57600.*k1-71./16695.*k3+71./1920.*k4-17253./339200.*k5+22./525.*k6-1./40.*k7);
        /* Tangents are unit-normalized: use their vector norm, not individual
           nearly zero components, to set a meaningful absolute error scale. */
        double oldnorm=0.,newnorm=0.;for(int j=5;j<10;j++){oldnorm+=q[j]*q[j];newnorm+=next[j]*next[j];}
        double tangent_scale=atol+rtol*sqrt(fmax(oldnorm,newnorm));
        double err=0.;for(int j=0;j<5;j++)err=fmax(err,fabs(error[j])/(atol+rtol*fmax(fabs(q[j]),fabs(next[j]))));
        for(int j=5;j<10;j++)err=fmax(err,fabs(error[j])/tangent_scale);
        if(!all(isfinite(next))||!isfinite(err))err=INFINITY;
        if(err<=1.){
            double divergence=next.sa;q=next;q.sa=0.;t+=h;cnt.s0++;
            if(target-t<1e-13*fmax(1.,target))t=target;
            double norm=0.;for(int j=5;j<10;j++)norm+=q[j]*q[j];norm=sqrt(norm);
            if(!(norm>0.)||!isfinite(norm)){status[id]=2;break;}
            if(win>=0)l[win]+=log(norm);
            for(int j=5;j<10;j++)q[j]/=norm;
            for(int j=0;j<5;j++){if(q[j]<-1e-8)status[id]=5;if(q[j]>1e6)status[id]=1;}
            if(status[id])break;
            if(fmin(q.s0,q.s1)<1e-18){status[id]=6;break;}
            double actual=w0*q.s0+w1*q.s1+w2*q.s2+w3*q.s3+q.s4;
            double expected=supply+(st.sc-supply)*exp(-t);
            st.sb=fmax(st.sb,fabs(actual-expected)/(1.+fabs(expected)));
            if(st.sb>1e-6){status[id]=7;break;}
            if(win>=0){for(int j=0;j<5;j++){st[2*j]=fmin(st[2*j],q[j]);st[2*j+1]=fmax(st[2*j+1],q[j]);}st.sa+=divergence;}
        }else cnt.s1++;
        h*=err==0.?5.:clamp(.9*pow(err,-.2),.2,5.);
    }
    state[id]=q;clock[id]=(double2)(t,h);logs[id]=l;stats[id]=st;counts[id]=cnt;
}
