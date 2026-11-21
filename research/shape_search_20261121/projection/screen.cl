#ifdef SMOOTH_LINK
typedef struct {
    float3 p, azimuth, meridian;
    float gain1, gain2, axis, centerline, dist2;
} torus_frame;
torus_frame frame_at(float3 x,int second,float tube) {
    if (second) x=(float3)(x.x-1,x.z,x.y);
    torus_frame f;
    f.axis=length(x.xy);
    float cs=x.x/f.axis, sn=x.y/f.axis;
    f.centerline=hypot(f.axis-1,x.z);
    float nr=(f.axis-1)/f.centerline, nz=x.z/f.centerline;
    f.gain1=(1+tube*nr)/f.axis; f.gain2=tube/f.centerline;
    f.dist2=(f.centerline-tube)*(f.centerline-tube);
    f.p=(float3)((1+tube*nr)*cs,(1+tube*nr)*sn,tube*nz);
    f.azimuth=(float3)(-sn,cs,0); f.meridian=(float3)(nz*cs,nz*sn,-nr);
    if (second) {
        f.p=f.p.xzy+(float3)(1,0,0);
        f.azimuth=f.azimuth.xzy; f.meridian=f.meridian.xzy;
    }
    return f;
}
float3 deriv(torus_frame f,float3 v) {
    return f.gain1*dot(f.azimuth,v)*f.azimuth + f.gain2*dot(f.meridian,v)*f.meridian;
}
float3 blend_deriv(torus_frame a,torus_frame b,float weight,float temp,float3 v) {
    float3 delta=a.p-b.p;
    return weight*deriv(a,v)+(1-weight)*deriv(b,v)
           + (2*weight*(1-weight)/temp)*dot(delta,v)*delta;
}
#endif
__kernel void screen(__global const float4 *data, __global float *out,
                     int stride, int burn, int window
#ifdef SAVE_CLOUD
                     , __global float4 *cloud, int saved_count
#endif
                     ) {
    const int i = get_global_id(0);
    float4 cfg = data[i];
    float beta=cfg.x, inv_tau=cfg.y, tube=cfg.z;
    int nsites=(int)cfg.w;
#ifdef SMOOTH_LINK
    float mix_temp=data[2*stride+i].w;
#endif
    float3 sites[8];
    float offsets[8];
    for (int j=0;j<nsites;j++) {
        sites[j]=data[(j+2)*stride+i].xyz;
        offsets[j]=-dot(sites[j],sites[j]);
    }
    float3 x=data[stride+i].xyz;
    float3 t=normalize((float3)(0.271f,0.534f,-0.801f));
    float sums[4]={0,0,0,0};
#ifdef SECOND_TANGENT
    float3 s=(float3)(0,1,0);
    float second_sums[4]={0,0,0,0};
#ifdef THIRD_TANGENT
    float third_sums[4]={0,0,0,0};
#endif
#endif
    float3 mean=0, m2=0, cross_moments=0;
    float min_axis=1e10f, min_tube=1e10f;
#ifdef LINKED
    int second_count=0;
    float min_switch=1e10f;
#endif
    int invalid=0;
    for (int step=0;step<burn+4*window;step++) {
        float w[8]; float maxi=-1e30f;
        for (int j=0;j<nsites;j++) {
            w[j]=(2*dot(x,sites[j])+offsets[j])*inv_tau;
            maxi=fmax(maxi,w[j]);
        }
        float weight=0; float3 m=0, raw=0;
#ifdef SECOND_TANGENT
        float3 raw2=0;
#ifdef THIRD_TANGENT
        float3 v3=normalize(cross(t,s)), raw3=0;
#endif
#endif
        for (int j=0;j<nsites;j++) {
            float a=exp(w[j]-maxi);
            float3 b=sites[j];
            weight+=a; m+=a*b; raw+=a*dot(b,t)*b;
#ifdef SECOND_TANGENT
            raw2+=a*dot(b,s)*b;
#ifdef THIRD_TANGENT
            raw3+=a*dot(b,v3)*b;
#endif
#endif
        }
        m/=weight;
        float3 dm=2*inv_tau*(raw/weight-m*dot(m,t));
        float3 v=2*m-x, dv=2*dm-t;
#ifdef SECOND_TANGENT
        float3 dm2=2*inv_tau*(raw2/weight-m*dot(m,s));
        float3 dv2=2*dm2-s;
#ifdef THIRD_TANGENT
        float3 dm3=2*inv_tau*(raw3/weight-m*dot(m,v3));
        float3 dv3=2*dm3-v3;
#endif
#endif
#ifdef SMOOTH_LINK
        torus_frame a=frame_at(v,0,tube), b=frame_at(v,1,tube);
        float mix=1/(1+exp(clamp((a.dist2-b.dist2)/mix_temp,-80.f,80.f)));
        float3 p=mix*a.p+(1-mix)*b.p;
        float3 dp=blend_deriv(a,b,mix,mix_temp,dv);
        float rho=fmin(a.axis,b.axis), l=fmin(a.centerline,b.centerline);
        if (step>=burn) {
            second_count+=(b.dist2<a.dist2);
            min_switch=fmin(min_switch,fabs(sqrt(a.dist2)-sqrt(b.dist2)));
        }
#ifdef SECOND_TANGENT
        float3 dp2=blend_deriv(a,b,mix,mix_temp,dv2);
#ifdef THIRD_TANGENT
        float3 dp3=blend_deriv(a,b,mix,mix_temp,dv3);
#endif
#endif
#else
#ifdef LINKED
        // Move the second circle (center [1,0,0] and xz plane) into
        // the first circle's coordinate frame to reuse its projector.
        float3 second=(float3)(v.x-1,v.z,v.y);
        float la=hypot(length(v.xy)-1,v.z);
        float lb=hypot(length(second.xy)-1,second.z);
        int use_second=lb<la;
        if (step>=burn) {
            second_count+=use_second;
            min_switch=fmin(min_switch,fabs(la-lb));
        }
        if (use_second) {v=second; dv=dv.xzy;}
#ifdef SECOND_TANGENT
        if (use_second) dv2=dv2.xzy;
#ifdef THIRD_TANGENT
        if (use_second) dv3=dv3.xzy;
#endif
#endif
#endif
        float rho=length(v.xy);
        float3 c=(float3)(v.x,v.y,0)/rho;
        float3 dc=((float3)(dv.x,dv.y,0)-c*dot(c,dv))/rho;
        float3 h=v-c;
        float l=length(h);
        float3 u=h/l, du=dv-dc;
        float3 p=c+tube*u, dp=dc+tube/l*(du-u*dot(u,du));
#ifdef SECOND_TANGENT
        float3 dc2=((float3)(dv2.x,dv2.y,0)-c*dot(c,dv2))/rho;
        float3 du2=dv2-dc2;
        float3 dp2=dc2+tube/l*(du2-u*dot(u,du2));
#ifdef THIRD_TANGENT
        float3 dc3=((float3)(dv3.x,dv3.y,0)-c*dot(c,dv3))/rho;
        float3 du3=dv3-dc3;
        float3 dp3=dc3+tube/l*(du3-u*dot(u,du3));
#endif
#endif
#ifdef LINKED
        if (use_second) {p=p.xzy+(float3)(1,0,0); dp=dp.xzy;}
#ifdef SECOND_TANGENT
        if (use_second) dp2=dp2.xzy;
#ifdef THIRD_TANGENT
        if (use_second) dp3=dp3.xzy;
#endif
#endif
#endif
#endif
        x=beta*(x+p)+(1-2*beta)*m;
        t=beta*(t+dp)+(1-2*beta)*dm;
        float k=length(t);
        if (!isfinite(k) || !isfinite(x.x) || rho<1e-7f || l<1e-7f) {invalid=1;break;}
        t=k>1e-12f?t/k:normalize((float3)(0.271f,0.534f,-0.801f));
#ifdef SECOND_TANGENT
        s=beta*(s+dp2)+(1-2*beta)*dm2;
        s-=dot(s,t)*t;
        s-=dot(s,t)*t;
        float k2=length(s);
        if (!isfinite(k2) || k2<1e-12f) {invalid=1;break;}
        s/=k2;
#ifdef THIRD_TANGENT
        // Residual third stretch after projecting onto the new normal.
        float3 next3=beta*(v3+dp3)+(1-2*beta)*dm3;
        float k3=fabs(dot(cross(t,s),next3));
        if (!isfinite(k3) || k3<1e-13f) {invalid=1;break;}
#endif
#endif
        if (step>=burn) {
            int j=step-burn;
            sums[j/window]+=log(fmax(k,1e-12f));
#ifdef SECOND_TANGENT
            second_sums[j/window]+=log(k2);
#ifdef THIRD_TANGENT
            third_sums[j/window]+=log(k3);
#endif
#endif
            float3 delta=x-mean;
            mean+=delta/(j+1);
            float3 tail=x-mean;
            m2+=delta*tail;
            cross_moments+=(float3)(delta.x*tail.y,delta.x*tail.z,delta.y*tail.z);
            min_axis=fmin(min_axis,rho); min_tube=fmin(min_tube,l);
#ifdef SAVE_CLOUD
            if (j>=4*window-saved_count) cloud[(j-(4*window-saved_count))*stride+i]=(float4)(x,0);
#endif
        }
    }
    for (int j=0;j<4;j++) out[j*stride+i]=invalid?NAN:sums[j]/window;
    int n=4*window;
    out[4*stride+i]=m2.x/n; out[5*stride+i]=m2.y/n; out[6*stride+i]=m2.z/n;
    out[7*stride+i]=cross_moments.x/n; out[8*stride+i]=cross_moments.y/n; out[9*stride+i]=cross_moments.z/n;
    out[10*stride+i]=min_axis; out[11*stride+i]=min_tube;
#ifdef LINKED
    out[12*stride+i]=(float)second_count/n; out[13*stride+i]=min_switch;
#endif
#ifdef SECOND_TANGENT
#ifdef LINKED
    int offset=14;
#else
    int offset=12;
#endif
    for (int j=0;j<4;j++) out[(offset+j)*stride+i]=invalid?NAN:second_sums[j]/window;
#ifdef THIRD_TANGENT
    for (int j=0;j<4;j++) out[(offset+4+j)*stride+i]=invalid?NAN:third_sums[j]/window;
#endif
#endif
}
