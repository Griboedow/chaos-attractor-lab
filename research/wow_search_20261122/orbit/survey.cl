/* Saturable dispersive cavity with an independently relaxing detuning. */
__kernel void survey(__global const float *pars, __global float *out,
                     __global float4 *cloud, int count, int burn, int win, int kept) {
    int j=get_global_id(0);
    const float p=pars[j*6], a=pars[j*6+1], d=pars[j*6+2];
    const float k=pars[j*6+3], gain=pars[j*6+4], h=pars[j*6+5];
    float x=.173f,y=-.341f,z=.271f;
    float3 v=(float3)(.26726124f,.53452248f,.80178374f);
    float rates[4]={0,0,0,0}, mn=1e10f,mx=-1e10f;
    float3 sum=(float3)(0),sq=(float3)(0),cross=(float3)(0);
    float shortsum=0,lows=0,highs=0;
    int total=burn+4*win;
    for(int i=0;i<total;i++){
        float I=x*x+y*y, recip=1.f/(1.f+I), f=I*recip;
        float sine,cosine; sine=sincos(d+k*f+z,&cosine);
        float ux=a*(cosine*x-sine*y),uy=a*(sine*x+cosine*y);
        float df=2.f*recip*recip*(x*v.x+y*v.y),angle=k*df+v.z;
        float3 vv=(float3)(a*(cosine*v.x-sine*v.y)-uy*angle,a*(sine*v.x+cosine*v.y)+ux*angle,(1.f-h)*v.z+h*gain*df);
        float mag=length(vv);v=vv/mag;
        float zz=(1.f-h)*z+h*gain*f;
        x=p+ux;y=uy;z=zz;
        if(i>=burn){
            int l=i-burn;float stretch=log(mag);rates[l/win]+=stretch;
            float3 u=(float3)(x,y,z);
            sum+=u;sq+=u*u;cross+=(float3)(x*y,x*z,y*z);mn=fmin(mn,z);mx=fmax(mx,z);
            shortsum+=stretch;
            if((l%24)==23){lows+=(shortsum<-.5f);highs+=(shortsum>1.5f);shortsum=0;}
        }
        if(kept>0 && i>=total-kept)cloud[(i-total+kept)*count+j]=(float4)(x,y,z,0);
    }
    float denom=4*win;
    for(int i=0;i<4;i++)out[j*16+i]=rates[i]/win;
    float3 m=sum/denom,var=sq/denom-m*m,cv=cross/denom-(float3)(m.x*m.y,m.x*m.z,m.y*m.z);
    out[j*16+4]=m.x;out[j*16+5]=m.y;out[j*16+6]=m.z;
    out[j*16+7]=var.x;out[j*16+8]=var.y;out[j*16+9]=var.z;
    out[j*16+10]=cv.x;out[j*16+11]=cv.y;out[j*16+12]=cv.z;
    out[j*16+13]=mx-mn;out[j*16+14]=lows/denom*24;out[j*16+15]=highs/denom*24;
}
