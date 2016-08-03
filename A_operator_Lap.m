function u=A_operator(s,f,U,nei)

n=sum(s.len);
t.x=zeros(8*n,1);
t.nx=t.x;
ls=0;
for i=-nei:nei
    for j=-nei:nei
        if i||j
            t.x(ls+1:ls+n)=s.x+U.e1*i+U.e2*j;
            t.nx(ls+1:ls+n)=s.nx;
            ls=ls+n;
        end
    end
end


NPt=FindNearPt(s,1,t,1,5*max(s.h));
u_close=Lap_SLP_Close_Apply(NPt,s,t,f);
u_far=Lap_SLP_Far_Apply(s,1,t,1,f);

u=u_close.sn+u_far.sn;
ls=0;
for j=1:8
   u=u+(u_close.tx(ls+1:ls+n)+u_far.tx(ls+1:ls+n)).*real(s.nx)+(u_close.ty(ls+1:ls+n)+u_far.ty(ls+1:ls+n)).*imag(s.nx);
   ls=ls+n;
end