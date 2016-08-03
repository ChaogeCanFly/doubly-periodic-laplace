function TestLapInteractions
%% SLP test
clear
clf
rng(1)
addpath(genpath('..'))
n=64;
alpha=linspace(0,2*pi,n+1)';
alpha=alpha(2:end);
s.x=(cos(alpha)+1i*(sin(alpha)+0.5*cos(2*alpha)))*1e-2;
s.len(1)=n;
% s.x=[s.x;cos(alpha)+1i*exp(sin(alpha))+0.23*1i];
% s.len(2)=n;
s=Quad(s);
f=[exp(sin(alpha));exp(cos(alpha))];
f=f(1:end/2);
% t.x=(4*pi*rand(40,2)-2*pi)*[1;1i];
% t.x=cos(alpha)+1i*exp(sin(alpha))+0.23*1i;
plot(real(s.x),imag(s.x))
hold on
% scatter(real(t.x),imag(t.x),'r')
% t.x=s.x+.2288708980178568*s.nx;
t.x=s.x+1e-14*s.nx;
t.nx=s.nx;
NPt=FindNearPt(s,1,t,1,5*s.h);
scatter(real(t.x),imag(t.x),'g')
%U=Lap_SLP_Close_Apply(NPt,s,t,f);
hold on
scatter(real(t.x(NPt.t{1})),imag(t.x(NPt.t{1})),'b');
[u,ux,uy]=lapSevalclose(t.x,s,f,'e');


U_Far=Lap_SLP_Far_Apply(s,1,t,1,f);
U_Close=Lap_SLP_Close_Apply(NPt,s,t,f);
S=Lap_SLP_Self_Matrix(s);
% U_Self_Far=Lap_SLP_Far_Apply(s,1,s,1,f);
U_Self=Lap_SLP_Self_Apply(s,f,S);
SLP=LapSLPmatrix(s,s,0);
SLP(diag(true(size(SLP,1),1)))=0;

clc
disp('Use FMM code to compute Laplace potentials.')
disp('Accuracy of self eval:'),disp(max(abs(U_Self.s+U_Far.s-S{1}*f)))
disp('Accuracy of self eval2:'),disp(max(abs(U_Self.s+SLP*f-S{1}*f)))
disp('Accuracy of target eval:'),disp(max(abs(U_Close.t+U_Far.t-u)))
disp('Accuracy of target x partial der:'),disp(max(abs(U_Far.tx+U_Close.tx-ux)))
disp('Accuracy of target y partial der:'),disp(max(abs(U_Far.ty+U_Close.ty-uy)))

% test self interaction
[u1,ux1,uy1]=lapSevalclose(s.x*(1+1e-15),s,f,'e'); % exterior closeeval not accurate for self, thus perturb target pt a little
[u2,ux2,uy2]=lapSevalclose(s.x,s,f,'i');

disp('Use closeeval to approx. int&ext limits of u and du/dn')
disp('diff in |u_in-u_ex| ='),disp(max(abs(u1 - u2)))
disp('diff in |u_in-u_self| ='),disp(max(abs(u2 - LapSLPselfmatrix(s)*f)))
disp('diff in |du/dn_in-du/dn_ex| - |density| (jump relation):'), disp(max(abs(real((ux2-1i*uy2).*s.nx)-real((ux1-1i*uy1).*s.nx))-abs(f)))

%%
SM=Lap_SLP_Far_Matrix(s,1,t,1);
[SM.sn*f-S.sn]


%% DLP test
clc
clear
clf
rng(1)
addpath(genpath('..'))
n=256;
alpha=linspace(0,2*pi,n+1)';
alpha=alpha(2:end);
s.x=(cos(alpha)+1i*(sin(alpha)+0.5*cos(2*alpha)))*1e-1;
s.len(1)=n;
% s.x=[s.x;cos(alpha)+1i*exp(sin(alpha))+0.23*1i];
% s.len(2)=n;
s=Quad(s);
f=[exp(sin(alpha));exp(cos(alpha))];
f=f(1:end/2);
plot(real(s.x),imag(s.x))
hold on
% t.x=(4*pi*rand(40,2)-2*pi)*[1;1i];
% t.x=s.x+.2288708980178568/9.999999999999955*s.nx;
t.x=s.x+1e-4*s.nx;
t.nx=s.nx;
NPt=FindNearPt(s,1,t,1,5*s.h);
scatter(real(t.x),imag(t.x),'g')
hold on
scatter(real(t.x(NPt.t{1})),imag(t.x(NPt.t{1})),'b');
axis([-1,1,-1,1])

%1 source to target

%1.1 u = DLPmatrix eval
    [u,ux,uy]=lapDevalclose(t.x,s,f,'e');

%1.2 U.pot = DLP_Far_Apply
    %There are two ways:(1) use rfmm/lfmm, (2) use zfmm DO(2)!
    U_Far = Lap_DLP_Far_Apply(s,1,t,1,f);
    U_Close=Lap_DLP_Close_Apply(NPt,s,t,f);
    U_Close_demo=Lap_DLP_Close_Apply_demo(NPt,s,t,f);
    D=Lap_DLP_Self_Matrix(s);
    
%1.3 compare the two
    
    disp('Use FMM code to compute Laplace potentials.')
    disp('Accuracy of self eval:'),disp(max(abs(U_Far.s-D{1}*f)))
    disp('Accuracy of source norm der:'),disp(max(abs(U_Far.sn+U_Close.sn-real((U_Far.sx-1i*U_Far.sy).*s.nx))))
    disp('Accuracy of target eval:'),disp(max(abs(U_Far.t+U_Close.t-u)))
    disp('Accuracy2 of target eval:'),disp(max(abs((U_Far.t-U_Close_demo.sub.t)+U_Close_demo.add.t-u)))
    disp('Accuracy of target y partial der:'),disp(max(abs(U_Far.ty+U_Close.ty-uy)))
    disp('Accuracy2 of target y partial der:'),disp(max(abs((U_Far.ty-U_Close_demo.sub.ty)+U_Close_demo.add.ty-uy)))
    disp('Diff in target y partial der:'),disp(max(abs(U_Far.ty-U_Close_demo.sub.ty)))
    disp('Accuracy of target x partial der:'),disp(max(abs(U_Far.tx+U_Close.tx-ux)))
    disp('Accuracy2 of target x partial der:'),disp(max(abs((U_Far.tx-U_Close_demo.sub.tx)+U_Close_demo.add.tx-ux)))
    disp('Diff in target x partial der:'),disp(max(abs(U_Far.tx-U_Close_demo.sub.tx)))
    %[U_Far.ty-U_Close_demo.sub.ty,U_Close_demo.add.ty-uy]
    %uy

%2 source to self


% U.sn
% [U.tx,U.ty]
% pause(100)
% size(s.x)
% size(f)
% size(t.x)
% 
% t.x=s.x;
% t.nx=s.nx*1i
% [A,An]=LapSLPmatrix(t,s,1);
% figure
% surf(An)
% pause(100)
% 
% [~,ux,uy]=lapDevalclose(t.x,s,f,'e');
% [Uy.t,uy]
% V=Lap_DLP_Close_Apply(NPt,s,t,f);
% 
% 
% %S=Lap_SLP_Self_Matrix(s);
% %U_self=Lap_DLP_Self_Apply(s,f,S);
% U_far_M=Lap_SLP_Far_Matrix(s,1,t,1);
% U_far=Lap_SLP_Far_Apply(s,1,t,1,f);
% U_close=Lap_DLP_Close_Apply(NPt,s,t,f);
% ls=0;
% U_exact=zeros(length(t.x),1);
% for k=1:length(s.len)
% ss.x=s.x(ls+1:ls+s.len(k));
% ss=quadr(ss);
% %tt.x=s.x+1e2*eps*s.nx;
% U_exact=U_exact+lapDevalclose(t.x,ss,f(ls+1:ls+s.len(k)),'e'); 
% ls=ls+s.len(k);
% end
% 
% U_far_M.t*f-U_far.t
% U=U_far.t+U_close.t;
% %plot(U,'b')
% 
% % 
% k=1;
% ls=0;
% [~,ux,uy]=lapSevalclose(t.x,ss,f(ls+1:ls+s.len(k)),'e');
% [~,uux,uuy]=lapDevalclose(t.x,ss,f(ls+1:ls+s.len(k)),'e');
% t.nx=0*t.x+1;
% [~,Anx]=LapSLPmatrix(t,ss);
% [ux-Anx*f(ls+1:ls+s.len(k))]
% [~,AAnx]=LapDLPmatrix(t,ss);
% uux-AAnx*f(ls+1:ls+s.len(k))