function PeriodicIsland2
% Doubly periodic Stokes flow with "island" geom using circle of SLP proxy 
% sources. Gary Marple using code from Alex Barnett 6/5/15.

clear all
warning off
addpath(genpath('.'))
uc.d = 2*pi; % unitcell, d=period in space
[I,N]=MakeI;
MI=length(I);
for k=1:MI
    I{k}=quadr(I{k});
end
mu = 0.7; % viscosity
uc.nei = 1; % how many nei copies either side (use eg 1e3 to test A w/o AP)
n = 32; % pts per side
[x,w] = gauss(n); x = (1+x)/2; w = w'/2; % quadr on [0,1]
U.x=uc.d*x-uc.d/2+1i*uc.d/2; U.w=uc.d*w; U.nx=0*U.x-1i;
D=U; D.x=U.x-1i*uc.d;
H = uc.d; L.x = -uc.d/2+ 1i*H*x-uc.d/2*1i; L.nx = 0*L.x+1; L.w = H*w; % left side
R = L; R.x = L.x+uc.d; % right side
M = 128; % # proxy pts in periodizing basis (2 force comps per pt, so 2M dofs)
b.x = 1.1*uc.d*exp(2i*pi*(1:M)'/M);        % the proxy pts
%b.nx = exp(2i*pi*(1:M)'/M);                % only needed if DLP proxy basis
b.w = 1+0*b.x; % unit quadr weights (dummy)

% computes f
h=.2; ue = @(x) h*[imag(x).^2;0*x]; % horiz Poisseuil flow (pres drop)
if 1
    f=zeros(2*sum(N),1); % No-slip BCs
else
    f=[];
    for k=1:MI
        f=[f;ue(I{k}.x)]; % Driving: I bdry vels, each stacked as [u_1;u_2]
    end
end

% computes g
alpha=1; % traction jump in x-direction
beta=0; % traction jump in y-direction
g=[zeros(2*n,1);alpha*ones(n,1);zeros(4*n,1);beta*ones(n,1)];



s.x=zeros(sum(N),1);
ls=0;
h=0;
for k=1:MI
   s.x(ls+1:ls+N(k))=I{k}.x;
   ls=ls+N(k);
   h0=ArcLength(real(I{k}.x([end,1:end])),imag(I{k}.x([end,1:end])),2*pi)/N(k);
   h=max(h0,h);
end
s.len=N;
s=Quad(s);
t.x=zeros(8*sum(N),1);
ls=0;
for j=-1:1
    for i=-1:1
        if i~=0||j~=0
            for k=1:MI
            t.x(ls+1:ls+s.len(k))=I{k}.x+i*uc.d+j*1i*uc.d;
            ls=ls+s.len(k);
            end
        end
    end
end


% SM=SLP_Self_Matrix(s);
% SDMinv=SM;
% for k=1:MI
%    SDMinv{k}=inv(eye(2*N(k))/2+SM{k}/mu+DLPmatrix(I{k},I{k},mu));
% end

%Aop=@(x)A_SD_operator(x,I,N,s,SM,t,h,uc,mu);
Aop=@(x)A_SD_operator(x,I,N,s,[],t,h,uc,mu);
Precond=@(x)PreCond(x,I,N,[],mu);


% Aop(RM)
% [Aop(RM),Precond(RM)]
%fill system matrices A B C Q (subblock ordering always U then D), same
% A2 = eye(2*sum(N))/2; % A's jump relation part.
% ls=0;
% for k=1:MI
%     for j=-1:1
%         for i=-1:1
%             I0=I{k};
%             I0.x=I{k}.x+uc.d*1i*j+uc.d*i;
%             A2(ls+1:ls+2*N(k),ls+1:ls+2*N(k))=A2(ls+1:ls+2*N(k),ls+1:ls+2*N(k))+SLPmatrix_scf(I{k}.x,I0,mu)+DLPmatrix_scf(I{k}.x,I0,mu);
%             for kk=[1:k-1,k+1:MI]
%                 I1.x=I{kk}.x+uc.d*1i*j+uc.d*i;
%                 I1=quadr(I1);
%                 lk=2*sum(N(1:kk-1));
%                 A2(ls+1:ls+2*N(k),lk+1:lk+2*N(kk))=...
%                     A2(ls+1:ls+2*N(k),lk+1:lk+2*N(kk))+...
%                     SLPmatrix_scf(I{k}.x,I1,mu)+DLPmatrix_scf(I{k}.x,I1,mu);
%             end
%         end
%     end
%    ls=ls+2*N(k);
% end

% single layer (monopoles) on proxy points:

% B2=[];
% for k=1:MI
% B2 = [B2;SLPmatrix(I{k},b,mu)]; % maps 2M peri dofs to vels on I
% end

% B1=zeros(2*sum(N),2*M);
% ls=0;
% for k=1:MI
%    B1(ls+1:ls+2*N(k),:)=SLPmatrix(I{k},b,mu);
%    ls=ls+2*N(k);
% end
B=@(x)ApplyB(x,I,N,b,mu);

% RM=rand(2*M,1);
% max(max(abs(B*RM-ApplyB(RM,I,N,b,mu))))
% pause(100)

% % computes C
% C1=zeros(8*n,2*sum(N));
% ls=0;
% for k=1:MI
%     for j=-1:1
%         a=uc.nei*uc.d;
%         I0=I{k}; I0.x=I{k}.x+uc.d*1i*j; IU=I{k}; ID=I{k};
%         IU.x=I{k}.x+uc.d*j-uc.d*1i; ID.x=I{k}.x-uc.d*j+uc.d*1i;
%         [S_RI,S_RIn] = SLPmatrix(R,I0,mu,-a); [S_LI,S_LIn] = SLPmatrix(L,I0,mu,a);
%         [S_UI,S_UIn] = SLPmatrix(U,IU,mu); [S_DI,S_DIn] = SLPmatrix(D,ID,mu);
%         [D_RI,D_RIn] = DLPmatrix(R,I0,mu,-a); [D_LI,D_LIn] = DLPmatrix(L,I0,mu,a);
%         [D_UI,D_UIn] = DLPmatrix(U,IU,mu); [D_DI,D_DIn] = DLPmatrix(D,ID,mu);
%         C1(:,ls+1:ls+2*N(k)) = C1(:,ls+1:ls+2*N(k))+...
%             [S_RI-S_LI; S_RIn-S_LIn;S_UI-S_DI;S_UIn-S_DIn]+...
%             [D_RI-D_LI; D_RIn-D_LIn;D_UI-D_DI;D_UIn-D_DIn];
%         % maps cancelled densities to discrepancy
%     end
%     ls=ls+2*N(k);
% end

C=@(x)ApplyQdagCP(x,I,N,L,R,U,D,uc,mu);

% computes Q
[Rb,Rbn] = SLPmatrix(R,b,mu); [Lb,Lbn] = SLPmatrix(L,b,mu);
[Ub,Ubn] = SLPmatrix(U,b,mu); [Db,Dbn] = SLPmatrix(D,b,mu);
Q = [Rb-Lb; Rbn-Lbn;Ub-Db; Ubn-Dbn]; % maps periodizing dofs to discrepancy

[Z,W]=MakeProj(I,N);
% CZ=ApplyQdagCP(Z,I,N,L,R,U,D,uc,mu);
CZ=C(Z);

%CP=C-(C*Z)*W';
% AP=A-(A*Z)*W';
APop=@(x)Aop(x-Z*(W'*x));

tau1=Z*[-alpha;beta;0]*uc.d;

% % we form the new RHS
f1=f-Aop(tau1);
g1=g-C(tau1);

% % we form Q\CP and Q\g1
%QdagCP=Q\CP;
QdagCP=@(x)ApplyQdagCP(x,I,N,L,R,U,D,uc,mu,Q,CZ,W);
Qdagg1=Q\g1;

S=@(x)(APop(x)-B((QdagCP(x))));
%tau2=(AP-B*QdagCP)\(f1-B*Qdagg1);
tic,
tau2=gmres(S,f1-B(Qdagg1),[],eps,150,Precond);
toc
cc=Qdagg1-QdagCP(tau2);
tau2=tau2-Z*(W'*tau2);
tautau=tau1+tau2;
G=[tautau;cc];

tau=cell(1,MI);
ls=0;
for k=1:MI
    tau{k}=G(ls+1:ls+2*N(k));
    ls=ls+2*N(k);
end
c=G(ls+1:end);

% we compute tau and c
% tau=cell(1,MI);
% ls=0;
% for k=1:MI
% tau{k}=tau1(ls+1:ls+2*N(k))+tau2(ls+1:ls+2*N(k));
% ls=ls+2*N(k);
% end
% c=c2;

% displays residual for linear system [A,B;C,Q]
display(['Linear system has residual ',num2str(max(max(abs([Aop(G(1:2*sum(N)))+B(G(2*sum(N)+1:end))-f;C(G(1:2*sum(N)))+Q*G(2*sum(N)+1:end)-g]))))])
% display(['Linear system has residual ',num2str(max(max(abs([A,B;C,Q]*G-[f;g]))))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting is done starting here (you can comment all of this)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generates grid points
ta=linspace(-uc.d/2,uc.d/2,128)';
ta=ta(2:end);
[XX,YY]=meshgrid(ta);

% evaluates target points
z.x=reshape(XX,size(XX,1)*size(XX,2),1)+1i*reshape(YY,size(YY,1)*size(YY,2),1);

% Nodes of islands
Nodes=zeros(sum(N),2);
ls=0;
for k=1:MI
    Nodes(ls+1:ls+N(k),:)=[real(I{k}.x),imag(I{k}.x)];
    ls=ls+N(k);
end

% Edges of islands
Edges=zeros(sum(N),2);
ls=0;
for k=1:MI
    Edges(ls+1:ls+N(k),:)=[ls+1:ls+N(k);ls+2:ls+N(k),ls+1]';
    ls=ls+N(k);
end

% Determines which grid points are inside the islands
in=inpoly([real(z.x),imag(z.x)],Nodes,Edges);


% This part computes the velocity at the grid points outside of the islands

zn=size(z.x,1);
u=zeros(2*zn,1);
maxn=1e5;

z.x=z.x(~in);
u0=u([~in;~in]);

zn=size(z.x,1);
lmin=0;
lmax=min([maxn,zn]);

for g=1:ceil(zn/maxn)
    [g;ceil(zn/maxn)]
    lind=[lmin+1:lmax,zn+lmin+1:zn+lmax]';
    z0.x=z.x(lind(1:end/2));
    
    u0(lind) = SLPmatrix(z0,b,mu)*c; % eval @ test pt: first do MFS (proxy) contrib
%     for k=1:MI
%         for j=-1:1
%             for i=-1:1
%                 % Close evaluation
%                 I0=I{k};
%                 I0.x=I{k}.x+uc.d*i+uc.d*1i*j;
%                 
%                 % Trapezoidal Rule
%                 u0(lind) = u0(lind)+(DLPmatrix_scf(z0.x,I0,mu)+SLPmatrix_scf(z0.x,I0,mu))*tau{k};
%             end
%         end
%     end
    z00.x=[];
    for j=-1:1
        for i=-1:1
            z00.x=[z00.x;z0.x+uc.d*i+uc.d*1i*j;];
        end
    end
    
%     size(SD_operator(G(1:2*sum(N)),I,N,s,z00,h,mu))
%     size(u0(lind))
    u0(lind) = u0(lind)+SD_operator(G(1:2*sum(N)),I,N,s,z00,h,mu);
    
    lmin=lmax;
    lmax=min([lmin+maxn,zn]);
end

u([~in;~in])=u0;
%u([in;in])=nan;
u([in;in])=0;

VV1=reshape(u(1:end/2),size(XX,1),size(XX,2)); % x-comp of vel
VV2=reshape(u(end/2+1:end),size(XX,1),size(XX,2)); % y-comp of vel


% Plots results ----

% quiver(XX,YY,VV1,VV2)
pcolor(XX,YY,sqrt(VV1.^2+VV2.^2))
shading flat
colormap(jet(1000))
hold on
for k=1:MI
    fill(real([I{k}.x(end);I{k}.x]),imag([I{k}.x(end);I{k}.x]),[0.8 0.8 0.8])
end
axis equal
axis([-uc.d/2,uc.d/2,-uc.d/2,uc.d/2])

% % % % % surface plot of x-comp of vel
% figure
% surface(XX,YY,VV1)
% title('x-component of velocity')
% xlabel('x')
% ylabel('y')
% axis([-uc.d/2,uc.d/2,-uc.d/2,uc.d/2,min(min(VV1)),max(max(VV1))])
% grid on
% 
% % surface plot of y-comp of vel
% figure
% surface(XX,YY,VV2)
% title('y-component of velocity')
% xlabel('x')
% ylabel('y')
% zlabel('u_y')
% axis([-uc.d/2,uc.d/2,-uc.d/2,uc.d/2,min(min(VV2)),max(max(VV2))])
% grid on
% end