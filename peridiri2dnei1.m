function peridiri2dnei1
% Solve 2D Neumann Laplace doubly-periodic BVP, via proxy periodizing
% Barnett for Zhao project 9/5/14. nei=1 for CSE figs 3/15/15. Schur 3/29/15

addpath(genpath('.'))
v = 1; % verbosity: 0 no plot; 1 potential plot; 2 diagnostic plots
nei = 1;      % 0 for no direct neighbors; 1 for preferred 3x3 scheme
schur = 1;    % 0 to solve expanded system, 1 to eliminate proxy (for GMRES)
jumps = [0 0]; %[1 0]; % potential jumps across R-L and T-B
M_num_ves = 5;
N_num_node = 64;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_e = 1;                           %%%
mu_i = 10;                          %%%
eta = (mu_i - mu_e)/(mu_i + mu_e);  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [I,N]=MakeI(10,M_num_ves,N_num_node); % create particles in the domain
% s.x=zeros(sum(N),1);
% s.len=N;
% ls=0;
% for k=1:length(N)
%     s.x(ls+1:ls+s.len(k))=I{k}.x/2/pi;
%     ls=ls+s.len(k);
% end
s = Make_s;

% Build the periodic box (L,R,T,B)
U.e1 = 1*2*pi; U.e2 = 1i*2*pi; U.nei = 0;  % unit cell
M = 40; [x,w] = gauss(M); w=w/2;  % walls
L.x = (-U.e1 + U.e2*x)/2; L.nx = (-1i*U.e2)/abs(U.e2) + 0*L.x; L.w=w*abs(U.e2);
R = L; R.x = L.x + U.e1;
B.x = (-U.e2 + U.e1*x)/2; B.nx = (1i*U.e1)/abs(U.e1) + 0*B.x; B.w=w*abs(U.e1);
T = B; T.x = B.x + U.e2;
% Proxy points for all periodic copies outside of the closest 3x3 boxes
Rp = 0.8; P=128; if nei==1, Rp = 1.27*abs(U.e1); P = 120; end
p.x = Rp * exp(1i*(1:P)'/P*2*pi); p = quadr(p); % proxy pts

% figure
% scatter(real(L.x),imag(L.x))
% hold on
% scatter(real(R.x),imag(R.x))
% scatter(real(T.x),imag(T.x))
% scatter(real(B.x),imag(B.x))
% scatter(real(p.x),imag(p.x))
% scatter(real(s.x),imag(s.x))
% axis equal
% pause(100)


% build matrices
% A = -eye(N)/2;   % exterior Neumann self matrix, starting w/ jump rel
% for i=-nei:nei
%     for j=-nei:nei  % direct sum
%     [~,Aij] = SLPmatrix(s,s,U.e1*i+U.e2*j); A = A + Aij;
%     end
% end
%s.len=N;
s=Quad(s);
%SM=Lap_SLP_Self_Matrix(s);
Aop=@(x) x/2+A_operator_Diri(s,x,U,nei); % A_operator = close(closeeval) & far(fmm) vesicle interactions.
% Aop=@(x)-x/2+eta*A_operator_Lap(s,x,U,nei); % A_operator = close(closeeval) & far(fmm) vesicle interactions.

% [~,Bm] = LapSLPmatrix(s,p);    % neu data from proxies.  name clash with B seg
[Bm, ~] = LapSLPmatrix(s,p);    % neu data from proxies.  name clash with B seg
%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%
% Bm = eta*Bm;            %%%
                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C=MakeC_Lap(s,L,R,T,B,U,nei);
C=MakeC_Diri(s,L,R,T,B,U,nei);

% w = [0*L.w,L.w,0*B.w,B.w]'; % col vec enforcing discrep consistency
% z = [zeros(M,1);ones(M,1);zeros(M,1);ones(M,1)]/2;
% w = w/(w'*z);
% C = C - z*(w'*C);   % rank-1 left-side correction to C

[QL,QLn] = LapSLPmatrix(L,p); [QR,QRn] = LapSLPmatrix(R,p);
[QB,QBn] = LapSLPmatrix(B,p); [QT,QTn] = LapSLPmatrix(T,p);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];

%E = [A Bm; C Q];   % expanded system matrix (don't really need in schur case)

% load weierstrass, load weierstrass_dip
% load jacobicn, load jacobicn_char
% load jacobisn, load jacobisn_char
% load jacobidn, load jacobidn_char

% load jacobisd, load jacobisd_char
% load jacobicd, load jacobicd_char
% load jacobind, load jacobind_char

% load jacobisc, load jacobisc_char
% load jacobinc, load jacobinc_char
% load jacobidc, load jacobidc_char

% load jacobins, load jacobins_char
load jacobids, load jacobids_char
% load jacobics, load jacobics_char

rhs = [0*s.x+charge(:); jumps(1)+0*L.x; 0*L.x; jumps(2)+0*B.x; 0*B.x]; % box driving
  Qdagg = Q\rhs(sum(s.len)+1:end); %norm(Qdagg) % for Schur rhs
  QdagC = Q\C;  % bkw stable
  %[UQ S V] = svd(Q,0); QdagC = V*(diag(min(1e12,1./diag(S)))*(UQ'*C)); % alt
  fprintf('norm QdagC = %.3g\n',norm(QdagC)) % should be O(1)
 % Ap = A - Bm*QdagC;  % Schur (in dense not fast apply form)
  %sig = Ap\(-Bm*Qdagg);  % direct solve, or...
  %sig = gmres(Ap, -Bm*Qdagg, [], 1e-12, N);  % iterative solve on fixed Ap
  %sig = gmres(@(x) A*x - Bm*(QdagC*x), -Bm*Qdagg, [], 1e-14, N); % "fast" form
  sig = gmres(@(x) Aop(x)-Bm*(QdagC*x), charge(:)-Bm*Qdagg,[],1e-14,50);
  %sig = gmres(@(x) A*x - Bm*(Q\(C*x)), -Bm*Qdagg, [], 1e-12, N);%stagnates why?
  %sig = gmres(@(x) A*x - Bm*(V*(diag(min(1e12,1./diag(S)))*(UQ'*(C*x)))), -Bm*Qdagg, [], 1e-12, N); %stagnates why?
  psi = Qdagg - QdagC*sig; co = [sig;psi]; % get remaining unknowns
%fprintf('resid norm = %.3g\n',norm(rhs - E*co))
fprintf('resid norm = %.3g\n',norm(rhs - [Aop(sig)+Bm*psi;C*sig+Q*psi]))
fprintf('sigma norm = %.3g, proxy strength norm = %.3g\n',norm(sig), norm(psi))
if v>1, figure; plot(E*co - rhs); title('residual vector'); end
fprintf('integral of sigma = %.3g\n',s.w'*sig)


% n=256;
% [X,Y]=meshgrid(-0.5:1/(n-1):0.5);
[X,Y]=meshgrid(-pi:h:pi);
n=length(X);

t.x=reshape(X,n^2,1)+1i*reshape(Y,n^2,1);
u=zeros(n^2,1);
ux=u;
uy=u;
ls=0;
while ls<n^2
    lsmax=min(ls+2e5,n^2);
    tt.x=t.x(ls+1:lsmax);
    [uu,uux,uuy]=At_operator_Diri(tt,s,sig,U,nei);
%     [uu,~,~]=At_operator_Diri(tt,s,sig,U,nei);
    u(ls+1:lsmax)=u(ls+1:lsmax)+LapSLPmatrix(tt,p)*psi+uu;
    tt.nx=0*tt.x+1;
    [~,UX]=LapSLPmatrix(tt,p);
    ux(ls+1:lsmax)=ux(ls+1:lsmax)+UX*psi+uux;
    tt.nx=0*tt.x+1i;
    [~,UY]=LapSLPmatrix(tt,p);
    uy(ls+1:lsmax)=uy(ls+1:lsmax)+UY*psi+uuy;
    ls=lsmax;
end

u=reshape(u,n,n);
ux=reshape(ux,n,n);
uy=reshape(uy,n,n);

% figure
% pcolor(X,Y,u)
% hold on
% colormap(jet(1000))
% axis equal
% axis([-0.5 0.5 -0.5 0.5]*2*pi)
% shading flat
% minu=min(min(u));
% maxu=max(max(u));
% LCurv=linspace(minu,maxu,100);
% LCurv=LCurv(2:end-1);
% contour(X,Y,u,LCurv,'k','LineWidth',1)
% ls=0;
% for k=1:s.M
%     for i=-1:1
%         for j=-1:1
%             fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+i*abs(U.e1),imag(s.x([ls+1:ls+s.len(k),ls+1]))+j*abs(U.e2),[0.8 0.8 0.8],'LineWidth',1)
%         end
%     end
%    ls=ls+s.len(k);
% end

% Plot error of the test function
poles = pi/2*[1+1i,1-1i,-1-1i,-1+1i];
for i = 1:4
    ind = (X-real(poles(i))).^2+(Y-imag(poles(i))).^2<=.5^2;
    u(ind) = NaN; jFun(ind)=NaN;
end
% fac = (u(2)-u(1))/(wFun(2)-wFun(1));
% u = u/fac;
figure(1), subplot(1,3,1), imagesc(X(1,:),Y(:,1),log10(abs(u-jFun))), colorbar, title('u error'), axis square off

% Plot error of the test function's derivatives
for i = 1:4
    ind = (X-real(poles(i))).^2+(Y-imag(poles(i))).^2<=.5^2;
    ux(ind) = NaN; jxFun(ind)=NaN;
    uy(ind) = NaN; jyFun(ind)=NaN;
end
figure(2), subplot(2,2,1), surf(X,Y,ux,'EdgeAlpha',0), title('ux'), axis square
        subplot(2,2,2), surf(X,Y,jxFun,'EdgeAlpha',0), title('jxFun'), axis square
        subplot(2,2,3), surf(X,Y,uy,'EdgeAlpha',0), title('uy'), axis square
        subplot(2,2,4), surf(X,Y,jyFun,'EdgeAlpha',0), title('jyFun'), axis square
figure(1), subplot(1,3,2), imagesc(log10(abs(ux-jxFun))), colorbar, title('ux error'), axis square off
        subplot(1,3,3), imagesc(log10(abs(uy-jyFun))), colorbar, title('uy error'), axis square off

keyboard

% figure
% pcolor(X,Y,sqrt(ux.^2+uy.^2))
% hold on
% colormap(jet(1000))
% shading flat
% axis equal
% axis([-0.5 0.5 -0.5 0.5]*2*pi)
% ls=0;
% for k=1:s.M
%     for i=-1:1
%         for j=-1:1
%             fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+i*abs(U.e1),imag(s.x([ls+1:ls+s.len(k),ls+1]))+j*abs(U.e2),[0.8 0.8 0.8],'LineWidth',1)
%         end
%     end
%     ls=ls+s.len(k);
% end

end