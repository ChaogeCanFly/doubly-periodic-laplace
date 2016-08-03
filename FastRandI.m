function [I,N]=FastRandI(clup,M,n)

% n=128;
if nargin==0
M=4;
clup=10;
end

addpath(genpath('.'))
S=[];
ls=1;
while size(S,1)/n<M
    R=0.5*1/1.1^ls/4;
    MR=round(5/R^2);
    S0=MakeSofR(R,MR,n);
    %S0=MakeBarnettShape(R,MR,n);
    S0=IslandIntersect(S0,S,n);
    S=IslandIntersect(S,S0,n);
    S0=RecursiveIntersect(S0,n);
    S0=RecursiveIntersect(S0,n);
    S=[S;S0];
    ls=ls+1;
end
S=S(1:M*n,:);

% size(S,1)/n


% ls=0;
% for k=1:M
%    ss=S([ls+1:ls+n,ls+1],:);
%    area=Area(ss(:,1),ss(:,2));
%    xc=Area(ss(:,1).^2/2,ss(:,2))/area;
%    yc=Area(ss(:,1).*ss(:,2),ss(:,2))/area;
%    S(ls+1:ls+n,:)=[CScale*(S(ls+1:ls+n,1)-xc)+xc,CScale*(S(ls+1:ls+n,2)-yc)+yc];
%    ls=ls+n;
% end
% Error=ClosenessFactor(S,n,M)-clup

maxp=zeros(M,1);
ls=0;
for k=1:M
   maxp(k)=ArcLength(S([ls+1:ls+n,ls+1],1),S([ls+1:ls+n,ls+1],2),2*pi);
   ls = ls+n;
end
maxp=max(maxp);
maxd=maxp/clup;

CScale=FindCScale(S,n,M,clup,maxd);
ClupErrorFun=@(x)ClupError(x,S,n,M,clup,maxd);
CScale=fzero(ClupErrorFun,[0 CScale]);

ls=0;
%rs=0;
for k=1:M
    I{k}.x=S(ls+1:ls+n,:)*[1;1i]-pi-1i*pi;
    area=Area(real(I{k}.x([1:end,1])),imag(I{k}.x([1:end,1])));
    xc=Area(real(I{k}.x([1:end,1])).^2/2,imag(I{k}.x([1:end,1])))/area;
    yc=Area(real(I{k}.x([1:end,1])).*imag(I{k}.x([1:end,1])),imag(I{k}.x([1:end,1])))/area;
    I{k}.x=CScale*(I{k}.x-xc-1i*yc)+xc+1i*yc;
    ls=ls+n;
end
N=n*ones(length(I),1);


for k=1:M
    for j=-1:1
        for i=-1:1
            plot(real(I{k}.x([1:end,1]))+2*pi*i,imag(I{k}.x([1:end,1]))+2*pi*j)
            hold on
        end
    end
end
axis equal
axis([-pi pi -pi pi])

end

function S=MakeSofR(R,M,n)
t=linspace(0,2*pi,n+1)';
t=t(1:end-1);
R=R*ones(M,1);
x0=2*pi*rand(M,1);
y0=2*pi*rand(M,1);
Sx=reshape(cos(t)*R'+ones(n,1)*x0',n*M,1);
Sy=reshape(sin(t)*R'+ones(n,1)*y0',n*M,1);
S=[Sx,Sy];
end

function S=MakeBarnettShape(R,M,n)
t=linspace(0,2*pi,n+1)';
t=t(1:end-1);
t=t*ones(1,M);
alpha=0.5*ones(n,1)*rand(1,M);
beta=2*pi*ones(n,1)*rand(1,M);
w=floor(7*ones(n,1)*rand(1,M)+2);
R=R*(1+alpha.*cos(w.*t+beta));
x=cos(t).*R;
y=sin(t).*R;
x0=2*pi*ones(n,1)*rand(1,M);
y0=2*pi*ones(n,1)*rand(1,M);
Sx=reshape(x+x0,n*M,1);
Sy=reshape(y+y0,n*M,1);
S=[Sx,Sy];
end


function clup=ClosenessFactor(S,n,M,maxd)
SP=zeros(n*M*9,2);
ls=0;
for i=-1:1
    for j=-1:1
        SP(ls+1:ls+n*M,:)=[S(:,1)+i*2*pi,S(:,2)+j*2*pi];
        ls=ls+n*M;
    end
end
Mdl = KDTreeSearcher(S);
[cn,cd] = knnsearch(Mdl,SP,'k',n+1);
% maxd=10;
% [idx,D]= rangesearch(Mdl,SP,maxd);
% lenidx=0;
% for k=1:length(idx)
%    lenidx=lenidx+length(idx{k});
% end
% lenidx
cl=zeros(9*M,1);
ls=0;
for j=1:9
    rs=0;
    for k=1:M
          cnk=cn(ls+1:ls+n,:);
         cdk=cd(ls+1:ls+n,:);
%         cnk=idx(ls+1:ls+n)
%         cdk=D(ls+1:ls+n)
        
%         dmin=[];
%         ns=1;
%         for i=1:n
%             if j==5
%                 cdkin=cnk{i}<rs+1|cnk{i}>rs+n;
%             else
%                 cdkin=true(size(cdk{i}));
%             end
%             if ~isequal(sum(cdkin),0)
%                 dmin(ns)=min(cdk{i}(cdkin));
%                 ns=ns+1;
%             end
%         end
%         
%         if isempty(dmin)
%             cl(M*(j-1)+k)=0;
%         else
%             cl(M*(j-1)+k)=ArcLength(S([rs+1:rs+n,rs+1],1),S([rs+1:rs+n,rs+1],2),2*pi)/min(dmin);
%         end
        
        
        if j==5
            cdkin=cnk<rs+1|cnk>rs+n;
        else
            cdkin=true(size(cdk));
        end
        
        if ~isempty(cdkin)
        [cdmin,cdloc]=min(cdk(cdkin));
        cnin=cnk(cdkin);
        cdv=ceil(cnin(cdloc)/n);
        ss1=SP(ls+1:ls+n,:);
        ss2=S(n*(cdv-1)+1:n*cdv,:);
        x0=2*pi*[(cdloc-1)/n-ceil(cdloc/n)+1;(cnin(cdloc)-1)/n-(cdv-1)];
        [xmin,dmin]=fminsearch(@(x)abs((Itp(ss1([1:end,1],:),mod(x(1),2*pi))-Itp(ss2([1:end,1],:),mod(x(2),2*pi)))*[1;1i]),x0);
        end
        
        if isempty(cdkin)
            cl(M*(j-1)+k)=0;
        else
            %cl(M*(j-1)+k)=ArcLength(S([rs+1:rs+n,rs+1],1),S([rs+1:rs+n,rs+1],2),2*pi)/min(min(cdk(cdkin)));
            cl(M*(j-1)+k)=ArcLength(S([rs+1:rs+n,rs+1],1),S([rs+1:rs+n,rs+1],2),2*pi)/dmin;
        end
        ls=ls+n;
        rs=rs+n;
    end
end
clup=max(cl);
end

function Error=ClupError(CScale,S,n,M,clup,maxd)
ls=0;
for k=1:M
    ss=S([ls+1:ls+n,ls+1],:);
    area=Area(ss(:,1),ss(:,2));
    xc=Area(ss(:,1).^2/2,ss(:,2))/area;
    yc=Area(ss(:,1).*ss(:,2),ss(:,2))/area;
    S(ls+1:ls+n,:)=[CScale*(S(ls+1:ls+n,1)-xc)+xc,CScale*(S(ls+1:ls+n,2)-yc)+yc];
    ls=ls+n;
end
Error=ClosenessFactor(S,n,M,2*CScale*maxd)-clup;
end

function CScale=FindCScale(S,n,M,clup,maxd)
CS=[1 1];
S0=0*S;
Error=-1;



while Error<=0
    CScale=CS(2);
    ls=0;
    hold off
    for k=1:M
        ss=S([ls+1:ls+n,ls+1],:);
        area=Area(ss(:,1),ss(:,2));
        xc=Area(ss(:,1).^2/2,ss(:,2))/area;
        yc=Area(ss(:,1).*ss(:,2),ss(:,2))/area;
        S0(ls+1:ls+n,:)=[CScale*(S(ls+1:ls+n,1)-xc)+xc,CScale*(S(ls+1:ls+n,2)-yc)+yc];
%         plot(S0(ls+1:ls+n,1),S0(ls+1:ls+n,2))
%         hold on
        ls=ls+n;
    end
%     axis equal
%     axis([0 2*pi 0 2*pi])
%     pause(0.1)
    S1=RecursiveIntersect(S0,n);
    S1=RecursiveIntersect(S1,n);
    if isequal(S0,S1)
        CS=[CS(2) 2*CS(2)];
        Error=ClupError(CScale,S,n,M,clup,maxd);
    else
        CS=[CS(1) 0.5*(CS(1)+CS(2))];
    end
end
end