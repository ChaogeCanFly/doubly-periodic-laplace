function u=Lap_DLP_Self_Apply(s,f,S)
% DLP normal derivative self-evaluation correction for Lap_DLP_Far_Apply.m
% Bowei Wu, 5/25/16
M=length(s.len);
u.sn=zeros(s.cumm(end),1);
ls=0;
df = reshape(DmFT(reshape(f,[],s.M)),[],1)./s.sp;
for k=1:M
    ss.x=s.x(ls+1:ls+s.len(k));
    ss.sp=s.sp(ls+1:ls+s.len(k));
    ss.w=s.w(ls+1:ls+s.len(k));
    ss.nx=s.nx(ls+1:ls+s.len(k));
    ss.cur=s.cur(ls+1:ls+s.len(k));
    ss.t=s.t(ls+1:ls+s.len(k));
    ff=df(ls+1:ls+s.len(k));
    if nargin==2
        zz=LapSLPselfmatrix(ss)*ff;
        zz=D1FT(zz)./ss.sp;
    else
        zz=S{k}*ff;
        zz=D1FT(zz)./ss.sp;
    end
    [~,DLPn]=LapDLPmatrix(ss,ss,0);
    DLPn(logical(eye(s.len(k))))=0;
    DLPn=DLPn*f(ls+1:ls+s.len(k));
    u.sn(ls+1:ls+s.len(k))=zz-DLPn;
    ls=ls+s.len(k);
end

function df=D1FT(f)

persistent K
%
% This program computes the derivative of a function of the form 
% f(t)=a*t+p(t), where a is a constant and p(t) is a periodic function on 
% the interval [0,2*pi]. Both endpoints must be included.
%
% Example:
%   t=linspace(0,2*pi,33)';
%   y=sin(t);
%   dy=Dt(y);
%

[N,p]=size(f);

if isempty(K)||size(K,1)~=N||size(K,2)~=p
    if (-1)^N==1
        K=1i*[(0:N/2-1)';0;(-N/2+1:-1)']*ones(1,p);
    else
        K=1i*[(1:(N-1)/2)';(-(N-1)/2:0)']*ones(1,p);
    end
end

df=ifft(K.*fft(f));

