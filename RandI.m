function [I,N]=RandI

n=64;
M=1;

t=linspace(0,2*pi,n+1)';
t=t(1:end-1);
%r=n*real(ifft(rand(n,1).*8.^(-(1:n)')));
r=ones(n,1);
r=rand*(r-min(r))./max(r)+0.4*(rand/2+0.5);
x=r.*cos(t);
y=r.*sin(t);

x0=pi*rand-pi/2;
y0=pi*rand-pi/2;

X=[x0+x+1i*(y0+y);nan+1i*nan];
% plot(x0+x,y0+y)
% hold on
% axis equal
% axis([-pi pi -pi pi])
% pause(0.1)

ls=1;
while ls<M
    %r=n*real(ifft(rand(n,1).*8.^(-(1:n)')));
    r=ones(n,1);
    r=rand*(r-min(r))./max(r)+0.4*(rand/2+0.5);
    x=r.*cos(t);
    y=r.*sin(t);
    x0=2*pi*rand-pi;
    y0=2*pi*rand-pi;
    in1=inpolygon(x0+x,y0+y,real(X(1:end-1)),imag(X(1:end-1)));
    in2=inpolygon(real(X(~isnan(X))),imag(X(~isnan(X))),x0+x,y0+y);
    if max(in1)+max(in2)==0&&min(x0+x)>-pi&&max(x0+x)<pi&&min(y0+y)>-pi&&max(y0+y)<pi
        X=[X;x0+x+1i*(y0+y);nan+1i*nan];
%         plot(x0+x,y0+y)
%         hold on
%         axis equal
%         axis([-pi pi -pi pi])
%         pause(0.01)
        ls=ls+1;
    end
end

I=cell(M,1);
N=n*ones(M,1);
ls=0;
for k=1:M
   I{k}.x=X(ls+1:ls+n);
   ls=ls+n+1;
end
        