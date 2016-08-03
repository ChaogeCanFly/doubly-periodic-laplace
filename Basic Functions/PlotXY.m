function PlotXY(X,Y,n,a)

if nargin==3
    a='r-';
end
for k=1:length(X{1})
plot(X{n}{k},Y{n}{k},a)
hold on
plot(X{n}{k}+2*pi,Y{n}{k},a)
plot(X{n}{k}-2*pi,Y{n}{k},a)
end
[U,D]=MakeUD(4000);
plot(real(U.x),imag(U.x),'k','LineWidth',1)
plot(real(D.x),imag(D.x),'k','LineWidth',1)
axis equal
H=max(imag(U.x))-min(imag(D.x));
axis([0 2*pi -1.1*H/2 1.1*H/2])
xlabel('x')
ylabel('y')