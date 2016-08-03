function [X,Y,SIG]=Reparameterize(X,Y,SIG,k)
%
% This function reparameterizes the vesicles so that discretization points
% are approximately equally spaced. 'X', 'Y', 'SIG', and 'S' are the same
% as in VWM.m. 'k' is the current time-step within VWM.m.
%
% See also:
%   Dt, INTt, It, VWM
%

% Variable index:
%   M = Number of vesicles.
%   n = The number of discretization points for the jth vesicle.
%   t = Vesicle parameterization variable. This is usually alpha in
%       literature.
%   ds = The derivative of arc length with respect to alpha for the jth
%       vesicle.
%   s = Arc length of the jth vesicle and also the desired arc length
%       values for the new discretization.
%   st = Integral of ds with respect to alpha.
%   ti = Alpha values for new discretization points.
%

M=length(Y{1});
for j=1:M
    n=length(X{k}{j});
    ds=sqrt(Dtp(X{k}{j}).^2+Dtp(Y{k}{j}).^2);
    s=sum(ds(1:n-1))*2*pi/(n-1);
    s=linspace(0,s,n)';
    st=INTtp(ds);
    
    ti=zeros(n,1);
    ti(1)=0;
    ti(n)=2*pi;
    jj=1;
    for ii=2:n-1
        while s(ii)<st(jj)||s(ii)>=st(jj+1)
            jj=jj+1;
        end
        ti(ii)=2*pi/(n-1)*(jj-1+(s(ii)-st(jj))/(st(jj+1)-st(jj)));
    end
    
    % Interpolates 'X', 'Y', and 'SIG' at these values
    xx=Itp([X{k}{j},Y{k}{j}],ti);
    X{k}{j}=xx(:,1);
    Y{k}{j}=xx(:,2);
    if k>1
        SIG{k-1}{j}=Itp(SIG{k-1}{j},ti);
    end
end
end