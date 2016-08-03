function Z=P_Matrix(x,y)
%
% Generates matrix for the inextensibility operator.
%

% Initializes variables
L=length(x)-1;

% Computes the derivative of arc length with respect to the
% parameterization variable.
dsdt=Dsdt(x,y);

% Applies the D operator on each vesicle.
EYE=eye(L);
EYE=EYE([1:end,1],:);
[dxds,dyds]=Ds(x,y,dsdt);
Sfs=Dtp(EYE)./(dsdt*ones(1,L));
Z=[(dxds*ones(1,L)).*Sfs,(dyds*ones(1,L)).*Sfs];
Z=Z(1:end-1,:);