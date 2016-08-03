function Z=fb_Matrix(x,y,kb)
%
% Generates matrix for the bending force.
%

L=length(x)-1;
dsdt=Dsdt(x,y);
Z=eye(L);
Z=Z([1:end,1],:);
for k=1:4
Z=Dtp(Z);
Z=Z./(dsdt*ones(1,L));
end
Z=Z(1:end-1,:);
% We multiply by the bending modulus.
Z=-kb*[Z,0*Z;0*Z,Z];
end