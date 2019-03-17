%%
%Afunction to caculate the N in equation (52)
function lambdaK_general = lambdaK_generalfbUN(u,z)
global k0 k1 k2 k3 visc perm b bdim I oo Kg mug Kf porosity alpha Kdry x2 a;
global mudry taudry L Kstar M HBiot Mdim HBiotdim Ldim rhodrydim rhofdim citaw;
k0=oo./sqrt(Ldim);
c1=sqrt(HBiotdim);
k1=oo./c1;
k2=sqrt(I.*oo.*bdim.*HBiotdim./(Ldim.*Mdim));
k3=oo;
s1=k1./oo;
s2=k2./oo;
s3=k3./oo;
X1=-(HBiotdim.*s1.^2-1)./(alpha.*Mdim.*s1.^2-rhofdim);
X2=-HBiotdim./(alpha.*Mdim);
X3=I.*(rhofdim.*oo./bdim);
cigma1=(HBiotdim-alpha.*Mdim)./2;
cigma2=Ldim./(2.*alpha);
cigma3=(HBiotdim-alpha.*Mdim)./2;
cigma4=Mdim.*Ldim./(2.*HBiotdim.*rhofdim);
if u<=k1
q1=-I.*sqrt(k1.^2-u.^2);
else q1=sqrt(u.^2-k1.^2);
end
q2=-I.*sqrt(k2.^2-u.^2);

if u<=k3
q3=-I.*sqrt(k3.^2-u.^2);
else q3=sqrt(u.^2-k3.^2);
end
%%
%     UU=-1./4.*sin(2.*citaw).*(u+exp(I.*k1.*cos(citaw)).*(I.*k1.*cos(citaw).*sin(u)-u.*cos(u)))./(u.^2-k1.^2.*cos(citaw).^2);
    UU=-1./4.*sin(2.*citaw).*(-1+exp(I.*u.*k1.*cos(citaw)))/(I.*u+I.*k1.*cos(citaw));
% if u==-k1*cos(citaw);
%     UU=2;
% else
% UU=-1/8.*sin(2.*citaw).*2.*sinh(I.*(u+k1.*cos(citaw)))./(I.*(u+k1.*cos(citaw)));
% end
Mx=-4./pi.*(((u.^2-cigma2.*k2.^2)-2.*q3.*q2.*u.^2./(2.*u.^2-k3.^2)).*(1+X3)./(2.*(X3-X2).*(u.^2-cigma3.*k1.^2).*q2)+q3./(2.*u.^2-k3.^2)).*I.*u.*UU.*exp(I.*u.*z);
lambdaK_general=Mx;%*u/z;%the u/z is because of dimensionless
%/z
end