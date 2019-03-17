%%
%A function to calculate the S in equation (52)
function lambdaK_general = lambdaK_generalfbUS(u,z)
global k0 k1 k2 k3 visc perm b bdim I oo Kg mug Kf porosity alpha Kdry x2 a;
global mudry taudry L Kstar M HBiot Mdim HBiotdim Ldim rhodrydim rhofdim citaw;
k0=oo/sqrt(Ldim);
c1=sqrt(HBiotdim);
k1=oo/c1;
k2=sqrt(I.*oo.*bdim.*HBiotdim./(Ldim.*Mdim));
k3=oo;
s1=k1./oo;
s2=k2./oo;
s3=k3./oo;
X1=0;
X2=-HBiotdim./(alpha.*Mdim);
X3=I.*(rhofdim.*oo./bdim);
cigma1=(HBiotdim-alpha.*Mdim)./2;
cigma2=Ldim./(2.*alpha);
cigma3=(HBiotdim-alpha.*Mdim)./2;
cigma4=Mdim.*Ldim./(2.*HBiotdim.*rhofdim);
E=(2*HBiotdim*alpha*Mdim*(alpha-Ldim)/(Mdim*Ldim*alpha*HBiotdim)+2/(alpha*Mdim)*(1-HBiotdim+alpha*Mdim))*k3^2;%The E value in equation (56) in Song (2017), which is related to calculation and velocity
q1=sqrt(u.^2-k1.^2);
q2=-I.*sqrt(k2.^2-u.^2);
q3=sqrt(u.^2-k3.^2);
%%
H1=4./(E.*u);
H2=(X3-X2).*(u.^2-cigma1.*k1.^2).*(u.^2-cigma3.*k1.^2)./q1;
H3=(X1-X3).*(u.^2-cigma2.*k2.^2).*(u.^2-cigma4.*k2.^2)./q2;
H4=(X1-X2).*u.^2.*q3;
FH=H2+H3-H4;
H=H1.*FH-1;
%%
%     UU=-1/4.*sin(2.*citaw).*(u+exp(I.*k1.*cos(citaw)).*(I.*k1.*cos(citaw).*sin(u)-u.*cos(u)))./(u.^2-k1.^2.*cos(citaw).^2);
% UU=-1./4.*sin(2.*citaw).*(-1+exp(I.*u.*k1.*cos(citaw)));
UU=-1./4.*sin(2.*citaw).*(-1+exp(I.*u.*k1.*cos(citaw)))/(I.*u+I.*k1.*cos(citaw));
% if u==-k1*cos(citaw);
%     UU=2;
% else
% UU=-1/8.*sin(2.*citaw).*2.*sinh(I.*(u+k1.*cos(citaw)))./(I.*(u+k1.*cos(citaw)));
% end
W=-I.*E./pi./(X2-X1)./k3.^2.*u.*UU./(u.^2-cigma3.*k1.^2);
HKK11=u.*(H+1).*W.*exp(I.*u.*z);
%%
lambdaK_general=HKK11;%*u/z;%the u/z is because of dimensionless
end