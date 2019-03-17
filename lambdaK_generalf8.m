%%
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% this is a function usd in scatteringproblem, it can calculate the Kb in
% Song (2017) in equation (A13)
% v.0.1, 18/05/2018, Boye Fu & Boris Gurevich, Curtin University and CRGC
%%
function lambdaK_general = lambdaK_generalf8(u,z)
global k0 k1 k2 k3 visc perm b bdim I oo Kg mug Kf porosity alpha Kdry;
global mudry taudry L Kstar M HBiot Mdim HBiotdim Ldim rhodrydim rhofdim;

if z==u
K=1/2*(besselj(1,u)^2+besselj(0,u)*(2*besselj(1,u)/u-besselj(2,u)));% The function about Besselj function in (A14) in Song (2017)
else
K=(u*besselj(-1,u)*besselj(0,z)-z*besselj(-1,z)*besselj(0,u))/(z^2-u^2);

end
%%
%Now we will caculate the H in (A17) in Song (2017)
k0=oo/sqrt(Ldim);
c1=sqrt(HBiotdim);
k1=oo/c1;
k2=sqrt(I*oo*bdim*HBiotdim/(Ldim*Mdim));
k3=oo;
s1=k1/oo;
s2=k2/oo;
s3=k3/oo;
X1=0;
X2=-HBiotdim/(alpha*Mdim);
X3=I*(rhofdim*oo/bdim);
cigma1=(HBiotdim-alpha*Mdim)/2;
cigma2=Ldim/(2*alpha);
cigma3=(HBiotdim-alpha*Mdim)/2;
cigma4=Mdim*Ldim/(2*HBiotdim*rhofdim);
E=k3^2*(1*(1-2*cigma2)/cigma4-(1+X2)*(1-2*cigma1)/cigma3);
if u<=k1
q1=-I*sqrt(k1^2-u^2);
else q1=sqrt(u^2-k1^2);
end
q2=-I*sqrt(k2^2-u^2);

if u<=k3
q3=-I*sqrt(k3^2-u^2);
else q3=sqrt(u^2-k3^2);
end
%%
H1=4/(E*u);
H2=(u^2-cigma1*k1^2)/q1-2*u^2*q3/(2*u^2-k3^2);
H3=(u^2-cigma2*k2^2)/q2-2*u^2*q3/(2*u^2-k3^2);
H4=(X3)/(X3-X2)*(u^2-cigma4*k2^2)/(u^2-cigma3*k1^2);
FH=H2-H3*H4;
H5=u;
H=H5*(H1*(X3-X2)*FH*(u^2-cigma3*k1^2)-1);
%%
lambdaK_general=K*H;%*u/z;%the u/z is because of dimensionless
%/z
end