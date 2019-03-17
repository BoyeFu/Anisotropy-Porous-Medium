%%
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% The code is got according to Galvin's thesis Galvin, R. (2007), Elastic wave attenuation, dispersion and anisotropy in fractured porous media. Doctoral dissertation, Curtin University.
% This is a extension of the code of Galvin and Gurevich, we extent the
% normal incidence to oblique incidence.
%
%   The function contain some self-defining function
%   lambdaK_generalf8: Calculate the KB in equation A13 in Song (2017)
%   gff1:calculate the right hand of A13 in Song (2007)
%   lambdaK_generalf:Calculate the R(x,y)T(y) in equation (29) in Galvin
%   (2009):Galvin, R. J., & Gurevich, B. (2009). Effective properties of a poroelastic medium containing a distribution of aligned cracks. Journal of Geophysical Research: Solid Earth, 114(B7).
%   gff: calculate the -p0S(x) in equation (29)in Galvin (2009)
%   Unit of The elastic moduli in the code is Pa, and we do  normalization by the shear modulus of dry rock
%   the unit of density in the code are kg/m3
%%
clear all;
t0=cputime;
global kk_e h tau L_e Ldim_e mu_e k0 k1 k2 k3 kk k2_g visc perm b bdim I citaw;
global rrr oo Kg mug Kf porosity alpha Kdry mudry taudry L Kstar M HBiot Kfdim Kgdim;
global Mdim HBiotdim Ldim rhofdim rhodrydim A1_g A2_g q1 q2 q3;%the modulus for dimensionless 
%%
%Some basic parameters we used in this paper
I=sqrt(-1);% the  imaginary unit 
n0=0.01;% the slit crack density
Kg=37*10^9;%Bulk modulus of rock matrix
mug=44*10^9;%shear modulus of rock matrix
Kf=2.25*10^9;%bulk modulus of pore fluid

porosity=0.1;%porosity
visc=1e-3;%viscous cofficient
perm=1e-13;%permeability
rhof=1e3;%the density of pore fluid
rhos=2.65e3;%density of rock matrix
rhodry=(1-porosity)*rhos;%density of dry rock
rho=(1-porosity)*rhos+porosity*rhof;%density used in Biot equation
b=visc/perm;%the ratio between viscous cofficient and permeability
alpha=1-(1-porosity)^(3/(1-porosity));%the Biot coefficient
Kdry=Kg*(1-alpha);%the dry rock bulk modulus
mudry=mug*(1-alpha);%the dry rock shear modulus
bdim=b/sqrt(rho*mudry);%%I a middle varible-the b for dimensionless
L=Kdry+4*mudry/3;%The P-wave modulus for dry rocks;
taudry=sqrt(mudry/L);%?? a middle variable dimensionless The "g" in Galvin (2009)
M=Kg/((1-Kdry/Kg)-porosity*(1-Kg/Kf));% the Biot modulus of the rock The M in Gassmann equation
HBiot=L+alpha^2*M;%Fast P-wave modulus
Mdim=M/mudry;%The M value for dimensionless
HBiotdim=HBiot/mudry;% the value for dimensionless-Biot Fast p wave modulus
Ldim=L/mudry;% the dimensionless P-wave modulus
Kfdim=Kf/mudry;
rhodrydim=rhodry/rho;%the dimensionless dry rock density
rhofdim=rhof/rho;%the dimensionless dry rock density of rock
Kgdim=Kg/mudry;%the dimensionless bulk modulus of rock matrix
%%
cita=[pi./2 pi./4 pi];
for qqq=1:length(cita);
citaw=cita(qqq);
%%


omegadim=[10^-4.5 0.0001 10^-3.5 0.001 10^-2.5 0.01 10^-1.5 0.1 10^-0.5 10^-0.25 1 10];% frequency
%%
% Now we will calculate the velocity and attenuation in different
% frequency, the method have used by Kawahara and Yamashita (1992) in
% equation (17):Kawahara, J., & Yamashita, T. (1992). Scattering of elastic waves by a fracture zone containing randomly distributed cracks. pure and applied geophysics, 139(1), 121-144.
%%
% oo=omegadim(3);% the frequency
% c0=sqrt(Ldim/rhodrydim);%the dimensionless p-wave velcocity for dry rock
% k0=oo/sqrt(Ldim);%the dimensionless P-wave number
% c1=sqrt(HBiotdim);%the dimensionless Fast P-wave velocity
% k1=oo/c1;%The dimensionless Fast P-wave number
% k2=sqrt(I*oo*bdim*HBiotdim/(Ldim*Mdim));%The dimensionless Slow P-wave number
% k3=oo;%the demensionless S-wave number(the demensionless is relation to shear modulus, therefore, the demensionless shear modulus is 1)
% kk=abs(k2);%the abslute value of slow P-wave number
% h=k1;
% x1=0.1;% the low limit of the integration
% if kk <= 1
% spac=2*kk;
% x2=4;
% x2=200*spac;
% spac=0.1*kk;
% else
% spac=0.002*kk;
% x2=40;
% end
% x2=max(1000,kk*50);%the high limit of the integration, becuase the integration region is (0,infinity) in equation (A13),we use the value x2 represent infinity
% 
% NMAX=500;
% % omk=eye(NMAX);%the number of discreted net, the detail is in Kawahara and Yamashita (1992) in
% % equation (17)
% % CALCULATE GAUSS-LEGENDRE WEIGHTS, In this paper, we use GAUSS-LEGENDRE
% % integration to deal with the integration in equation (A13)
% m=(NMAX+1)/2;%%£¿
% xm=(x2+x1)/2;%midpoint 
% xl=(x2-x1)/2;%midpoint
% for i=1:m
% Z=cos(pi*(i-1/4)/(NMAX+1/2));%1>Z>0
% Z1=Z-1;%0>Z1>-1
% ZZ(floor(i))=Z;
% ZZZ(i)=Z-Z1;
% while (Z-Z1) > 3d-14
% p1=1;
% p2=0;
% for j=1:NMAX
% p3=p2;
% p2=p1;%=((2*(j-1)-1)*Z*p2-(j-2)*p3)/j
% p1=((2*j-1)*Z*p2-(j-1)*p3)/j;%p1=((2*j-1)*Z*((2*(j-1)-1)*Z*p2-(j-2)*p2)/j-(j-1)*p2)/j
% end
% pp=NMAX*(Z*p1-p2)/(Z^2-1);%z=1 pp=inf z=0, NMAX((j-1)(p2-p1)/j)
% Z1=Z;
% Z=Z1-p1/pp;
% end
% x(i)=xm-xl*Z;%x1<x<(x2+x1)/2
% x(NMAX+1-i)=xm+xl*Z;%x2>x>(x2+x1)/2
% end
% for j=1:NMAX
% z(j)=x(j);
% end
% for j=1:NMAX
% %     gmat11(j)=real(gff1(z(j)));
%     gmat11(j)=(gff111(z(j)));
% end

%%
for q=1:length(omegadim)
    %The basic parameters for calculation 
oo=omegadim(q);% the frequency
c0=sqrt(Ldim/rhodrydim);%the dimensionless p-wave velcocity for dry rock
k0=oo/sqrt(Ldim);%the dimensionless P-wave number
c1=sqrt(HBiotdim);%the dimensionless Fast P-wave velocity
k1=oo/c1;%The dimensionless Fast P-wave number
k2=sqrt(I*oo*bdim*HBiotdim/(Ldim*Mdim));%The dimensionless Slow P-wave number
k3=oo;%the demensionless S-wave nuber(the demensionless is relation to shear modulus, therefore, the demensionless shear modulus is 1)
kk=abs(k2);%the abslute value of slow P-wave number
kk2(q)=kk;
omeg(q)=kk^2;% the the abslute value of slow P-wave number sqrt
h=k1;
x1=0;% the low limit of the integration
if kk <= 1
spac=2*kk;
x2=4;
x2=200*spac;
spac=0.1*kk;
else
spac=0.002*kk;
x2=40;
end
x2=max(1000,kk*50);%the high limit of the integration, becuase the integration region is (0,infinity) in equation (A13),we use the value x2 represent infinity
p=-(HBiotdim-alpha*Mdim-2.*(1-sin(citaw).^2));%the f in equation A13
NMAX=500;
omk=eye(NMAX);%the number of discreted net, the detail is in Kawahara and Yamashita (1992) in
% equation (17)
% CALCULATE GAUSS-LEGENDRE WEIGHTS, In this paper, we use GAUSS-LEGENDRE
% integration to deal with the integration in equation (A13)
m=(NMAX+1)/2;%%£¿
xm=(x2+x1)/2;%midpoint 
xl=(x2-x1)/2;%midpoint
NN(q)=xm;
MM(q)=xl;
for i=1:m
Z=cos(pi*(i-1/4)/(NMAX+1/2));%1>Z>0
Z1=Z-1;%0>Z1>-1
ZZ(floor(i))=Z;
ZZZ(i)=Z-Z1;
while (Z-Z1) > 3d-14
p1=1;
p2=0;
for j=1:NMAX
p3=p2;
p2=p1;%=((2*(j-1)-1)*Z*p2-(j-2)*p3)/j
p1=((2*j-1)*Z*p2-(j-1)*p3)/j;%p1=((2*j-1)*Z*((2*(j-1)-1)*Z*p2-(j-2)*p2)/j-(j-1)*p2)/j
end
pp=NMAX*(Z*p1-p2)/(Z^2-1);%z=1 pp=inf z=0, NMAX((j-1)(p2-p1)/j)
Z1=Z;
Z=Z1-p1/pp;
end
x(i)=xm-xl*Z;%x1<x<(x2+x1)/2
x(NMAX+1-i)=xm+xl*Z;%x2>x>(x2+x1)/2
w(i)=2*xl/((1-Z^2)*pp^2);%(x2-x1)/2
w(NMAX+1-i)=w(i);
end
for i=1:NMAX
u(i)=x(i);
end
for j=1:NMAX
z(j)=x(j);
end
for i=1:NMAX
for j=1:NMAX
omkk_g(i,j)=omk(i,j)+lambdaK_generalf8(u(i),z(j))*w(i);%the coefficient matrix of equation (A13), a discretization of the integration of KB in equation (A13) in Song (2017)
end
end
for j=1:NMAX
% gmat1(q,j)=gmat11(j);
gmat1(q,j)=gff111(j);
    gmat2(q,j)=gff222(z(j));
gmatg(q,j)= gmat1(q,j)+ gmat2(q,j);%The discretization of the right hand in equation (A13) in Song (2017)
end
f_g=linsolve(omkk_g.',gmatg(q,:).');
B0=f_g(1);% The B value in equation (A13) in Song (2017)
%%
%The slowness of elastic wave
s1(q)=k1/oo;
s2(q)=k2/oo;
s3(q)=k3/oo;
%%
%X1 X2 X3 is the low frequency limit of the X1 X2 X3 equation (24)-(25)
X1(q)=0;
X2(q)=-HBiotdim/(alpha*Mdim);
X3(q)=I*(rhofdim*oo/bdim);
%%
%cigma1 cigma2 cigma3 cigma4 is the low frequency limit of the cigma1
%cigma2 cigma3 cigma4 in equation (30)-(33)
cigma1(q)=(HBiotdim-alpha*Mdim)/2;
cigma2(q)=Ldim/(2*alpha);
cigma3(q)=(HBiotdim-alpha*Mdim)/2;
cigma4(q)=Mdim*Ldim/(2*HBiotdim*rhofdim);
E(q)=(2*HBiotdim*alpha*Mdim*(alpha-Ldim)/(Mdim*Ldim*alpha*HBiotdim)+2/(alpha*Mdim)*(1-HBiotdim+alpha*Mdim))*k3^2;%The E value in equation (49)
A0=(3.14*(-X2(q))*(0-cigma3(q)*k1^2)*B0/E(q)/(sqrt(-1)));%The amplitude of P wave according to equation (45)
Breco(q)=B0;
ff(q)=((HBiotdim-alpha*Mdim-2.*(1-sin(citaw).^2)))/((HBiotdim-alpha*Mdim))*(I^3)/pi*A0;%;g according to equation (60)
Q1(qqq,q)=4*pi*n0*imag(ff(q));%The attenuation in the whole frequency region
%%
%The parameters used in Appendix C
citaww=pi/2-citaw;
Csath11=HBiot;
Csath22=HBiot;
Csath12=HBiot-2.*mudry;
Csath66=(mudry).*(1-Ldim.*1.*pi.*n0./(2.*(HBiotdim-alpha.^2.*Mdim-1).*(1)));
MMouh=((Csath11-Csath66).*sin(citaww)^2-(Csath22-Csath66).*cos(citaww).^2).^2+(Csath12+Csath66).^2.*sin(2.*citaww).^2;
Csath0(qqq)=(Csath11.*sin(citaww).^2+Csath22.*cos(citaww).^2+Csath66+sqrt(MMouh))./2;
CCsats(qqq,q)=(Csath0(qqq)/rho.*(1-4*pi*n0*real(ff(q))));
CCsatsc(qqq,q)=(Csath0(qqq)/rho.*(1-4*pi*n0*ff(q)));
v(qqq,q)=sqrt(Csath0(qqq)/rho.*(1-4*pi*n0*real(ff(q))));%The velocity in the whole frequency region
v2(qqq,q)=(real(sqrt((Csath0(qqq)/rho.*(1-4*pi*n0*(ff(q)))))).^(-1)).^(-1);%The velocity in the whole frequency region
v3(qqq,q)=((sqrt((Csath0(qqq)/rho.*(1-4*pi*n0*(ff(q)))))).^(-1)).^(-1);%The velocity in the whole frequency region
Q2(qqq,q)=imag(CCsatsc(qqq,q))./real(CCsatsc(qqq,q));
end

%%
%The parameters used in interpolation method
T=2*(HBiotdim-alpha*Mdim)^2*(2-4*alpha*taudry^2+3*alpha^2*taudry^4)*n0*bdim/(15*taudry^2*(1-taudry^2)^2*HBiotdim*Ldim);
G=2*pi*n0*(HBiotdim-alpha*Mdim)^2*sqrt((1/(HBiotdim*Mdim*Ldim*bdim)));
Csatl11=HBiot.*(1-(HBiotdim-alpha.^2.*Mdim).*(HBiotdim-alpha.*Mdim-2).^2.*pi.*n0./(2.*(HBiotdim-alpha^2.*Mdim-1).*HBiotdim));
Csatl22=HBiot.*(1-(HBiotdim-alpha.*Mdim).^2.*(HBiotdim-alpha.^2.*Mdim).*pi.*n0./(2.*HBiotdim.*(HBiotdim-alpha.^2.*Mdim-1)));
Csatl12=(HBiot-2.*mudry).*(1-(HBiotdim.^2-alpha.*(1+alpha).*HBiotdim.*Mdim+alpha.^3.*Mdim.^2).*(HBiotdim-alpha.*Mdim-2).*pi.*n0./(2.*(HBiotdim-alpha.^2.*Mdim-1).*(HBiotdim-2)));
Csatl66=(mudry).*(1-Ldim.*1.*pi.*n0./(2.*(HBiotdim-alpha.^2.*Mdim-1).*(1)));
tao11=((Csath22-Csatl22)/(Csath22*G))^2;
kesisi11=(Csath22-Csatl22)^3/(2*Csatl22*Csath22^2*T*G^2);
tao22=((Csath22-Csatl22)/(Csath22*G))^2;
kesisi22=(Csath22-Csatl22)^3/(2*Csatl22*Csath22^2*T*G^2);
tao12=((Csath22-Csatl22)/(Csath22*G))^2;
kesisi12=(Csath22-Csatl22)^3/(2*Csatl22*Csath22^2*T*G^2);
tao66=((Csath22-Csatl22)/(Csath22*G))^2;
kesisi66=(Csath22-Csatl22)^3/(2*Csatl22*Csath22^2*T*G^2);
for q=1:length(omegadim)
Csat11=Csath11./(1+(((Csath11-Csatl11)./Csatl11)./(1-kesisi11+kesisi11.*sqrt(1-1i.*omegadim(q).*tao11./(kesisi11^2)))));
Csat22=Csath22./(1+(((Csath22-Csatl22)./Csatl22)./(1-kesisi22+kesisi22.*sqrt(1-1i.*omegadim(q).*tao22./(kesisi22^2)))));
Csat12=Csath12./(1+(((Csath12-Csatl12)./Csatl12)./(1-kesisi12+kesisi12.*sqrt(1-1i.*omegadim(q).*tao12./(kesisi12^2)))));
Csat66=Csath66./(1+(((Csath66-Csatl66)./Csatl66)./(1-kesisi66+kesisi66.*sqrt(1-1i.*omegadim(q).*tao66./(kesisi66^2)))));
MMou=((Csat11-Csat66).*sin(citaww)^2-(Csat22-Csat66).*cos(citaww).^2).^2+(Csat12+Csat66).^2.*sin(2.*citaww).^2;
Csatt(qqq,q)=(Csat11.*sin(citaww).^2+Csat22.*cos(citaww).^2+Csat66+sqrt(MMou))./2;
vsatt(qqq,q)=(real((1./sqrt((Csatt(qqq,q)./rho))))).^(-1);% The velocity of Branch function
QBranch(qqq,q)=imag(Csatt(qqq,q))./real(Csatt(qqq,q));% The attenuation of Branch function
Csattt(qqq,q)=vsatt(qqq,q).^2.*rho+1i.*vsatt(qqq,q).^2.*rho.*QBranch(qqq,q);
vsattt(qqq,q)=(((1./sqrt((Csatt(qqq,q)./rho))))).^(-1);% The velocity of Branch function
end
end
%%
%anisotropic parameters for P-wave
% Theoretical method
ees=real(CCsats(3,:)-CCsats(1,:))./(2.*real(CCsats(1,:)));
eeQs=1./2.*(Q1(3,:)-Q1(1,:));
detav=4.*((v3(2,:)./v3(1,:))-1)-((v3(3,:)./v3(1,:))-1);
detareal=real(detav);
detaimag=imag(detav);
%Interpolation method
eeb=real(Csatt(3,:)-Csatt(1,:))./(2.*real(Csatt(1,:)));
eeQb=1./2.*(QBranch(3,:)-QBranch(1,:));
detavb=4.*((vsattt(2,:)./vsattt(1,:))-1)-((vsattt(3,:)./vsattt(1,:))-1);
detabreal=real(detavb);
detabimag=imag(detavb);
%%
%%
figure(1); 
plot(log10(kk2.^2),log10(Q2(1,:)),'.')
hold on;
plot(log10(kk2.^2),log10(QBranch(1,:)),'k')
hold on
plot(log10(kk2.^2),log10(Q2(2,:)),'.')
hold on;
plot(log10(kk2.^2),log10(QBranch(2,:)),'b')
hold on
plot(log10(kk2.^2),log10(Q2(3,:)),'.')
hold on;
plot(log10(kk2.^2),log10(QBranch(3,:)),'r')
legend('Numerical simulations (  )','Interpolation approach (  )','Numerical simulations (  )','Interpolation approach (  )','Numerical simulations (  )','Interpolation approach (  )')
xlabel ('Frequency')
ylabel ('Log Attenuation')
grid on;
grid on;
figure(2); 
% plot(log10(kk2.^2),vG,'*');
% hold on;
plot(log10(kk2.^2),v2(1,:),'.')
hold on;
plot(log10(kk2.^2),vsatt(1,:),'b')
hold on
plot(log10(kk2.^2),v2(2,:),'.')
hold on;
plot(log10(kk2.^2),vsatt(2,:),'b')
hold on
plot(log10(kk2.^2),v2(3,:),'.')
hold on;
plot(log10(kk2.^2),vsatt(3,:),'b')
legend('Numerical simulations (  )','Interpolation approach (  )','Numerical simulations (  )','Interpolation approach (  )','Numerical simulations (  )','Interpolation approach (  )')
xlabel ('Frequency')
ylabel ('Velocity')
grid on;
grid on;
figure(3)
plot(log10(kk2.^2),ees,'.')
hold on;
plot(log10(kk2.^2),eeb,'b')
hold on
plot(log10(kk2.^2),detareal,'*')
hold on;
plot(log10(kk2.^2),detabreal,'r')
xlabel ('Frequency')
ylabel ('Anisotropy parameters')
legend('Numerical simulations (      )','Interpolation approach (      )','Numerical simulations (      )','Interpolation approach (      )')
figure(4)
plot(log10(kk2.^2),log10(eeQs),'.')
hold on;
plot(log10(kk2.^2),log10(eeQb),'b')
hold on
plot(log10(kk2.^2),log10(detaimag),'*')
hold on;
plot(log10(kk2.^2),log10(detabimag),'r')
xlabel ('Frequency')
ylabel ('Anisotropy parameters (lg)')
legend('Numerical simulations (      )','Interpolation approach (      )','Numerical simulations (      )','Interpolation approach (      )')
t1=cputime-t0;

