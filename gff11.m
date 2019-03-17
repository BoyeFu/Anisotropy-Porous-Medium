%%
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% this is a function usd in scatteringproblem, it can calculate the N and S
% in equation (52)
% v.0.1, 18/10/2018, Boye Fu & Boris Gurevich, Curtin University and CRGC

%%
function g = gff11(z)
global I kk Kg mug Kf porosity alpha Kdry mudry taudry L Kstar M HBiot k1 citaw k2;
global Mdim HBiotdim Ldim;
x2=max(1000,kk*10);%the high limit of the integration, becuase the integration region is (0,infinity) in equation (A13),we use the value x2 represent infinity
x1=0.0;
NMAX=200;
% CALCULATE GAUSS-LEGENDRE WEIGHTS, In this paper, we use GAUSS-LEGENDRE
m=(NMAX+1)/2;%%£¿
xm=(x2+x1)/2;%midpoint 
xl=(x2-x1)/2;%midpoint
NN=xm;
MM=xl;
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
%%
FF=zeros(1,NMAX);
for i=1:NMAX;   
FF(i)=(lambdaK_generalfbUN(u(i),z)-lambdaK_generalfbUS(u(i),z)).*w(i);%The value of integral of the difference between N and S in equation (52)
end
%%
p12=sum(FF);
p=p12;
g=p;
end