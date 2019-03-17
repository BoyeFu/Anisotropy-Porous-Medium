%%
% Scattering problem: Seismic dispersion and attenuation in saturated porous rock with aligned slit cracks
% this is a function usd in scatteringproblem, it can calculate the f(w)
% combining with gf222
% v.0.1, 18/05/2018, Boye Fu & Boris Gurevich, Curtin University and CRGC

%%
function g = gff22(z)
global I kk Kg mug Kf porosity alpha Kdry mudry taudry L Kstar M HBiot k1 citaw k2;
global Mdim HBiotdim Ldim;
p=-(HBiotdim-alpha*Mdim-2.*(1-sin(citaw).^2)).*real(exp(I.*(k1.*(z).*cos(citaw))));%the f in equation A13

g=p;
end