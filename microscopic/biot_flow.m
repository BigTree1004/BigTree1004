function [Vp,Qinv]=biot_flow(K,G,phi,Kf,Kd,Ks,tao,yeta,perm_0,Rhof,Rho,om)
% calculate the Biot effect
% 
% Inputs:
% K = effective bulk modulus which contains the mesoscopic/microscopic dispersion
% G = effective shear modulus
% phi =  prosity of the rock
% Kf = effective fluid modulus
% tao = Tortuosity of pores
% yeta = effective viscosity
% perm_0 = low frequency permeability
% Rhof = effective fluid density
% Rho = effective bulk density
% om = 2*pi*frequency, is the circular frequency
% 
% Outputs:
% Vp = complex P wave velocity containing Biot effect
% Qinv = P wave attenuation

K=conj(K);  G=conj(G);

alpha=1-Kd/Ks;% Biot coefficient
beta=phi/Kf+(alpha-phi)/Ks;

temp4=sqrt(1-1i*1/2*tao*perm_0*Rhof*om/(yeta*phi));
temp5=1i*tao*perm_0*Rhof*om/(yeta*phi);
perm=perm_0./(temp4-temp5);
theta=1i.*perm./(yeta*om);
bs=Rhof.*theta.*om.^2;
b0=-beta.*(Kd+4*G/3+alpha^2/beta)/alpha;
c=(alpha-bs.*Rho./(Rhof*b0))./(alpha+bs);
temp6=sqrt((K+4/3*G)/Rho);
kp0=om./temp6;

b=1/2*b0.*(c-sqrt(c.^2-4*alpha*(1-c)./b0));
temp_11=(1+b*Rhof/Rho)./(1-b/b0);
kp_11=kp0.*sqrt(temp_11);
Vp=om./kp_11;
Qinv=2*imag(kp_11)./real(kp_11); 

end
