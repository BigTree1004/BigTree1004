function [Echange,K,Vp,Qpinv]=DVS_of_Micro(w,Ks,Gs,Rhos,perm,phi,tao,Kd,Gd,Pc,yeta,Kf,Rhof,scenario,all_alpha_c,cr,V_alpha_c)
% calculate the microscopic dispersion characteristics
% 
% Inputs:
% w = 2*pi*frequency, is the angular frequency
% Ks =  bulk modulus of the grain
% Gs = shear modulus of the grain
% Rhos =  density of the grain
% perm =  permeability of the rock
% phi =  prosity of the rock
% tao = Tortuosity of pores
% Kd = bulk modulus of dry rock
% Gd = shear modulus of dry rock
% Pc = confining pressure
% yeta = viscosity of fluid
% Kf = bulk modulus of fluid
% Rhof = density of fluid
% scenario = 1,2,3
%     1: single crack aspect ratio
%     2: unimodel crack aspect ratio spectrum
%     3: bimodel crack aspect ratio spectrum
% all_alpha_c = crack aspect ratio or crack aspect ratio spectrum
% cr = crack density
% V_alpha_c = volume fraction for spectrum
% 
% Outputs:
% DVS = dynamic volumetric strain
% K = complex bulk modulus
% Vp = complex P wave velocity
% Qpinv = P wave attenuation

if scenario==1,V_alpha_c=1;end

if nargin<17,V_alpha_c=1;end

nu=(3*Ks-2*Gs)/(6*Ks+2*Gs); %  Poisson ratio

alpha=1-Kd/Ks;% Biot coefficient
K_Gass=Kd+alpha^2/((alpha-phi)/Ks+phi/Kf);% Gassmann's bulk modulus
rho=Rhos*(1-phi)+phi*Rhof;% effective bulk density

Echange=zeros(1,length(w));

for i=1:length(all_alpha_c)
    
alpha_c=all_alpha_c(i);

zeta=sqrt((3*1i*w*yeta)/(alpha_c^2*Kf)); 
temp1=(1/Kd-1/Ks)/(1/Kd-1/K_Gass);
temp2=(1-nu)/Gs;
temp3=4*(1-nu)*Kf/(3*Gs*alpha_c);
 
a11=1/Kd;
a12=alpha/Kd;
M1=(alpha-phi)/Ks+phi/Kf;
B1=alpha/(Kd*M1+alpha^2);
a22=a12/B1;
 
 for j=1:1:length(w)
f(j)=2*besselj(1,zeta(j))/(zeta(j)*besselj(0,zeta(j)));
S(j)=8/3*pi*cr*temp2*f(j)*(temp1-f(j))/(1+temp3*(1-f(j)));
dinvK(j)=a11-a12^2./(a22-S(j));
 end
   
dinvK=conj(dinvK);
DVS=dinvK*Pc;
temp_Echange=DVS(1)-DVS;

Echange=Echange+V_alpha_c(i)*temp_Echange;
end

K=1./(1/K_Gass+Echange/Pc);% effective bulk modulus which contains the microscopic dispersion
G=1./(1/Gd+4/15*Echange/Pc);% effective shear modulus 

[Vp,Qpinv]=biot_flow(K,G,phi,Kf,Kd,K_Gass,tao,yeta,perm,Rhof,rho,w);

end
