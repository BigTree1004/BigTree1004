function [Echange,K,Vp,Qpinv]=DVS_of_Meso(w,K0,Rho0,perm,phi,tao,Kd,Mud,Pc,n1,Kf1,Rhof1,n2,Kf2,Rhof2,scenario,all_a,all_Sg,V_a)
% calculate the mesoscopic dispersion characteristics
% 
% Inputs:
% w = 2*pi*frequency, is the angular frequency
% K0 =  bulk modulus of the grain
% Rho0 =  density of the grain
% perm =  permeability of the rock
% phi =  prosity of the rock
% tao = Tortuosity of pores
% Kd = bulk modulus of dry rock
% Mud = shear modulus of dry rock
% Pc = confining pressure
% n1, n2 = viscosity of fluid 1, gas, and fluid 2, water
% Kf1, Kf2 = bulk modulus of gas and water
% Rhof1, Rhof2 = density of gas and water
% scenario = 1,2,3
%     1: single gas radius and single gas saturation
%     2: gas radius spectrum and single gas saturation
%     3: gas radius spectrum and gas saturation spectrum
% all_a = gas radius or gas radius spectrum
% all_Sg = gas saturation or gas saturation spectrum
% V_a = volume fraction for spectrum
% 
% Outputs:
% DVS = dynamic volumetric strain
% K = complex bulk modulus
% Vp = complex P wave velocity
% Qpinv = P wave attenuation


% set all_Sg and V_a for scenario 1 and 2
switch scenario
    case 1, V_a=1;
    case 2, all_Sg=all_Sg*ones(1,length(all_a));
end

if nargin<17,V_a=1;end

% calculate the basic rock physics parameters
tem1=(1-Kd/K0)^2; tem2=phi/Kf1+(1-phi)/K0-Kd/K0^2;
Ks1=Kd+tem1/tem2;  Mus1=Mud;
tem1=(1-Kd/K0)^2; tem2=phi/Kf2+(1-phi)/K0-Kd/K0^2;
Ks2=Kd+tem1/tem2;  Mus2=Mud;

y1=1-Kd/K0;
D1=K0/2/(y1+phi*(K0-Kf1)/Kf1);
y2=1-Kd/K0;
D2=K0/2/(y2+phi*(K0-Kf2)/Kf2);

H1=Ks1+4*Mus1/3;
H2=Ks2+4*Mus2/3;
la1=Kd-2*Mus1/3+2*y1^2*D1;
la2=Kd-2*Mus2/3+2*y2^2*D2;

KA1=(phi/Kf1+(1-phi)/K0-Kd/K0^2)^-1;
KA2=(phi/Kf2+(1-phi)/K0-Kd/K0^2)^-1;
KE1=(1-(Kf1*(1-Ks1/K0)*(1-Kd/K0))/(phi*(Ks1)*(1-Kf1/K0)))*KA1;
KE2=(1-(Kf2*(1-Ks2/K0)*(1-Kd/K0))/(phi*(Ks2)*(1-Kf2/K0)))*KA2;
Q1=(1-Kd/K0)*KA1/(Ks1);
Q2=(1-Kd/K0)*KA2/(Ks2);

alpha1=(1i*w*n1/perm/KE1).^0.5;
alpha2=(1i*w*n2/perm/KE2).^0.5;

Echange=zeros(1,length(w));

for i=1:length(all_a)

a=all_a(i)*1e-2; %m
b=a/all_Sg(i)^(1/3);

A=[
    a            -a          -1/a^2;
    H1+2*la1     -H2-2*la2   2*(H2-la2)/a^3;
    0            H2+2*la2    -2*(H2-la2)/b^3;
];

B=[0;0;-Pc];    X=A\B;
A1=X(1);    A2=X(2);    B2=X(3);

% calculate the high limit pore pressure in region 1 and 2
ph1=-6*y1*D1*A1;
ph2=-6*y2*D2*A2;

Z1=n1*a/perm*(1-exp(-2*alpha1*a))./((alpha1*a-1)+(alpha1*a+1).*exp(-2*alpha1*a));
Z2=-1*n2*a/perm*((alpha2*b+1).*exp(-2*alpha2*(b-a))+(alpha2*b-1))./...
    ((alpha2*b+1).*(alpha2*a-1).*exp(-2*alpha2*(b-a))-(alpha2*b-1).*(alpha2*a+1));

v=(ph1-ph2)./(Z1+Z2);

DVS=4*pi*a^2*v*(Q2-Q1)./(1i*w)/(4*pi*b^3/3);
DVS(isnan(DVS))=0;

temp_Echange=DVS(1)-DVS;

Echange=Echange+V_a(i)*temp_Echange;

end

Sg=sum(all_Sg.*V_a);

Kf=1/(Sg/Kf1+(1-Sg)/Kf2); % Wood's equation
alpha=1-Kd/K0; %Biot coefficient
K_Gass=Kd+alpha^2/((alpha-phi)/K0+phi/Kf); % Gassmann's bulk modulus

K=1./(1/K_Gass+Echange/Pc);% effective bulk modulus which contains the mesoscopic dispersion

n=n2*(n1/n2)^Sg; % effective viscosity
Rhof=Rhof1*Sg+Rhof2*(1-Sg); % effective fluid density
rho=(Rho0*(1-phi)+Rhof1*phi)*Sg+(Rho0*(1-phi)+Rhof2*phi)*(1-Sg); % effective bulk density

[Vp,Qpinv]=biot_flow(K,Mud,phi,Kf,Kd,K0,tao,n,perm,Rhof,rho,w);

end
