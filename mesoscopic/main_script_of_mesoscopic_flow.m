% Script description
% Script to calculate: Complex Patchy Saturation Characteristics at Mesoscopic Scale 
% The file contains the basic rock physics parameters is data.mat,
% and can be changed manually.

clear all; close all;
addpath(genpath(pwd));

% load data
load('data_meso.mat');
load('color.mat');

d1=-4; d2=10; nf=201;   f=logspace(d1,d2,nf);   w=2*pi*f; %frequency

%------------------------------scenario 1----------------------------------%
scenario=1; a=8.3;   Sg=0.1;
[Echange1,K1,Vp1,Qpinv1]=DVS_of_Meso(w,K0,Rho0,perm,phi,tao,Kd,Mud,Pc,n1,Kf1,Rhof1,n2,Kf2,Rhof2,scenario,a,Sg);

%------------------------------scenario 2----------------------------------%
scenario=2; 

sig=0.7;    peak=8.3;%cm
[all_a,V_a]=Spectrum(peak,sig,-1,2,51);

[Echange2,K2,Vp2,Qpinv2]=DVS_of_Meso(w,K0,Rho0,perm,phi,tao,Kd,Mud,Pc,n1,Kf1,Rhof1,n2,Kf2,Rhof2,scenario,all_a,Sg,V_a);

%------------------------------scenario 3----------------------------------%
scenario=3; 

all_Sg=log10(all_a)+1;
all_Sg=all_Sg*Sg/sum(all_Sg.*V_a);

[Echange3,K3,Vp3,Qpinv3]=DVS_of_Meso(w,K0,Rho0,perm,phi,tao,Kd,Mud,Pc,n1,Kf1,Rhof1,n2,Kf2,Rhof2,scenario,all_a,all_Sg,V_a);


%------------------------------plot the results----------------------------%
line_width=4;

hfig=figure;set(hfig,'Color','w');
semilogx(f,real(Vp1),'color',color(1,:),'LineWidth',line_width);hold on;
semilogx(f,real(Vp2),'color',color(2,:),'LineWidth',line_width);hold on;
semilogx(f,real(Vp3),'color',color(3,:),'LineWidth',line_width);hold on;
ylabel('Vp (m/s)','FontSize',24,'FontWeight','demi');
set(gca,'FontSize',24,'LineWidth',2);
set(gca,'FontWeight','demi','box','on');
xlabel('Frequency (Hz)','FontSize',24,'FontWeight','demi');
xlim([f(1) f(nf)]);
set(gca,'xtick',[1e-4 1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
grid on;
legend('Scenario 1','Scenario 2','Scenario 3');

hfig=figure;set(hfig,'Color','w');
semilogx(f,real(Qpinv1),'color',color(1,:),'LineWidth',line_width);hold on;
semilogx(f,real(Qpinv2),'color',color(2,:),'LineWidth',line_width);hold on;
semilogx(f,real(Qpinv3),'color',color(3,:),'LineWidth',line_width);hold on;
ylabel('1/Q','FontSize',24,'FontWeight','demi');
set(gca,'FontSize',24,'LineWidth',2);
set(gca,'FontWeight','demi','box','on');
xlabel('Frequency (Hz)','FontSize',24,'FontWeight','demi');
xlim([f(1) f(nf)]);
ylim([-0.001 0.151]);
set(gca,'xtick',[1e-4 1e-2 1e0 1e2 1e4 1e6 1e8 1e10]);
grid on;


hfig=figure;set(hfig,'Color','w');
[hAx,hLine1,hLine2]=plotyy(all_a,V_a*100,all_a,all_Sg*100);
set(hLine1,'color','b','LineWidth',line_width);
set(hLine2,'color','g','linestyle','--','linewidth',line_width);
xlabel('Radius of gas patch a (cm)','fontsize',24,'color','k');
set(get(hAx(1),'ylabel'),'string', 'Volume Fraction (%)','fontsize',24,'color','k');
set(get(hAx(2),'ylabel'),'string', 'Gas saturation (%)','fontsize',24,'color','k');
set(hAx(1),'xscale','log','FontSize',20,'LineWidth',2);
set(hAx(2),'xscale','log','FontSize',20,'LineWidth',2);
grid on;
