% Script description
% Script to calculate: Complex Crack Characteristics at Microscopic Scale
% The file contains the basic rock physics parameters is data.mat,
% and can be changed manually.

clear all;  close all;
addpath(genpath(pwd));

% load data
load('data_micro.mat');
load('color.mat');

d1=-4; d2=10; nf=201;   f=logspace(d1,d2,nf);   w=2*pi*f; %frequency
line_width=4;


%------------------------------scenario 1----------------------------------%
scenario=1;
alpha_c=0.0005; cr=0.2;

[Echange1,K1,Vp1,Qpinv1]=DVS_of_Micro(w,Km,Gm,Rhom,perm,phi,tao,Kd,Gd,Pc,yeta,Kf,Rhof,scenario,alpha_c,cr);
%------------------------------scenario 2----------------------------------%
scenario=2;
sig=0.3;    peak=alpha_c;
[all_alpha_c,V_alpha_c]=Spectrum(peak,sig,-5,-2,101); % unimodel spectrum

[Echange2,K2,Vp2,Qpinv2]=DVS_of_Micro(w,Km,Gm,Rhom,perm,phi,tao,Kd,Gd,Pc,yeta,Kf,Rhof,scenario,all_alpha_c,cr,V_alpha_c);

% plot the crack aspect ratio spectrum of scenario 2
hfig=figure;set(hfig,'Color','w');
semilogx(all_alpha_c,V_alpha_c*100,'k','LineWidth',line_width);hold on;
xlabel('Crack aspect ratio','fontsize',24,'color','k');
ylabel('Volume Fraction (%)','fontsize',24,'color','k');
xlim([all_alpha_c(1) all_alpha_c(end)]);
set(gca,'FontSize',20,'LineWidth',2);
grid on;

%------------------------------scenario 3----------------------------------%
scenario=3;
sig=0.3;    peak=alpha_c;   [all_alpha_c,V_alpha_c1]=Spectrum(peak,sig,-5,-2,101); % the first peak of bimodel spectrum
sig=0.3;    peak=0.00005;     [~,V_alpha_c2]=Spectrum(peak,sig,-5,-2,101);% the second peak of bimodel spectrum

V_alpha_c=V_alpha_c1/max(V_alpha_c1)+V_alpha_c2/max(V_alpha_c2);
V_alpha_c = V_alpha_c./sum(V_alpha_c);

[Echange3,K3,Vp3,Qpinv3]=DVS_of_Micro(w,Km,Gm,Rhom,perm,phi,tao,Kd,Gd,Pc,yeta,Kf,Rhof,scenario,all_alpha_c,cr,V_alpha_c);

% plot the crack aspect ratio spectrum of scenario 3
hfig=figure;set(hfig,'Color','w');
semilogx(all_alpha_c,V_alpha_c*100,'k','LineWidth',line_width);hold on;
xlabel('Crack aspect ratio','fontsize',24,'color','k');
ylabel('Volume Fraction (%)','fontsize',24,'color','k');
xlim([all_alpha_c(1) all_alpha_c(end)]);
set(gca,'FontSize',20,'LineWidth',2);
grid on;

%------------------------------plot the results----------------------------%
hfig=figure;set(hfig,'Color','w');
semilogx(f,real(Vp1),'color',color(1,:),'LineWidth',line_width);hold on;
semilogx(f,real(Vp2),'color',color(2,:),'LineWidth',line_width);hold on;
semilogx(f,real(Vp3),'color',color(3,:),'LineWidth',line_width);hold on;
ylabel('Vp (m/s)','FontSize',24,'FontWeight','demi');
set(gca,'FontSize',24,'LineWidth',2);
set(gca,'FontWeight','demi','box','on');
xlabel('Frequency (Hz)','FontSize',24,'FontWeight','demi');
xlim([f(1) 1e7]);
set(gca,'xtick',[1e-4 1e-2 1e0 1e2 1e4 1e6 1e8 1e9]);
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
xlim([f(1) 1e7]);
ylim([0 0.07]);
set(gca,'xtick',[1e-4 1e-2 1e0 1e2 1e4 1e6 1e8 1e9]);
grid on;

