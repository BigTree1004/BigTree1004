function [x,V]=Spectrum(peak,sig,logx1,logx2,nx)
% calculate the radius specturm and the volume fraction
% 
% Inputs:
% peak = the radius which has the biggest volume fraction
% sig = the degree of dispersion
% loga1 = logarithm of the minimum radius
% loga2 = logarithm of the maximum radius
% 
% Outputs:
% all_a = the specturm of radius a
% V_a = the volume fraction of radius a 

x=logspace(logx1,logx2,nx);
V = 10.^normpdf(log10(x),log10(peak),sig);
V = V-min(V); V = V./sum(V(:));

end