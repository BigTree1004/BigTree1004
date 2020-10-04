This program package is used to study the effects of complex heterogeneities on wave dispersion and attenuation at micro- and meso-scopic scales, which is based on the model proposed in the paper "Extended Gassmann Equation with Dynamic Volumetric Strain: Modeling Wave Dispersion and Attenuation of Complex Heterogenous Porous Rocks".
The detailed descriptions of each program are as follows:

1. mesoscopic：

main_script_of_mesoscopic_flow.m studies the influence of complex patchy saturation characteristics at mesoscopic scale, which corresponds to Figures 6 and 7 in the paper. Note that the data file data_meso.mat contains the basic rock physics parameters, and the data file color.mat contains the colors which are used for plotting.
DVS_of_Meso.m is a subprogram which is used to calculate the DVS, bulk modulus, P wave velocity and attenuation results at mesoscopic scale.
Spectrum.m is a subprogram which is used to calculate gas radius spectrum at mesoscopic scale.
biot_flow.m is a subprogram which is used to incorporate the influence of macroscopic Biot flow with the mesoscopic (or microscopic) flow.

2. microscopic：

main_script_of_microscopic_flow.m studies the influence of complex crack characteristics at microscopic scale, which corresponds to Figures 8 and 9 in the paper. Note that the data file data_micro.mat contains the basic rock physics parameters, and the data file color.mat contains the colors which are used for plotting.
DVS_of_Micro.m is a subprogram which is used to calculate the DVS, bulk modulus, P wave velocity and attenuation results at microscopic scale.
Spectrum.m is a subprogram which is used to calculate crack aspect ratio spectrum at microscopic scale.
biot_flow.m is a subprogram which is used to incorporate the influence of macroscopic Biot flow with the microscopic (or mesoscopic) flow.
