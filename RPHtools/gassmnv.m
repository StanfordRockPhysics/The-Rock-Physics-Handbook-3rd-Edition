function [vp2,vs2,ro2,k2]=gassmnv(vp1,vs1,ro1,rofl1,kfl1,rofl2,kfl2,k0,phi)
% [VP2,VS2,RO2,K2]=GASSMNV(VP1,VS1,RO1,ROFL1,KFL1,ROFL2,KFL2,K0,PHI)
%
% Gassmann fluid substituion with velocities as input/outputs
% VP1, VS1, RO1: rock Vp, Vs, and density with fluid 1
% ROFL1, KFL1:   density and bulk modulus of initial fluid
% ROFL2, KFL2:   density and bulk modulus of new fluid
% K0, PHI:       mineral bulk modulus, and rock porosity
% VP2,VS2,RO2, K2:  Vp, Vs, density, and bulk modulus of rock with new fluid
% Giving VS1=0, and mineral P-wave modulus in place of K0 does approximate 
% Gassmann calculation.
%
% See also GASSMNK

%Written by T. Mukerji

ro2=ro1 - phi.*rofl1 +phi.*rofl2;
mu1=ro1.*vs1.^2; k1=ro1.*vp1.^2-(4/3)*mu1;
a= k1./(k0-k1) - kfl1./(phi.*(k0-kfl1)) + kfl2./(phi.*(k0-kfl2)); 
k2= k0.*a ./(1+a); mu2=mu1;
vp2=sqrt((k2+(4/3)*mu2)./ro2); vs2=sqrt(mu2./ro2);
