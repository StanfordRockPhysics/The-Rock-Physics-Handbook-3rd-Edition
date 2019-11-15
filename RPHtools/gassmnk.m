function k2=gassmnk(k1,kfl1,kfl2,k0,phi)
%K2=GASSMNK(K1,KFL1,KFL2,K0,PHI)
%
% Fluid substitution using Gassmann equation
% K1: original bulk modulus of rock saturated with fluid of bulk modulus KFL1
% K2: rock bulk modulus with new fluid of bulk modulus KFL2
% K0: mineral modulus
% PHI: porosity
%
% Giving the P-wave modulus in place of the bulk modulus does the approximate
% Gassmann calculation.
%
% See also GASSMNV

%Written by T. Mukerji

a= k1./(k0-k1) - kfl1./(phi.*(k0-kfl1)) + kfl2./(phi.*(k0-kfl2)); 
k2= (k0.*a ./(1+a)).*(phi~=0) + k1.*(phi==0); 
