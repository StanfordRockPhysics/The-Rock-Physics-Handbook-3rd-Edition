function [k,g,phi,c]=hertzmind(kmin,gmin,pres,phi,c)
%[K,G,PHI,C]=HERTZMIND(KMIN,GMIN,P,PHI,C)
%calculates bulk and shear moduli of dry sphere packs under hydrostatic 
%pressure condition using Hertz-Mindlin model.
%
%inputs:
%       KMIN    :       mineral bulk modulus.
%       GMIN    :       mineral shear modulus.
%       P       :       pressure in Pa.
%       PHI     :       porosity (optional), default [0.2:0.05:0.7].
%       C       :       coordination number (optional, see below).
%outputs:
%       K      :       	bulk modulus.
%       G      :	shear modulus. 
%       PHI     :       porosity.
%       C       :       coordination number
%
%       The program use the porosity-coordination number relation in page150 of 
%       the Rock Physics Handbook, if C is omitted.
%
%With no output arguments, plots modulus vs. porosity plot.

%Written by Isao Takahashi 4/12/00
%correction: denominator for k, g corrected to (1-nu)^2 (RPH pg 151);
%Gary/Jack/Tapan 7/2005


if nargin==3,
phi=[.2:.1:.7]';
phi=[.2:.05:.7]';
end
if nargin<=4,
ctemp=[14.007 12.336 10.843 9.5078 8.3147 7.2517 6.3108 5.4878 4.7826 4.1988 3.7440];
por=[.2:.05:.7];
c=interp1(por,ctemp,phi);
end;

numin=(3*kmin-2*gmin)/(6*kmin+2*gmin);
k=((c.^2.*(1-phi).^2.*gmin.^2)./(18*pi^2*(1-numin).^2).*pres).^(1/3);
g=(5-4*numin)/(5*(2-numin))*((3*c.^2.*(1-phi).^2.*gmin.^2)./(2*pi^2*(1-numin).^2).*pres).^(1/3);


if nargout==0
plot(phi,k,'-g', phi,g,'-c');
end;
