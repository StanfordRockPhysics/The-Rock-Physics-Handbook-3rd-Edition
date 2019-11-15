function [vp,vs,ro,phi,c]=hertzmindv(vpmin,vsmin,romin,pres,phi,c)
%[VP,VS,RO,PHI,C]=HERTZMINDV(VPMIN,VSMIN,ROMIN,P,PHI,C)
%calculates P and S velocity of dry sphere packs under hydrostatic 
%pressure condition using Hertz-Mindlin model.
%
%inputs:
%       VPMIN, VSMIN, ROMIN    :       mineral P and S velocities, and density
%       P       :       pressure in Pa.
%       PHI     :       porosity (optional), default [0.2:0.05:0.7].
%       C       :       coordination number (optional, see below).
%outputs:
%       VP, VS, RO      : dry bulk P and S velocities, and density 
%       PHI     :       porosity.
%       C       :       coordination number
%
%       The program use the porosity-coordination number relation in page150 of 
%       the Rock Physics Handbook, if C is omitted.
%
%With no output arguments, plots modulus vs. porosity plot.
%See also HERTZMIND

%Written by Isao Takahashi 4/12/00
%modifications for velocity by T. Mukerji 5/2000
%corrections: k, g, denominator (1-nu)^2 (RPH pg. 151) Gary/Jack/Tapan
%7/2005


if nargin==4,
phi=[.2:.1:.7]';
phi=[.2:.05:.7]';
end
if nargin<=5,
ctemp=[14.007 12.336 10.843 9.5078 8.3147 7.2517 6.3108 5.4878 4.7826 4.1988 3.7440];
por=[.2:.05:.7];
c=interp1(por,ctemp,phi);
end;

[kmin, gmin]=v2ku(vpmin, vsmin, romin);

numin=(3*kmin-2*gmin)/(6*kmin+2*gmin);
k=((c.^2.*(1-phi).^2.*gmin.^2)./(18*pi^2*(1-numin).^2).*pres).^(1/3);
g=(5-4*numin)/(5*(2-numin))*((3*c.^2.*(1-phi).^2.*gmin.^2)./(2*pi^2*(1-numin).^2).*pres).^(1/3);

ro = (1-phi).*romin;
[vp, vs]=ku2v(k,g, ro);

if nargout==0
plot(phi,vp,'-g', phi,vs,'-c');
end;
