function [cti]=echeng(c0,phi,a,kfl)
%function [cti]=echeng(c0,phi,a,kfl)
%Calculates the TI Cijkl, CTI, for cracked rock (single aligned crack set), 
%with crack porosity PHI, aspect ratio A. KFL is the modulus of 
%fluid in the cracks. The stiffnesses (5 values) are arranged as:
%c = [c11 c13 c33 c44 c66]
%C0 is the background matrix stiffnesses. For isotropic background
%c11=lambda+2*mu, c13=lambda, c33=c11, and c44=c66=mu
%The calculation uses the Eshelby-Cheng formulation which is valid for
%all aspect ratios.

%Written by T. Mukerji 1997
%Ref. "The Rock Physics Handbook"

c11=c0(:,1); c13=c0(:,2); c33=c0(:,3); c44=c0(:,4); c66=c0(:,5);
lam=c0(:,2); mu=c0(:,4); k=lam+(2/3)*mu; c=kfl./(3*(k-kfl));

sig=(3*k-2*mu)./(6*k+2*mu);
r=(1-2*sig)./(8*pi*(1-sig));q=3*r./(1-2*sig);
sa=sqrt(1-a.^2);
ia=2*pi*a.*(acos(a)-a.*sa)./sa.^3;
ic=4*pi-2*ia; iac=(ic-ia)./(3*sa.^2); iaa=pi-(3/4)*iac; iab=iaa./3;
s1313 = 0.5*q.*iac.*(1-a.^2) + 0.5*r.*(ia+ic);
s1212 = q.*iab+r.*ia;
s31 = q.*iac-r.*ic;
s13 = q.*iac.*a.^2 - r.*ia;
s12 = q.*iab - r.*ia;
s33 = q.*(4*pi/3 - 2*iac.*a.^2) + ic.*r;
s11 = q.*iaa + r.*ia;
e=s33.*s11-s31.*s13-(s33+s11-2*c-1)+c.*(s31+s13-s11-s33);
d=s33.*s11+s33.*s12-2*s31.*s13-(s11+s12+s33-1-3*c)-c.*(s11+s12+2*(s33-s13-s31));

cti11=lam.*(s31-s33+1)+2*mu.*e./(d.*(s12-s11+1));
cti13=( (lam+2*mu).*(s13+s31)-4*mu.*c+lam.*(s13-s12-s11-s33+2) )./(2*d);
cti33=( (lam+2*mu).*(-s12-s11+1)+2*lam.*s13+4*mu.*c )./d;
cti44=mu./(1-2*s1313);
cti66=mu./(1-2*s1212);

cti11=c11-phi.*cti11;
cti13=c13-phi.*cti13;
cti33=c33-phi.*cti33;
cti44=c44-phi.*cti44;
cti66=c66-phi.*cti66;

cti=[cti11 cti13 cti33 cti44 cti66];
