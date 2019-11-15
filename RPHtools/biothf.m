function [vp1,vp2,vs]=biothf(vpdry,vsdry,k0,mu0,ro0,rofl,kfl,por,alfa)
%[VP1,VP2,VS]=BIOTHF(VPDRY,VSDRY,K0,MU0,RO0,ROFL,KFL,POR,ALFA)
%	The high frequency limiting velocities predicted by Biot theory
%	as given by Johnson and Plona.
%Inputs: VPDRY, VSDRY: P and S velocities of dry rock
%	 K0, MU0, RO0: Mineral bulk and shear moduli, and density
%	 ROFL,KFL: Pore fluid density and bulk modulus
%	 POR: Porosity
%	 ALFA: Tortuosity parameter in Biot theory. Always > 1. Usually 1 to 3.
%	VPDRY, VSDRY can be single scalar values or a vector of inputs.
%	If they are vectors, then all other inputs should either be vectors
%	of the same length as VPDRY and VSDRY or they should be scalars.
%Outputs: VP1, VP2, VS: The high frequency limits of the fast and slow P-waves
%	  and shear wave velocities respectively
%
%	See also BIOT, BIOTHFB

%Written by T. Mukerji 1996

if length(k0)~=length(vpdry), k0=k0*ones(size(vpdry)); end
if length(mu0)~=length(vpdry), mu0=mu0*ones(size(vpdry)); end
if length(ro0)~=length(vpdry), ro0=ro0*ones(size(vpdry)); end
if length(rofl)~=length(vpdry), rofl=rofl*ones(size(vpdry)); end
if length(kfl)~=length(vpdry), kfl=kfl*ones(size(vpdry)); end
if length(por)~=length(vpdry), por=por*ones(size(vpdry)); end
if length(alfa)~=length(vpdry), alfa=alfa*ones(size(vpdry)); end

rodry=(1-por).*ro0; ro=(1-por).*ro0+por.*rofl;
mudry=rodry.*vsdry.^2; kdry=rodry.*vpdry.^2-(4/3)*mudry; b=kdry./k0;
robiot=ro0.*(1-por)+por.*rofl.*(1-1./alfa);
den=1-por-b+por.*k0./kfl;
p=((1-por).*(1-por-b).*k0+ por.*(k0./kfl).*kdry)./den + (4/3)*mudry;
q=((1-por-b).*por.*k0)./den;
r=por.^2.*k0./den;
ro12=(1-alfa).*por.*rofl; ro11=(1-por).*ro0-ro12; ro22=por.*rofl.*alfa;
vs=sqrt(mudry./robiot);
t1=p.*ro22+r.*ro11-2*q.*ro12; t2=p.*r-q.^2; t3=ro11.*ro22-ro12.^2;
vp1sqr=(t1+sqrt(t1.^2-4*t2.*t3))./(2*t3); vp1=sqrt(vp1sqr);
vp2sqr=(t1-sqrt(t1.^2-4*t2.*t3))./(2*t3); vp2=sqrt(vp2sqr);
