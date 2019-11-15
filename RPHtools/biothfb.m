function [vp1,vs]=biothfb(vpdry,vsdry,k0,mu0,ro0,rofl,kfl,por,alfa)
%[VP1,VS]=BIOTHFB(VPDRY,VSDRY,K0,MU0,RO0,ROFL,KFL,POR,ALFA)
%       The approximate high frequency limiting velocities predicted 
%	by Biot theory as given by Geertsma-Smit and recast by Bourbie et al.
%	This gives higher velocities than the actual high frequency limiting
%	velocities.
%Inputs: VPDRY, VSDRY: P and S velocities of dry rock
%        K0, MU0, RO0: Mineral bulk and shear moduli, and density
%        ROFL,KFL: Pore fluid density and bulk modulus
%        POR: Porosity
%        ALFA: Tortuosity parameter in Biot theory. Always > 1. Usually 1 to 3.
%       VPDRY, VSDRY can be single scalar values or a vector of inputs.
%       If they are vectors, then all other inputs should either be vectors
%       of the same length as VPDRY and VSDRY or they should be scalars.
%Outputs: VP1, VS: The high frequency limits of the fast P-wave
%         and shear wave velocities respectively
%
%	See also BIOT, BIOTHF 

% Written by T. Mukerji, 1996
 
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
t1=por.*ro./(rofl.*alfa)+(1-b).*(1-b-2*por./alfa); t2=(1-b-por)./k0 + por./kfl; 
vp1=sqrt((kdry+(4/3)*mudry+t1./t2)./robiot);
vs=sqrt(mudry./robiot);
