function [k,mu,kv,muv,porv]=dem1(k1,mu1,k2,mu2,asp,phic,por)
%DEM1 - Effective elastic moduli using Differential Effective Medium
%      formulation. Returns effective moduli at the porosity POR
%      specified in the input.
%
%[K,MU]=DEM1(K1,MU1,K2,MU2,ASP,PHIC,POR)
%
%	K1, MU1:	Bulk and shear moduli of background matrix
%	K2, MU2:	Bulk and shear moduli of inclusions
%	ASP:		Aspect ratio of inclusions
%			<1 for oblate spheroids; >1 for prolate spheroids
%	PHIC:		percolation porosity for modified DEM model
%			=1 for usual DEM
%	K, MU:		Effective bulk and shear moduli
%	POR:		Porosity, fraction of phase 2.
%			For the modified DEM, where phase two is the
%			critical phase, POR is the actual porosity.

%Written by T. Mukerji, 1997

global DEMINPT;
DEMINPT=ones(1,6);

%mu2=3*k2*(1-2*nu2)/(2-2*nu2);

DEMINPT(1)=k1; 
DEMINPT(2)=mu1; 
DEMINPT(3)=k2; 
DEMINPT(4)=mu2; 
DEMINPT(5)=asp; 
DEMINPT(6)=phic; 

tfinal=por./phic;

%save deminpt;

%[tout, yout]=ode45m('demyprime',0.00,tfinal,[k1; mu1],1e-10);
[tout, yout]=ode45m('demyprime',0.00,tfinal,[k1; mu1],1e-5);

n=length(tout);
 k=real(yout(n,1)); mu=real(yout(n,2));
kv=real(yout(:,1)); muv=real(yout(:,2));
porv=phic.*tout;
%if nargout==0
%plot(por,k,'-g',por,mu,'--y', 'linewidth', 1);
%end;
clear DEMINPT;
