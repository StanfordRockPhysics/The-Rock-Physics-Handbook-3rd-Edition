function [k,mu,por]=dem(k1,mu1,k2,mu2,asp,phic)
%DEM - Effective elastic moduli using Differential Effective Medium
%      formulation.
%
%[K,MU,POR]=DEM(K1,MU1,K2,MU2,ASP,PHIC)
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

% Written by T. Mukerji

global DEMINPT;
DEMINPT=ones(1,6);

%mu2=3*k2*(1-2*nu2)/(2-2*nu2);

DEMINPT(1)=k1; 
DEMINPT(2)=mu1; 
DEMINPT(3)=k2; 
DEMINPT(4)=mu2; 
DEMINPT(5)=asp; 
DEMINPT(6)=phic; 

%save deminpt;

[tout, yout]=ode45m('demyprime',0.00,0.99999,[k1; mu1],1e-10);

 k=real(yout(:,1)); mu=real(yout(:,2));
por=phic*tout;
if nargout==0
plot(por,k,'-g',por,mu,'--r', 'linewidth', 1);
end;
clear DEMINPT;
