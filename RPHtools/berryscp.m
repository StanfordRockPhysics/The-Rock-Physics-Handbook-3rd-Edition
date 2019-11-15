function [kbr, mubr]=berryscp(k, mu, asp, x, p)
%BERRYSCP - Effective elastic moduli vs. pressure for multi-component 
% composite using Berryman's Self-Consistent method. Pressure dependence
% is determined by thinning and closing any penny-shaped components that
% are filled with a liquid or gas (shear modulus = 0).  Components that are
% solids or whose aspect ratio is greater than 0.2 (arbitrary separation of 
% stiffer pores) are assumed to not change with stress.
%
%[KBR,MUBR]=BERRYSCP(K,MU,ASP,X,P)
%	K,MU:          Bulk and shear moduli of the N constituent
%		       phases (K, MU, vectors of length N)
%	ASP:           Aspect ratio for the inclusions of the N phases
%		       < 1 for oblate spheroids; >1 for prolate spheroids.
%       X:             Fraction of each phase. Sum(X) should be 1.
%	KBR,MUBR:      Effective bulk and shear moduli 
%       P:             range of effective pressure values to compute.  Should
%                      be in same units as moduli.
%
% Note: k,mu,asp,x must be vectors of the same length (N)
%       p is a vector of any length (Np)
%       kbr,mubr are output vectors of length Np (same as the pressure)
%
% The program plots kbr and mubr as a function of pressure when called
% without the output arguments (kbr and mubr).
%
% See also BERRYSC, BERRYSCM

% Written by G. Mavko, 1999
% Additions by Diana Sava, 1999


%******* Uses Berryman SC oblate and prolate spheroids *****************

%force inputs to all be column vectors
k=k(:); mu=mu(:); asp=asp(:); x=x(:);

%force aspect=1 to be .99, so that general spheroid eqn works
indx=find(asp==1); asp(indx)=0.99;

% assume that mineral modulus is the largest of the input moduli
indx = find(k==max(k));
kmin = k(indx);
mumin = mu(indx);
prmin = (3*kmin - 2*mumin)/(6*kmin + 2*mumin);

% loop over pressures
for kk = 1:length(p)
	%compute change in aspect ratio for fluid-filled cracks with asp <.2
    	delasp = p(kk)*2*(1 - prmin)/(pi*mumin);
	delasp = (asp < .2).*(mu == 0)*delasp;
	
	%make sure that change does not exceed original aspect ratio
	delasp = min(delasp, asp);
	
	%thin the cracks accordingly and reduce their concentrations
	aspk = asp - delasp;
	xk = x.*(1 - delasp./asp);
	
	%renormalize the concentrations
	xk = xk./sum(xk);

	%take only the constituents with non-zero aspect ratio
	indx = find(aspk ~= 0);
	kin = k(indx);
	muin = mu(indx);
	aspkin = aspk(indx);
	xkin   = xk(indx);
	%compute self consistent moduli
	[kbr(kk),mubr(kk)]=berryscm(kin,muin,aspkin,xkin);
end

if nargout==0
	plot(p, kbr,'r',p,mubr,'b','linewidth',1);
end
