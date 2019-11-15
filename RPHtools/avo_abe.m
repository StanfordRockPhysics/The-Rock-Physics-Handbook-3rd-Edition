function [A,B1,B2,E1,E2]=avo_abe(vp1,vs1,d1,vp2,vs2,d2);
%function [A,B1,B2,E1,E2]=avo_abe(vp1,vs1,d1,vp2,vs2,d2);
%
%Calculates AVO parameters:
%  A: Intercept (P-P) i.e. normal incidence reflectivity
%  B1, B2: P-P Gradient using different approximations
%  E1, E2: P-S Gradient using different approximations
%input parameters:
%  layer 1 (top): vp1, vs1, density1 (d1)
%  layer 2 (bottom): vp2, vs2, density2 (d2)
%  
%        B1 Shuey's paper (2terms->B Castag)
%        B2 Castagna's paper->Shuey
%(note: both are Shuey's approximation, but B2 is using Castagna's "way" to
%calculate them. The results are slightly different)
%        E1 Ezequiel Gonzalez approximation
%        E2 Alejandro & Reinaldo approx
%
% See also AVOPP, AVOPS

% written by Ezequiel Gonzalez (Oct,1999)
% modified T. Mukerji Feb 2001

da=(d1+d2)/2;     Dd=(d2-d1);
vpa=(vp1+vp2)/2;  Dvp=(vp2-vp1);
vsa=(vs1+vs2)/2;  Dvs=(vs2-vs1);
Ro=0.5.*((Dvp./vpa)+(Dd./da));
A=Ro;

%%   case 1,		%Shuey's paper (2terms->B Castag)
	poi1=((0.5.*(vp1./vs1).^2)-1)./((vp1./vs1).^2-1);
	poi2=((0.5.*(vp2./vs2).^2)-1)./((vp2./vs2).^2-1);
	poia=(poi1+poi2)./2;   Dpoi=(poi2-poi1);
	Bx=(Dvp./vpa)./((Dvp./vpa)+(Dd./da));
	Ax=Bx-(2.*(1+Bx).*(1-2.*poia)./(1-poia));
	B1=(Ax.*Ro)+(Dpoi./(1-poia).^2);
%%   case 2,		%Castagna's paper->Shuey
	B2=(-2.*vsa.^2.*Dd./(vpa.^2.*da)) + (0.5.*Dvp./vpa) - ...
	   (4.*vsa.*Dvs./(vpa.^2));

%% E1 Gonzalez approx.
E1=(-0.5.*Dd./da)-((vsa./vpa).*((Dd./da)+(2.*Dvs./vsa))) + ...
  (((vsa./vpa).^3).*((0.5.*Dd./da)+(Dvs./vsa)));

%% E2 Alejandro & Reinaldo
E2=-2.*(vs1./vp1).*((Dd./da.*(0.5+(0.25.*vpa./vsa)))+...
        (Dvs./vsa));

