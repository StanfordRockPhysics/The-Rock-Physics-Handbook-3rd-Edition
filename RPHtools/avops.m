function Rps=avops(vp1,vs1,d1,vp2,vs2,d2,ang,approx);
%Rps=AVOPS(vp1,vs1,d1,vp2,vs2,d2,ang,approx);
%
%Calculates P-to-S reflectivity (Rps) as a function of
%the angle of incidence (ang).
%input parameters:
%  layer 1 (top): vp1, vs1, density1 (d1)
%  layer 2 (bottom): vp2, vs2, density2 (d2)
% ang: vector with angles(DEG)
% approx: 1)Full Zoeppritz (A&R)
%	  2)Aki & Richards
%         3)Donati's 98 SEG paper (quadratic)
%         4)Donati's 98 SEGpaper (linear)
%         5)"max" (linear) simplification
%	  6)Ezequiel Gonzalez' approx
%         7)Alejandro & Reinaldo approx
%
% With no output arguments, plots Rps vs. angle.
%
% See also AVOPP, AVO_ABE

% written by Ezequiel Gonzalez (Oct,1999)

t=ang.*pi./180;	p=sin(t)./vp1;	ct=cos(t);
da=(d1+d2)/2;     Dd=(d2-d1);
vpa=(vp1+vp2)/2;  Dvp=(vp2-vp1);
vsa=(vs1+vs2)/2;  Dvs=(vs2-vs1);
cj1=sqrt(1-(sin(t).^2.*(vs1.^2./vp1.^2)));

switch approx
   case 1,		%FULL Zoeppritz (A&K)
	ct2=sqrt(1-(sin(t).^2.*(vp2.^2./vp1.^2)));
	cj2=sqrt(1-(sin(t).^2.*(vs2.^2./vp1.^2)));
	a=(d2.*(1-(2.*vs2.^2.*p.^2)))-(d1.*(1-(2.*vs1.^2.*p.^2)));
	b=(d2.*(1-(2.*vs2.^2.*p.^2)))+(2.*d1.*vs1.^2.*p.^2);
	c=(d1.*(1-(2.*vs1.^2.*p.^2)))+(2.*d2.*vs2.^2.*p.^2);
	d=2.*((d2.*vs2.^2)-(d1.*vs1.^2));
	E=(b.*ct./vp1)+(c.*ct2./vp2);
	F=(b.*cj1./vs1)+(c.*cj2./vs2);
	G=a-(d.*ct.*cj2./(vp1.*vs2));
	H=a-(d.*ct2.*cj1./(vp2.*vs1));
	D=(E.*F)+(G.*H.*p.^2);
	Rps=-2.*(ct./vp1).*...
		((a.*b)+(c.*d.*ct2.*cj2./(vp2.*vs2))).*p.*vp1./(vs1.*D);
   case 2,		%Aki & Richard (aprox)
%assuming (angles) i=i1, and j=j1
	Rps=(-p.*vpa./(2.*cj1)).*( ((Dd./da).*(1-(2.*vsa.^2.*p.^2)+ ...
		(2.*vsa.^2.*ct.*cj1./(vpa.*vsa)))) - ((Dvs./vsa).* ...
		((4.*vsa.^2.*p.^2)-(4.*vsa.^2.*ct.*cj1./(vpa.*vsa)))) );
   case 3,		%Donati's paper (quadratic)
	A0=-0.5.*(((Dd./da).*(1-(2.*vsa.^2./vpa.^2)))-(4.*vsa.*Dvs./vpa.^2));
	A1=-0.5.*(((Dd./da)+(2.*Dvs./vsa)).*((2.*vsa./vpa)-(vsa.^3./vpa.^3)));
	A2=-(vsa.^2./vpa.^2).*((Dd./da)+(2.*Dvs./vsa));
	Rps=sin(t).*(A0 + (A1.*ct) + (A2.*ct.^2));
   case 4,		%Donati's paper (linear)
	A0=-0.5.*(((Dd./da).*(1-(2.*vsa.^2./vpa.^2)))-(4.*vsa.*Dvs./vpa.^2));
	A1=-0.5.*(((Dd./da)+(2.*Dvs./vsa)).*((2.*vsa./vpa)-(vsa.^3./vpa.^3)));
	Rps=sin(t).*(A0 + (A1.*ct));
   case 5,		%Max. simplification
	Rps=-sin(t).*((Dd./(2.*da))-(((Dd./da)+(2.*Dvs./vsa)).*vsa/vpa));
   case 6,              %Ezequiel
	Rps=sin(t).*((-0.5.*Dd./da)-((vsa./vpa).*((Dd./da)+(2.*Dvs./vsa)))+...
	(((vsa./vpa).^3).*((0.5.*Dd./da)+(Dvs./vsa))));
    case 7,      %Alejandro & Reinaldo
        Rps=-2.*(vs1./vp1).*sin(t).*((Dd./da.*(0.5+(0.25.*vpa./vsa)))+...
        (Dvs./vsa));
   otherwise,	
end

if nargout==0
plot(ang,Rps)
end;
