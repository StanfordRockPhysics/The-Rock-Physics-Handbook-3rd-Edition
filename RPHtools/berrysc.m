function [kbr,mubr,por]=berrysc(k1,mu1,k2,mu2,asp1,asp2)
%BERRYSC - Effective elastic moduli using Berryman's Self-Consistent
%(Coherent Potential Approximation) method. 
%
%[KBR,MUBR,POR]=BERRYSC(K1,MU1,K2,MU2,ASP1,ASP2)
%	K1,MU1,K2,MU2: Bulk and shear moduli of the two constituent
%		       phases
%	ASP1, ASP2:    Aspect ratio for the inclusions of the two phases
%			< 1 for oblate spheroids; >1 for prolate spheroids.
%	KBR,MUBR:  	Effective bulk and shear moduli 
%	POR:		Porosity (fraction of phase 2)
%
% See also BERRYSCM

% Written by T. Mukerji, 1994
%
%******* Berryman SC oblate and prolate spheroids *****************

kbr=[]; mubr=[]; por=[];

if asp1==1.
   asp1=0.99;
end
if asp2==1.
   asp2=0.99;
end
  
%

if asp1 < 1.
theta1=(asp1/((1-asp1^2)^(3/2)))*(acos(asp1) -asp1*sqrt(1-asp1^2));
fn1=(asp1^2/(1-asp1^2))*(3*theta1 -2);
end

if asp2 < 1.
theta2=(asp2/((1-asp2^2)^(3/2)))*(acos(asp2) -asp2*sqrt(1-asp2^2));
fn2=(asp2^2/(1-asp2^2))*(3*theta2 -2);
end

if asp1 > 1.
theta1=(asp1/((asp1^2-1)^(3/2)))*(asp1*sqrt(asp1^2-1)-acosh(asp1));
fn1=(asp1^2/(asp1^2-1))*(2-3*theta1);
end

if asp2 > 1.
theta2=(asp2/((asp2^2-1)^(3/2)))*(asp2*sqrt(asp2^2-1)-acosh(asp2));
fn2=(asp2^2/(asp2^2-1))*(2-3*theta2);
end

epsilon=1e-7;
for x1=0:0.01:1
if x1==0, x1=x1+epsilon;end;
if x1==1, x1=x1-epsilon;end;
 x2=1-x1;
% for x2=0.:0.001:0.6
%    x1=1-x2;

ksc= x1*k1 + x2*k2;
musc= x1*mu1 + x2*mu2;
knew= 0.;
munew= 0.;
%tol=real(ksc)/50000.0;
tol=1e-6*k1;
del=abs(ksc-knew);
niter=0;

while( (del > abs(tol)) & (niter<3000) )
	nusc=(3*ksc-2*musc)/(2*(3*ksc+musc));
	a1=mu1/musc -1; a2=mu2/musc -1;
	b1=(1/3)*(k1/ksc -mu1/musc); b2=(1/3)*(k2/ksc -mu2/musc);
	r=(1-2*nusc)/(2*(1-nusc));

	f11=1+a1*((3/2)*(fn1+theta1)-r*((3/2)*fn1+(5/2)*theta1-(4/3)));
	f12=1+a2*((3/2)*(fn2+theta2)-r*((3/2)*fn2+(5/2)*theta2-(4/3)));

	f21=1+a1*(1+(3/2)*(fn1+theta1)-(r/2)*(3*fn1+5*theta1))+b1*(3-4*r);
	f21=f21+(a1/2)*(a1+3*b1)*(3-4*r)*(fn1+theta1-r*(fn1-theta1+2*theta1^2));
	f22=1+a2*(1+(3/2)*(fn2+theta2)-(r/2)*(3*fn2+5*theta2))+b2*(3-4*r);
	f22=f22+(a2/2)*(a2+3*b2)*(3-4*r)*(fn2+theta2-r*(fn2-theta2+2*theta2^2));

	f31=1+a1*(1-(fn1+(3/2)*theta1)+r*(fn1+theta1));
	f32=1+a2*(1-(fn2+(3/2)*theta2)+r*(fn2+theta2));

	f41=1+(a1/4)*(fn1+3*theta1-r*(fn1-theta1));
	f42=1+(a2/4)*(fn2+3*theta2-r*(fn2-theta2));

	f51=a1*(-fn1+r*(fn1+theta1-(4/3))) + b1*theta1*(3-4*r);
	f52=a2*(-fn2+r*(fn2+theta2-(4/3))) + b2*theta2*(3-4*r);

	f61=1+a1*(1+fn1-r*(fn1+theta1))+b1*(1-theta1)*(3-4*r);
	f62=1+a2*(1+fn2-r*(fn2+theta2))+b2*(1-theta2)*(3-4*r);

	f71=2+(a1/4)*(3*fn1+9*theta1-r*(3*fn1+5*theta1)) + b1*theta1*(3-4*r);
	f72=2+(a2/4)*(3*fn2+9*theta2-r*(3*fn2+5*theta2)) + b2*theta2*(3-4*r);

	f81=a1*(1-2*r+(fn1/2)*(r-1)+(theta1/2)*(5*r-3))+b1*(1-theta1)*(3-4*r);
	f82=a2*(1-2*r+(fn2/2)*(r-1)+(theta2/2)*(5*r-3))+b2*(1-theta2)*(3-4*r);

	f91=a1*((r-1)*fn1-r*theta1) + b1*theta1*(3-4*r);
	f92=a2*((r-1)*fn2-r*theta2) + b2*theta2*(3-4*r);


	p1=3*f11/f21; p2=3*f12/f22;
	q1=(2/f31) + (1/f41) +((f41*f51 + f61*f71 - f81*f91)/(f21*f41));
	q2=(2/f32) + (1/f42) +((f42*f52 + f62*f72 - f82*f92)/(f22*f42));

	p1=p1/3; p2=p2/3;
	q1=q1/5; q2=q2/5;

%------------------------------------------------------------------------

	knew= (x1*k1*p1 + x2*k2*p2)/(x1*p1 + x2*p2);
	munew= (x1*mu1*q1 + x2*mu2*q2)/(x1*q1 + x2*q2);
	
	del=abs(ksc-knew);
	ksc=knew;
	musc=munew;
	niter=niter+1;
end		
% x1, niter  
%	vpbr=sqrt((ksc+ (4*musc/3))/ro);
%	vsbr=sqrt(musc/ro);
	kbr=[kbr,ksc]; mubr=[mubr, musc]; por=[por, 1-x1];
end
%end
if nargout==0
plot(por,real(kbr),'-g',por,real(mubr),'-c','linewidth',1)
end
%plot(por,real(mubr))
