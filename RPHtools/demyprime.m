function yprime=demyprime(t,y)
%function yprime=demyprime(t,y)
%used by DEM

%Written by T. Mukerji

%load deminpt;
global DEMINPT;

k1=DEMINPT(1); mu1=DEMINPT(2); k2=DEMINPT(3); mu2=DEMINPT(4);
asp=DEMINPT(5); phic=DEMINPT(6);

krc=k1*k2/((1-phic)*k2 + phic*k1);
murc=mu1*mu2/((1-phic)*mu2 + phic*mu1);

ka=k2; mua=mu2;
%ka=krc; mua=murc;

k=y(1); mu=y(2);
yprime=zeros(2,1);

if asp==1.
   asp=0.99;
end
  
%******* P and Q *****************
%

if asp < 1.
theta=(asp/((1-asp^2)^(3/2)))*(acos(asp) -asp*sqrt(1-asp^2));
fn=(asp^2/(1-asp^2))*(3*theta -2);
end

%if asp2 < 1.
%theta2=(asp2/((1-asp2^2)^(3/2)))*(acos(asp2) -asp2*sqrt(1-asp2^2));
%fn2=(asp2^2/(1-asp2^2))*(3*theta2 -2);
%end

if asp > 1.
theta=(asp/((asp^2-1)^(3/2)))*(asp*sqrt(asp^2-1)-acosh(asp));
fn=(asp^2/(asp^2-1))*(2-3*theta);
end

%if asp2 > 1.
%theta2=(asp2/((asp2^2-1)^(3/2)))*(asp2*sqrt(asp2^2-1)-acosh(asp2));
%fn2=(asp2^2/(asp2^2-1))*(2-3*theta2);
%end


	nu=(3*k-2*mu)/(2*(3*k+mu));
	r=(1-2*nu)/(2*(1-nu));
	a=mua/mu -1; 
%	a2=mu2/musc -1;
	b=(1/3)*(ka/k -mua/mu); 
%	b2=(1/3)*(k2/ksc -mu2/musc);

	f1a=1+a*((3/2)*(fn+theta)-r*((3/2)*fn+(5/2)*theta-(4/3)));

%	f11=1+a1*((3/2)*(fn1+theta1)-r*((3/2)*fn1+(5/2)*theta1-(4/3)));
%	f12=1+a2*((3/2)*(fn2+theta2)-r*((3/2)*fn2+(5/2)*theta2-(4/3)));


	f2a=1+a*(1+(3/2)*(fn+theta)-(r/2)*(3*fn+5*theta))+b*(3-4*r);
	f2a=f2a+(a/2)*(a+3*b)*(3-4*r)*(fn+theta-r*(fn-theta+2*theta^2));

%	f21=1+a1*(1+(3/2)*(fn1+theta1)-(r/2)*(3*fn1+5*theta1))+b1*(3-4*r);
%	f21=f21+(a1/2)*(a1+3*b1)*(3-4*r)*(fn1+theta1-r*(fn1-theta1+2*theta1^2));
%	f22=1+a2*(1+(3/2)*(fn2+theta2)-(r/2)*(3*fn2+5*theta2))+b2*(3-4*r);
%	f22=f22+(a2/2)*(a2+3*b2)*(3-4*r)*(fn2+theta2-r*(fn2-theta2+2*theta2^2));

	f3a=1+a*(1-(fn+(3/2)*theta)+r*(fn+theta));

%	f31=1+a1*(1-(fn1+(3/2)*theta1)+r*(fn1+theta1));
%	f32=1+a2*(1-(fn2+(3/2)*theta2)+r*(fn2+theta2));

	f4a=1+(a/4)*(fn+3*theta-r*(fn-theta));

%	f41=1+(a1/4)*(fn1+3*theta1-r*(fn1-theta1));
%	f42=1+(a2/4)*(fn2+3*theta2-r*(fn2-theta2));

	f5a=a*(-fn+r*(fn+theta-(4/3))) + b*theta*(3-4*r);

%	f51=a1*(-fn1+r*(fn1+theta1-(4/3))) + b1*theta1*(3-4*r);
%	f52=a2*(-fn2+r*(fn2+theta2-(4/3))) + b2*theta2*(3-4*r);

	f6a=1+a*(1+fn-r*(fn+theta))+b*(1-theta)*(3-4*r);

%	f61=1+a1*(1+fn1-r*(fn1+theta1))+b1*(1-theta1)*(3-4*r);
%	f62=1+a2*(1+fn2-r*(fn2+theta2))+b2*(1-theta2)*(3-4*r);

	f7a=2+(a/4)*(3*fn+9*theta-r*(3*fn+5*theta)) + b*theta*(3-4*r);

%	f71=2+(a1/4)*(3*fn1+9*theta1-r*(3*fn1+5*theta1)) + b1*theta1*(3-4*r);
%	f72=2+(a2/4)*(3*fn2+9*theta2-r*(3*fn2+5*theta2)) + b2*theta2*(3-4*r);

	f8a=a*(1-2*r+(fn/2)*(r-1)+(theta/2)*(5*r-3))+b*(1-theta)*(3-4*r);

%	f81=a1*(1-2*r+(fn1/2)*(r-1)+(theta1/2)*(5*r-3))+b1*(1-theta1)*(3-4*r);
%	f82=a2*(1-2*r+(fn2/2)*(r-1)+(theta2/2)*(5*r-3))+b2*(1-theta2)*(3-4*r);

	f9a=a*((r-1)*fn-r*theta) + b*theta*(3-4*r);

%	f91=a1*((r-1)*fn1-r*theta1) + b1*theta1*(3-4*r);
%	f92=a2*((r-1)*fn2-r*theta2) + b2*theta2*(3-4*r);


	pa=3*f1a/f2a;
	qa=(2/f3a) + (1/f4a) +((f4a*f5a + f6a*f7a - f8a*f9a)/(f2a*f4a));
	pa=pa/3.; qa=qa/5.;

%	p1=3*f11/f21; p2=3*f12/f22;
%	q1=(2/f31) + (1/f41) +((f41*f51 + f61*f71 - f81*f91)/(f21*f41));
%	q2=(2/f32) + (1/f42) +((f42*f52 + f62*f72 - f82*f92)/(f22*f42));

%	p1=p1/3; p2=p2/3;
%	q1=q1/5; q2=q2/5;

%------------------------------------------------------------------------
%zeta=(mu/6)*(9*k+8*mu)/(k+2*mu);

krhs=(ka-k)*pa;
yprime(1)=krhs/(1-t);

murhs=(mua-mu)*qa;
yprime(2)=murhs/(1-t);
