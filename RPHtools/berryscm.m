function [kbr,mubr]=berryscm(k,mu,asp,x)
%BERRYSCM - Effective elastic moduli for multi-component composite
% using Berryman's Self-Consistent (Coherent Potential Approximation) method. 
%
%[KBR,MUBR]=BERRYSCM(K,MU,ASP,X)
%	K,MU:          Bulk and shear moduli of the N constituent
%		       phases (K, MU, vectors of length N)
%	ASP:           Aspect ratio for the inclusions of the N phases
%			< 1 for oblate spheroids; >1 for prolate spheroids.
%       X:              Fraction of each phase. Sum(X) should be 1.
%	KBR,MUBR:  	Effective bulk and shear moduli 
%
% See also BERRYSC

% Written by T. Mukerji, 1997
%
%******* Berryman SC oblate and prolate spheroids *****************

kbr=[]; mubr=[]; %por=[];

k=k(:); mu=mu(:); asp=asp(:); x=x(:);
indx=find(asp==1); asp(indx)=0.99*ones(size(indx));
theta=zeros(size(asp)); fn=zeros(size(asp));  
%

obdx=find(asp<1);
theta(obdx)=(asp(obdx)./((1-asp(obdx).^2).^(3/2))).*...
             (acos(asp(obdx)) -asp(obdx).*sqrt(1-asp(obdx).^2));
fn(obdx)=(asp(obdx).^2./(1-asp(obdx).^2)).*(3.*theta(obdx) -2);

%if asp2 < 1.
%theta2=(asp2/((1-asp2^2)^(3/2)))*(acos(asp2) -asp2*sqrt(1-asp2^2));
%fn2=(asp2^2/(1-asp2^2))*(3*theta2 -2);
%end

prdx=find(asp>1);
theta(prdx)=(asp(prdx)./((asp(prdx).^2-1).^(3/2))).*...
             (asp(prdx).*sqrt(asp(prdx).^2-1)-acosh(asp(prdx)));
fn(prdx)=(asp(prdx).^2./(asp(prdx).^2-1)).*(2-3.*theta(prdx));

%if asp2 > 1.
%theta2=(asp2/((asp2^2-1)^(3/2)))*(asp2*sqrt(asp2^2-1)-acosh(asp2));
%fn2=(asp2^2/(asp2^2-1))*(2-3*theta2);
%end


%for x1=0:0.01:1
% x2=1-x1;
% for x2=0.:0.001:0.6
%    x1=1-x2;

%ksc= x1*k1 + x2*k2;
ksc= sum(k.*x);
musc= sum(mu.*x);
knew= 0.;
munew= 0.;
%tol=real(ksc)/50000.0;
tol=1e-6*k(1);
del=abs(ksc-knew);
niter=0;

while( (del > abs(tol)) & (niter<3000) )
	nusc=(3*ksc-2*musc)/(2*(3*ksc+musc));
	a=mu./musc -1; 
	b=(1/3)*(k./ksc -mu./musc); 
	r=(1-2*nusc)/(2*(1-nusc));

	f1=1+a.*((3/2).*(fn+theta)-r.*((3/2).*fn+(5/2).*theta-(4/3)));
%	f12=1+a2*((3/2)*(fn2+theta2)-r*((3/2)*fn2+(5/2)*theta2-(4/3)));

	f2=1+a.*(1+(3/2).*(fn+theta)-(r/2).*(3.*fn+5.*theta))+b.*(3-4*r);
	f2=f2+(a/2).*(a+3.*b).*(3-4.*r).*(fn+theta-r.*(fn-theta+2.*theta.^2));
%	f22=1+a2*(1+(3/2)*(fn2+theta2)-(r/2)*(3*fn2+5*theta2))+b2*(3-4*r);
%	f22=f22+(a2/2)*(a2+3*b2)*(3-4*r)*(fn2+theta2-r*(fn2-theta2+2*theta2^2));

	f3=1+a.*(1-(fn+(3/2).*theta)+r.*(fn+theta));
%	f32=1+a2*(1-(fn2+(3/2)*theta2)+r*(fn2+theta2));

	f4=1+(a./4).*(fn+3.*theta-r.*(fn-theta));
%	f42=1+(a2/4)*(fn2+3*theta2-r*(fn2-theta2));

	f5=a.*(-fn+r.*(fn+theta-(4/3))) + b.*theta.*(3-4*r);
%	f52=a2*(-fn2+r*(fn2+theta2-(4/3))) + b2*theta2*(3-4*r);

	f6=1+a.*(1+fn-r.*(fn+theta))+b.*(1-theta).*(3-4.*r);
%	f62=1+a2*(1+fn2-r*(fn2+theta2))+b2*(1-theta2)*(3-4*r);

	f7=2+(a./4).*(3.*fn+9.*theta-r.*(3.*fn+5.*theta)) + b.*theta.*(3-4.*r);
%	f72=2+(a2/4)*(3*fn2+9*theta2-r*(3*fn2+5*theta2)) + b2*theta2*(3-4*r);

      f8=a.*(1-2.*r+(fn./2).*(r-1)+(theta./2).*(5.*r-3))+b.*(1-theta).*(3-4.*r);
%	f82=a2*(1-2*r+(fn2/2)*(r-1)+(theta2/2)*(5*r-3))+b2*(1-theta2)*(3-4*r);

	f9=a.*((r-1).*fn-r.*theta) + b.*theta.*(3-4.*r);
%	f92=a2*((r-1)*fn2-r*theta2) + b2*theta2*(3-4*r);


	p=3*f1./f2; %p2=3*f12/f22;
	q=(2./f3) + (1./f4) +((f4.*f5 + f6.*f7 - f8.*f9)./(f2.*f4));
%	q2=(2/f32) + (1/f42) +((f42*f52 + f62*f72 - f82*f92)/(f22*f42));

	p=p./3; %p2=p2/3;
	q=q./5; %q2=q2/5;

%------------------------------------------------------------------------

	knew= sum(x.*k.*p)/sum(x.*p);
	munew= sum(x.*mu.*q)/sum(x.*q);
	
	del=abs(ksc-knew);
	ksc=knew;
	musc=munew;
	niter=niter+1;
end		
kbr=ksc; mubr=musc;
% x1, niter  
%	vpbr=sqrt((ksc+ (4*musc/3))/ro);
%	vsbr=sqrt(musc/ro);
%	kbr=[kbr,ksc]; mubr=[mubr, musc]; por=[por, 1-x1];
%end
%end
%if nargout==0
%plot(por,real(kbr),'-g',por,real(mubr),'-c','linewidth',1)
%end
%plot(por,real(mubr))
