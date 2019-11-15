function [vp,k,atn,f,kinf,klf]=patchw(kdry,mudry,k0,mu0,ro0,phi,perm,fl,sg,a,f,plt)
%PATCHW - Velocity dispersion and attenuation for patchy saturation using
%	  White's model (incorporating Dutta-Ode correction).
%
%[VP,K,ATN,F,KINF,KLF]=PATCHW(KDRY,MUDRY,K0,MU0,RO0,PHI,PERM,FL,SG,A,F,PLT)
%
%KDRY, MUDRY:	Dry rock bulk and shear moduli
%K0, MU0, RO0:	Mineral bulk and shear moduli, and density
%PHI, PERM:	Porosity and absolute permeability
%FL, A, SG:	FL=[KF1 KF2;ROF1 ROF2;NU1 NU2] Bulk moduli (K), density (RO),
%		and viscosity (NU) of fluids 1 and 2. In the White model
%		fluid 1 occupies a central sphere of radius A, surrounded 
%		by a shell saturated with fluid 2. SG is overall saturation of
%		fluid 1.
%F:		Frequencies at which to calculate velocity and attenuation
%
%VP, K, ATN:	Velocity (P-wave), complex bulk modulus, and attenuation at
%		frequencies F.
%KINF, KLF:	Bulk moduli at high and low frequency limits.
%
%With no output arguments PATCHW plots VP, K, and ATN, versus frequency.
%PLT:		Optional plot parameter. PLOT=1 and no output arguments
%		plots only VP versus frequency. PLOT=2 and no
%		output arguments plots only the bulk modulus. 

%Written by T. Mukerji

kf1=fl(1,1);rof1=fl(2,1);nu1=fl(3,1);
kf2=fl(1,2);rof2=fl(2,2);nu2=fl(3,2);

i=sqrt(-1); om=2*pi*f; %sg=(a/b)^3; 
b=a./(sg.^(1/3));

k1=gassmnk(kdry,0.0,kf1,k0,phi); k2=gassmnk(kdry,0.0,kf2,k0,phi);
mu1=mudry; mu2=mudry;
m1=k1+4*mu1/3; m2=k2+4*mu2/3;

r1=(k1-kdry).*(3*k2+4*mu2);
t1=(1-kdry./k0).*(k2.*(3*k1+4*mu2)+4*mu2.*(k1-k2).*sg);
r1=r1./t1;

r2=(k2-kdry).*(3*k1+4*mu1);
t2=(1-kdry./k0).*(k2.*(3*k1+4*mu2)+4*mu2.*(k1-k2).*sg);
r2=r2./t2;

ka1=phi./kf1 + (1-phi)./k0 - kdry./k0.^2; ka1=1./ka1;
ka2=phi./kf2 + (1-phi)./k0 - kdry./k0.^2; ka2=1./ka2;
q1=(1-kdry./k0).*ka1./k1; q2=(1-kdry./k0).*ka2./k2;
ke1=kf1.*(1-k1./k0).*(1-kdry./k0)./(phi.*k1.*(1-kf1./k0)); ke1=(1-ke1).*ka1;
ke2=kf2.*(1-k2./k0).*(1-kdry./k0)./(phi.*k2.*(1-kf2./k0)); ke2=(1-ke2).*ka2;
alpha1=sqrt(i*om.*nu1./(perm.*ke1)); alpha2=sqrt(i*om.*nu2./(perm.*ke2));
z1=(1-exp(-2.*a.*alpha1)).*nu1./perm;
t1=(alpha1.*a-1) + (alpha1.*a+1).*exp(-2.*a.*alpha1);
z1=z1./t1;
z2=(alpha2.*b+1) + (alpha2.*b-1).*exp(2.*(b-a).*alpha2);
z2=z2.*nu2./perm;
t1=(alpha2.*b+1).*(alpha2.*a-1);
t2=(alpha2.*b-1).*(alpha2.*a+1).*exp(2*alpha2.*(b-a));
z2=z2./(t1-t2);
z2=-z2;

w=3.*a.*(r1-r2).*(q2-q1)./(i*b.^3.*om.*(z1+z2));

kinf=k2.*(3*k1+4*mu2)+4*mu2.*(k1-k2).*sg;
kinf=kinf./((3*k1+4*mu2)-3*(k1-k2).*sg)

klf=k2.*(k1-kdry) + sg.*kdry.*(k2-k1); klf=klf./((k1-kdry)+sg.*(k2-k1))

k=kinf./(1-kinf.*w);
kapp=kinf./(1+3*a*r2*q2*kinf./(i*b^3*om.*z2));

ro=(1-phi)*ro0+phi*sg*rof1+phi*(1-sg)*rof2;
m=k+4*mudry/3;
theta=atan(imag(m)./real(m));
vp=sqrt(abs(m)./ro)./cos(theta/2);
atn=om.*tan(theta/2)./vp;
vpinf=sqrt((kinf+4*mudry/3)/ro), vplf=sqrt((klf+4*mudry/3)/ro)

if nargout==0
 if nargin<12
subplot(2,2,1)
semilogx(f,real(k),min(f),klf,'*w',max(f),kinf,'*w');
title('Bulk Modulus (Re(K))')
subplot(2,2,2)
semilogx(f,imag(k))
title('Bulk Modulus (Im(K))')
subplot(2,2,3)
semilogx(f,vp,min(f),vplf,'*w',max(f),vpinf,'*w');
title('Velocity')
subplot(2,2,4)
semilogx(f,atn)
title('Attenuation')
 elseif plt==1
semilogx(f,vp,min(f),vplf,'*w',max(f),vpinf,'*w');
 elseif plt==2
semilogx(f,real(k),min(f),klf,'*w',max(f),kinf,'*w');
end; end;

