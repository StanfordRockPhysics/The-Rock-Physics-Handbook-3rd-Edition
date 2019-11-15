function [vp1,freq,vp2,vs,q1inv,q2inv,qsinv]=biot(vpdry,vsdry,k0,mu0,ro0,rofl,kfl,nu,por,perm,a,alfa,d1,d2,opt)
%[VP1,FREQ,VP2,VS,Q1INV,Q2INV,QSINV]
%      =BIOT(VPDRY,VSDRY,K0,MU0,RO0,ROFL,KFL,NU,POR,PERM,A,ALFA,D1,D2,'PLOTOPT')
%
%	calculates the complete Biot velocity dispersion and attenuation
%	curves for all frequencies. 
% Inputs: VPDRY: P-wave velocity of dry porous rock
%	  VSDRY: S-wave velocity of dry porous rock
%	  K0, MU0, RO0:  Mineral Bulk modulus, shear modulus, and density
%	  ROFL,KFL,NU: Pore fluid density, bulk modulus, and viscosity
%	  POR, PERM: Porosity and absolute permeability of porous rock
%	  A: pore size parameter; usually about 1/6-1/7 of grain diameter
%	  ALFA: Tortuosity parameter (always > 1), usually between 1-3
%	  D1,D2: frequency range from 10^D1 to 10^D2
%	  'PLOTOPT': optional input for plotting. If nothing is given, plots of
%	  all six outputs versus frequency are displayed. The following options
%	  selects only one of the outputs for plotting.
%	  'PLOTOPT'='VP1' displays VP1 (fast P-wave) dispersion curve
%		   ='VP2' displays VP2 (slow P-wave) dispersion curve
%		   ='VS'  displays VS (shear wave) dispersion curve
%		   ='QP1' displays P1 attenuation curve
%		   ='QP2' displays P2 attenuation curve
%		   ='QS'  displays S attenuation curve
% Outputs: VP1: Fast P-wave velocities at all frequencies
%	   FREQ: Frequencies at which velocity and attenuation are calculated
%	   VP2: Slow P-wave velocities at all frequencies
%	   VS:  S-wave velocities
%	   Q1INV, Q2INV, QSINV: Fast and slow P-wave and S-wave attenuations
%
%	   See also BIOTHF, BIOTHFB

% Written by T. Mukerji, 1996

if nargin<15, opt='all6'; end; 
ro=(1-por)*ro0+por*rofl; rodry=(1-por)*ro0; omc=por*nu/(rofl*perm);
mudry=rodry.*vsdry.^2; kdry=rodry.*vpdry.^2-(4/3)*mudry;
d=k0*(1+por*(k0/kfl-1)); h=kdry+(4/3)*mudry+(k0-kdry).^2./(d-kdry);
c=k0*(k0-kdry)./(d-kdry); m=k0.^2./(d-kdry);
freq=logspace(d1,d2,100); om=2*pi*freq; zeta=sqrt((rofl*a^2/nu)*om);
t=zeros(size(zeta)); idx=[zeta <= 1e3]; idxl=[zeta <= 1e-1];
t(idx)=exp(i*3*pi/4)*bessel(1,exp(-i*pi/4)*zeta(idx))./...
        bessel(0,exp(-i*pi/4)*zeta(idx));
t(~idx)=((1+i)/sqrt(2))*ones(sum(~idx),1);
f=0.25*zeta.*t./(1+i*2*t./zeta); f(idxl)=ones(sum(idxl),1);
qf=alfa*rofl/por - (i*nu/perm)*f./om;
%%%%%%%%%%%%test of qf
%qf=alfa*rofl/por + (i*nu/perm)*f./om; this equation is wrong; sign is -
%%%%%%%%%%%
b0=rofl^2-ro*qf; b1=h*qf+m*ro-2*c*rofl; b2=c^2-m*h;
sl1sqr=(-b1+sqrt(b1.^2-4*b2*b0))./(2*b2);
sl2sqr=(-b1-sqrt(b1.^2-4*b2*b0))./(2*b2);
slssqr=(ro*qf-rofl^2)./(mudry*qf);
vp1=1./real(sqrt(sl1sqr)); vp2=1./real(sqrt(sl2sqr)); vs=1./real(sqrt(slssqr));
q1inv=imag(1./sl1sqr)./real(1./sl1sqr); q2inv=imag(1./sl2sqr)./real(1./sl2sqr);
qsinv=imag(1./slssqr)./real(1./slssqr); 
if strcmp(opt,'all6')
subplot(2,3,1),semilogx(freq,vp1),title('Vp1');
subplot(2,3,2),semilogx(freq,vp2),title('Vp2');
subplot(2,3,3),semilogx(freq,vs),title('Vs');
subplot(2,3,4),semilogx(freq,q1inv),title('1/Qp1');
subplot(2,3,5),semilogx(freq,q2inv),title('1/Qp2');
subplot(2,3,6),semilogx(freq,qsinv),title('1/Qs');
subplot(2,3,1);
elseif strcmp(opt,'vp1')
semilogx(freq,vp1), title('Vp1');
elseif strcmp(opt,'vp2')
semilogx(freq,vp2), title('Vp2');
elseif strcmp(opt,'vs')
semilogx(freq,vs), title('Vs');
elseif strcmp(opt,'qp1')
semilogx(freq,q1inv), title('1/Qp1');
elseif strcmp(opt,'qp2')
semilogx(freq,q2inv), title('1/Qp2');
elseif strcmp(opt,'qs')
semilogx(freq,qsinv), title('1/Qs');
else
end
