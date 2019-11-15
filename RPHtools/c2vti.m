function [vp,vsh,vsv]=c2vti(cc,rho,zeta)
% [VP,VSH,VSV]=C2VTI(CC,RHO,ZETA)
% calculate velocities (VP,VSH,VSV) at angle ZETA from the symmetry axis of
% a TI medium from its elastic stiffnesses
%
% input:    CC=[c11,c33,c44,c66,c13], RHO, and ZETA (in degree)
%           ZETA may be a vector of angles.
% output:   VP,VSH,VSV at angle(s) ZETA

% Written by Frank Liu

c11=cc(:,1);
c33=cc(:,2);
c44=cc(:,3);
c66=cc(:,4);
c13=cc(:,5);

conv=pi/180;
zeta=conv*zeta;

s2=sin(zeta).*sin(zeta);
c2=cos(zeta).*cos(zeta);
s22=sin(2*zeta).*sin(2*zeta);

mm=((c11-c44).*s2-(c33-c44).*c2).^2+(c13+c44).^2.*s22;
mm=sqrt(mm);

mp2=(c11.*s2+c33.*c2+c44+mm)./2;
msv2=(c11.*s2+c33.*c2+c44-mm)./2;
msh2=c66.*s2+c44.*c2;

vp=sqrt(mp2./rho);
vsv=sqrt(msv2./rho);
vsh=sqrt(msh2./rho);

