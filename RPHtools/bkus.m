function [vv,cc,rho]=bkus(f,r,vp,vs)
% [VV,CC,RHO]=BKUS(F,R,VP,VS)
% Backus average for the effective anisotropic (TI) stiffness, CC, and 
% velocity, VV, of layered medium made of thin isotropic layers.
%
% input:   F - volume fraction (layer thickness/total thickness)
%          R - density
%          VP - p velocity
%          VS - s velocity
% output:  VV=[vp33 vp13 vp11 vs33 vsh11];  CC=[c11 c33 c44 c66 c13]
%

% Written by Xingzhou 'Frank' Liu

if abs(sum(f)-1)>.0001,
disp('The sum of fractions <> 1 !')
pause
end
mu=r.*vs.*vs;
lam=r.*vp.*vp-2*mu;

x=sum(f.*mu.*(lam+mu)./(lam+2*mu));
y=sum(f.*mu.*lam./(lam+2*mu));
z=sum(f.*lam./(lam+2*mu));
u=sum(f./(lam+2*mu));
v=sum(f./mu);
w=sum(f.*mu);

c11=4*x+z*z/u;
c12=2*y+z*z/u;
c33=1/u;
c13=z/u;
c44=1/v;
c66=w;
c661=(c11-c12)/2;

if abs(c661-c66)>.1
disp('There is a problem with the calculation')
disp('c66 is not equal to (c11-c12)/2')
pause
end

cc=[c11 c33 c44 c66 c13];

rho=sum(f.*r);

vp11=sqrt(c11/rho);
vp33=sqrt(c33/rho);
vsh11=sqrt(c66/rho);
vs33=sqrt(c44/rho);

vp13=c2vti(cc,rho,45);

vv=[vp33 vp13 vp11 vs33 vsh11];

