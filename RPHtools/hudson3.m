function   [C,den,ex,gx,dx,ey,gy,dy,gxy]=hudson3(ec,ar,Kfl,rofl,K,G,ro)
%function   [C,den,ex,dx,ey,dy,gxy]=hudson3(ec,ar,Kfl,rofl,K,G,ro)
% hudson3 - calculates the anisotropic elastic parameters for cracked rock 
%           using Hudson's first order weak inclusion theory valid for small
%           crack density and aspect ratios. Assumes 3 sets of cracks
%           with normals aligned with the 3 principal axes. 
%         
% input and output parameters:
%
%Outputs:
% C           - C6x6 stiffness matrix
% den         - density of cracked rock
% ex, gx, dx, ey, gy, dy, gxy    - epsilon,  gamma, delta, equivalent Thomsen's weakly anisotropic parameters
%
%Inputs:
% ec=[ec1,ec2,ec3]  crack densities with normal || to axes 1,2,3 respectively
% ar=[ar1,ar2,ar3]  aspect ratio of the 3 sets of cracks
% Kfl, K: bulk modulus of fluid in cracks, and bulk modulus of isotropic host
% G: shear modulus of isotropic host
% rofl, ro: bulk density of fluid and isotropic host rock
%
% See also HUDSON1, HUDSONCONE, HUDSONF,ECHENG
% References:
% Hudson (1980, 1981, 1990)
% L. Thomsen, 1986, Weak elastic anisotropy, Geophysics Vol 51 Oct 1986
%

% Written by Diana Sava
% 06/27/00


lam=K-2/3.*G;
mu=G;
%ec=[ec1 ec2 ec3];
kapa=Kfl.*(lam+2.*mu)./(pi.*ar.*mu.*(lam+mu));
u3=4/3*(lam+2*mu)./((lam+mu).*(1+kapa));
u1=16/3.*(lam+2*mu)./(3.*lam+4.*mu);

c11cor=-lam.^2.*ec.*u3./mu;
c13cor=-lam.*(lam+2*mu).*ec.*u3./mu;
c33cor=-(lam+2*mu).^2.*ec.*u3./mu;
c44cor=-mu.*ec.*u1;
c66cor=0*ec;
c12cor=c11cor;


c11=lam+2*mu+c33cor(1)+c11cor(2)+c11cor(3);
c12=lam     +c13cor(1)+c13cor(2)+c12cor(3);
c13=lam     +c13cor(1)+c12cor(2)+c13cor(3);
c22=lam+2*mu+c11cor(1)+c33cor(2)+c11cor(3);
c23=lam     +c12cor(1)+c13cor(2)+c13cor(3);
c33=lam+2*mu+c11cor(1)+c11cor(2)+c33cor(3);
c44=mu      +c44cor(2)+c44cor(3);
c55=mu      +c44cor(1)+c44cor(3);
c66=mu      +c44cor(1)+c44cor(2);

C=zeros(6,6);

C(1,1)=c11;
C(2,2)=c22;
C(3,3)=c33;
C(1,3)=c13;
C(3,1)=C(1,3);
C(1,2)=c12;
C(2,1)=C(1,2);
C(2,3)=c23;
C(3,2)=C(2,3);
C(4,4)=c44;
C(5,5)=c55;
C(6,6)=c66;

phic=4*pi/3*sum(ar.*ec);
den=(1-phic).*ro+phic.*rofl;

ey=(c11-c33)./(2*c33);   % or e2  in literature (defined in the symmetry plane with normal in y direction);
dy=((c13+c55).^2-(c33-c55).^2)./(2*c33.*(c33-c55));
gy=(c66-c44)./2./c44;

ex=(c22-c33)./(2*c33);  % or e1 in literature (defined in the symmetry plane with normal in x direction);
dx=((c23+c44).^2-(c33-c44).^2)./(2*c33.*(c33-c44));
gx=(c66-c55)./2./c55;

gxy=(c44-c55)./(2*c55);

% d3 defined in the horizontal plane
d3=((c12+c66).^2-(c11-c66).^2)./2./c11./(c11-c66);
