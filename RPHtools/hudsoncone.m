function [Vp0,Vs0,e,g,d,C]=hudsoncone(ec,ar,Kfl,K,G,den,theta,ax)
%function [Vp0,Vs0,e,g,d,Ctih]=hudsoncone(ec,ar,Kfl,K,G,den,theta,ax)
% hudsoncone - calculates the anisotropic elastic parameters for cracked rock 
%              using Hudson's first order weak inclusion theory valid for small
%              crack density and aspect ratios.
%            - cracks have the normals randomly distributed at a fixed
%              angle from the TI symmetry axis;
%             
% input and output parameters:
% Vp0, Vs0 - P and S velocities of the cracked rock along the symmetry axis 
% e, g, d - epsilon,  gamma, delta, Thomsen's weakly anisotropic parameter
% C  - C6x6 compliance matrix with symmetry axis along 1-axis (default)
%      
% ec - crack density
% ar - aspect ratio of cracks
% Kfl, K: bulk modulus of fluid in cracks, and bulk modulus of isotropic matrix
% G: shear modulus of isotropic matrix
% den: density of cracked rock
% theta: the constant angle between the crack normals and the symmetry axis of TI medium
%
%optional input parameter: ax - defines axis of symmetry; ax =1, 
% (default); ax= 3, 
%
% See also HUDSON1, ECHENG
%
% References:
% Hudson (1980, 1981, 1990)
% L. Thomsen, 1986, Weak elastic anisotropy, Geophysics Vol51 Oct 1986

% Written by    Diana Sava
% Date:         03/16/00

if nargin<7, ax=1; 
end;

lam=K-2/3.*G;
mu=G;
t=theta;
kapa = Kfl.*(lam+2.*mu)./(pi.*ar.*mu.*(lam+mu));
u3 = 4/3.*(lam+2.*mu)./((lam+mu).*(1+kapa));
u1 = 16/3.*(lam+2.*mu)./(3.*lam+4.*mu);

c11cor = -ec./mu./2.*(u3.*(2*lam.^2+4*lam.*mu.*sin(t).^2+3*mu.^2.*sin(t).^4)+u1.*mu.^2.*sin(t).^2.*(4-3*sin(t).^2));

c33cor = -ec./mu.*(u3.*(lam+2*mu.*cos(t).^2).^2+u1.*mu.^2.*4.*cos(t).^2.*sin(t).^2);

c12cor = -ec./mu/2.*(u3.*(2*lam.^2+4*lam.*mu.*sin(t).^2+mu.^2.*sin(t).^4)-u1.*mu.^2.*sin(t).^4);

c13cor = -ec./mu.*(u3.*(lam+mu.*sin(t).^2).*(lam+2*mu.*cos(t).^2)-u1.*mu.^2*2.*sin(t).^2.*cos(t).^2);

c44cor = -ec/2.*mu.*(u3*4.*sin(t).^2.*cos(t).^2+u1.*(sin(t).^2+2*cos(t).^2-4*sin(t).^2.*cos(t).^2));

c66cor = -ec/2.*mu.*(u3*4*sin(t).^4+u1.*sin(t).^2.*(2-sin(t).^2));

c22cor = c11cor;
c23cor = c13cor;
c55cor = c44cor;

c11 = lam + 2*mu + c11cor;
c12 = lam        + c12cor;
c13 = lam        + c13cor;
c33 = lam + 2*mu + c33cor;
c44 = mu         + c44cor;
c66 = mu         + c66cor;


C=zeros(6,6);

if ax==1
C(1,1)=c33;
C(2,2)=c11;
C(3,3)=c11;
C(1,3)=c13;
C(3,1)=c13;
C(1,2)=c13;
C(2,1)=c13;
C(2,3)=c11-2.*c66;
C(3,2)=C(2,3);
C(4,4)=c66;
C(5,5)=c44;
C(6,6)=c44;
end;
if ax==3
C(1,1)=c11;
C(2,2)=c11;
C(3,3)=c33;
C(1,3)=c13;
C(3,1)=c13;
C(1,2)=c11-2.*c66;
C(2,1)=C(1,2);
C(2,3)=c13;
C(3,2)=C(2,3);
C(4,4)=c44;
C(5,5)=c44;
C(6,6)=c66;
end;


Vp0=sqrt(c33./den);
Vs0=sqrt(c44./den);
e=(c11-c33)./(2.*c33);
g=(c66-c44)./(2.*c44);
d=((c13+c44).^2-(c33-c44).^2)./(2.*c33.*(c33-c44));
