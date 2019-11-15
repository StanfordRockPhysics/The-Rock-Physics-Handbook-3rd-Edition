function [Ctih,den]=hudson(ec,ar,Kfl,rhofl,K,G,rho,ax)
%function [Ctih,den]=hudson(ec,ar,Kfl,rhofl,K,G,rho,ax)
% hudson - calculate the anisotropic elastic parameters for cracked rock 
%           using Hudson's first order weak inclusion theory valid for small
%           crack density and aspect ratios. Assumes a 
%           single crack set with all normals aligned along 1 or 3-axis.
%           For dry cracks use fluid bulk modulus 0
%
% input and output parameters:
% den - bulk density of cracked rock with fluid
% Ctih - C6x6 matrix with symmetry axis along aligned crack normal
%       tih stands for transversely isotropic with a horizontal symmetry axis
%        If inputs are nx1 vectors, Ctih is a 6x6xn array.
% ec - crack density
% ar - aspect ratio of cracks
% Kfl, K: bulk modulus of fluid in cracks, and bulk modulus of isotropic matrix
% rhofl:  density of fluid in cracks
% G: shear modulus of isotropic matrix
% rho: density of matrix
%
% optional input parameter: ax - defines axis of symmetry; ax =1, for crack
% normals aligned along 1-axis (default); ax= 3, for crack normals along
% 3-axis.
%
% See also ECHENG, MCHUDSON, CTI2V

% References:
% Hudson (1980, 1981, 1990)
% L. Thomsen, 1986, Weak elastic anisotropy, Geophysics Vol51 Oct 1986
%
% Author:	Li Teng
% Date:		7/20/98
% axis option added 2/5/98 T. Mukerji

if nargin<8, ax=1; end;

lam=K-2/3.*G;
mu=G;
kapa=Kfl.*(lam+2.*mu)./(pi.*ar.*mu.*(lam+mu));
u3=4/3.*(lam+2.*mu)./((lam+mu).*(1+kapa));
u1=16/3.*(lam+2.*mu)./(3.*lam+4.*mu);
c11=lam+2.*mu-lam.^2.*ec.*u3./mu;
c13=lam-lam.*(lam+2.*mu).*ec.*u3./mu;
c33=lam+2.*mu-(lam+2.*mu).^2.*ec.*u3./mu;
c44=mu-mu.*ec.*u1;
c66=mu;

%Ctih=zeros(6,6);

if ax==1
Ctih(1,1,:)=c33;
Ctih(2,2,:)=c11;
Ctih(3,3,:)=c11;
Ctih(1,3,:)=c13;
Ctih(3,1,:)=c13;
Ctih(1,2,:)=c13;
Ctih(2,1,:)=c13;
Ctih(2,3,:)=c11-2.*c66;
Ctih(3,2,:)=Ctih(2,3,:);
Ctih(4,4,:)=c66;
Ctih(5,5,:)=c44;
Ctih(6,6,:)=c44;
end;
if ax==3
Ctih(1,1,:)=c11;
Ctih(2,2,:)=c11;
Ctih(3,3,:)=c33;
Ctih(1,3,:)=c13;
Ctih(3,1,:)=c13;
Ctih(1,2,:)=c11-2.*c66;
Ctih(2,1,:)=Ctih(1,2,:);
Ctih(2,3,:)=c13;
Ctih(3,2,:)=Ctih(2,3,:);
Ctih(4,4,:)=c44;
Ctih(5,5,:)=c44;
Ctih(6,6,:)=c66;
end;


phi=(4*pi/3).*ar.*ec;
den=(1-phi).*rho + phi.*rhofl;

Vp0=sqrt(c33./den);
Vs0=sqrt(c44./den);
e=(c11-c33)./(2.*c33);
g=(c66-c44)./(2.*c44);
d=((c13+c44).^2-(c33-c44).^2)./(2.*c33.*(c33-c44));

