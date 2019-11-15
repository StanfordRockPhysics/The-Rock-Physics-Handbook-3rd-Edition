function [vp,vs]=ku2v(k,u,rho)
% [VP,VS]=KU2V(K,U,RHO)
% get velocities from moduli and density for isotropic media
% input:    K,U,RHO bulk and shear moduli, and density
% output:   VP,VS 

% Written by Frank Liu

c=4/3;
vp=sqrt((k+c*u)./rho);
vs=sqrt(u./rho);

if nargout==0, [vp, vs], end;
